func_adj <- function(data_chart, data_emr, data_phone, data_all3, raking=FALSE, formula) {
  
  
  data_emr$fptype[grep('UNK Implant', data_emr$fptype)] <- 'UNK implant'
  if (any(data_emr$fptype %in% 'Eto implant') & any(data_emr$fptype %in% 'Depo')){
    data_emr$fptype <- as.character(data_emr$fptype)
    if (!is.null(grep('Implant', data_emr$fptype)))
      data_emr$fptype[grep('Implant', data_emr$fptype)] <- 'Implant'
    data_emr$fptype[grep('implant', data_emr$fptype)] <- 'Implant'
  }
  if (any(data_chart$fptype %in% 'Eto implant') & any(data_emr$fptype %in% 'Depo')){
    data_chart$fptype <- as.character(data_chart$fptype)
    data_chart$fptype[grep('implant', data_chart$fptype)] <- 'Implant'
  }
  if (any(data_phone$fptype %in% 'Eto implant') & any(data_emr$fptype %in% 'Depo')){
    data_phone$fptype <- as.character(data_phone$fptype)
    data_phone$fptype[grep('implant', data_phone$fptype)] <- 'Implant'
  }
  
  if (any(data_chart$fptype %in% 'Other MEFP')){
    data_chart$fptype <- as.character(data_chart$fptype)
    data_chart$fptype[grep('Other MEFP', data_chart$fptype)] <- 'MEFP'
  }
  if (any(data_phone$fptype %in% 'Other MEFP')){
    data_phone$fptype <- as.character(data_phone$fptype)
    data_phone$fptype[grep('Other MEFP', data_phone$fptype)] <- 'MEFP'
  }
  if (any(data_chart$fptype %in% 'OCP')){
    data_chart$fptype <- as.character(data_chart$fptype)
    data_chart$fptype[grep('OCP', data_chart$fptype)] <- 'Oral'
  }
  if (any(data_phone$fptype %in% 'OCP')){
    data_phone$fptype <- as.character(data_phone$fptype)
    data_phone$fptype[grep('OCP', data_phone$fptype)] <- 'Oral'
  }
  
  data_emr$fptype   <- as.factor(data_emr$fptype)
  data_chart$fptype <- as.factor(data_chart$fptype)
  data_phone$fptype <- as.factor(data_phone$fptype)
  
  if (any(data_emr$fptype %in% 'Depo')) {
    data_emr$arttype   <- relevel(data_emr$arttype,   ref='Nevirapine')
    data_chart$arttype <- relevel(data_chart$arttype, ref='Nevirapine')
    data_phone$arttype <- relevel(data_phone$arttype, ref='Nevirapine')
    data_emr$fptype   <- relevel(data_emr$fptype,   ref='Implant')
    data_chart$fptype <- relevel(data_chart$fptype, ref='Implant')
    data_phone$fptype <- relevel(data_phone$fptype, ref='Implant')
  } else {
    data_emr$fptype   <- relevel(data_emr$fptype,   ref='Eto implant')
    data_chart$fptype <- relevel(data_chart$fptype, ref='Eto implant')
    data_phone$fptype <- relevel(data_phone$fptype, ref='Eto implant')
  }
  data_phone$arttype <- factor(data_phone$arttype, levels=levels(data_chart$arttype))  
  
  ## EMR
  if_design <- svydesign(id = ~id, weights = ~ NULL, data = data_emr)
  res_full  <- svyglm(formula, design=if_design, family=poisson(link=log))
  
  ## Complete case chart review
  if_design <- svydesign(id = ~id, weights = ~ NULL, data = data_chart)
  res_2ph   <- svyglm(formula, design=if_design, family=poisson(link=log))
  
  ## Complete case phone interview
  if_design <- svydesign(id = ~ id, weights = ~ NULL, data = data_phone)
  res_3ph   <- svyglm(formula, design=if_design, family=poisson(link=log))
  
  ## Get the first probability, from phase-1 to phase-2
  if (!exists('Prob_keep')){
    Prob <<- cal_prob(data_all3, data_samp, data_chart, data_emr)
    Prob <- Prob_keep
  }
  Prob <- Prob_keep
  data_chart <- data_chart %>% left_join(Prob, by='id')
  
  ## IPW  
  if_design  <- svydesign(id = ~id, weights = ~ wgt, data = data_chart)
  resIPW_2ph <- svyglm(formula, design=if_design, family=poisson(link=log))
  
  ## 3-phase
  data_fullph2 <- data_emr %>% filter(id %in% data_chart$id) %>%
    mutate(cat = ifelse(fptype == 'Implant' & pregnant == 1, 4,
                        ifelse(fptype == 'Implant' & pregnant == 0, 3,
                               ifelse(fptype == 'Depo' & pregnant == 1, 2, 1))))
  
  data_ph3 <- data_fullph2 %>% 
    dplyr::group_by(id) %>%
    mutate(R = ifelse(id %in% data_phone$id, 1, 0))
  
  data_ph3b = aggregate(data_ph3$cat, list(data_ph3$id), max)
  colnames(data_ph3b) <- c('id', 'comb_ph3')
  data_ph3 <- data_ph3 %>% left_join(data_ph3b, by='id') %>% group_by(id) %>% slice(1)
  
  ## Get second probability, from phase-2 to phase-3
  fit_ps3 <- glm(R ~ comb_ph3, data=data_ph3, family=binomial)
  data_ph3$ps3 <- predict(fit_ps3, newdata=data_ph3, type='response')
  data_phase3 <- data_ph3 %>% dplyr::select(id, ps3, comb_ph3)
  data_phone2 <- data_phone %>%
    left_join(data_phase3, by='id') %>% filter(!is.na(ps3)) %>%
    left_join(Prob, by='id') %>% ungroup() %>% droplevels() %>% 
    mutate(ps = Probs, ps_prod = ps*ps3, wgt = 1/ps_prod, wgt2 = 1/ps3,
           fudiffyears2 = as.numeric(as.Date(end_date) - as.Date(start_date))/365.25)
  
  ## IPW
  if_design  <- svydesign(id = ~ id, weights = ~ wgt, data = data_phone2)
  resIPW_3ph <- svyglm(formula, design=if_design, family=poisson(link=log))
  
  ## Raking 
  if (raking) {
    data_emr    <- data_emr %>% filter(arttype %in% c("Nevirapine", "EFV", "Not on ART", "PI")) %>% droplevels()
    data_chart  <- data_chart %>% filter(arttype %in% c("Nevirapine", "EFV", "Not on ART", "PI")) %>% droplevels()
    data_phone2 <- data_phone2 %>% filter(arttype %in% c("Nevirapine", "EFV", "Not on ART", "PI")) %>% droplevels()
    datasets_if <- rak_func(data_emr, data_chart, data_phone2, formula)

    ## 2 and 3 phase: calibrate the product of weights
    resCAL_2ph <- datasets_if$rescal2ph
    resCAL_3ph <- datasets_if$rescal3ph

    ## 2 and 3 phase: calibrate each weight separately
    raking_3ph  <- datasets_if$raking3phb
    varsRak_3ph <- datasets_if$varsRak3phb
  } else (
    resCAL_2ph <- resCAL_3ph <- raking_3phb <- varsRak_3phb <- NULL
  )
  
  return(list(res_full=res_full, res_2ph=res_2ph, res_3ph=res_3ph,
              resIPW_2ph=resIPW_2ph, resIPW_3ph=resIPW_3ph, resCAL_2ph=resCAL_2ph, resCAL_3ph=resCAL_3ph,
              raking_3phb=raking_3ph, varsRak_3phb=varsRak_3ph))
}

probs <- function(p, cat) {
  zeros <- !((1:32) %in% cat)
  p[zeros] <- 0
  p <- p[p > 0]
  p1 <- (1-p)[-length(p)]
  P <- c(1,cumprod(p1))*p
  sum(P)
}

cal_prob <- function(data_all, data_samp, data_Rain, data_full) {
  ord <- c(
    "Implant, Efavirenz, 1", "Implant, Nevirapine, 1", "Implant, PI-based, 1", "Implant, No ARV, 1",
    "MEFP, Efavirenz, 1"   , "MEFP, Nevirapine, 1"   , "MEFP, PI-based, 1"   , "MEFP, No ARV, 1",
    "Depo, Efavirenz, 1"   , "Depo, Nevirapine, 1"   , "Depo, PI-based, 1"   , "Depo, No ARV, 1",
    "No FP, Efavirenz, 1"  , "No FP, Nevirapine, 1"  , "No FP, PI-based, 1"  , "No FP, No ARV, 1",
    
    "Implant, Efavirenz, 0", "Implant, Nevirapine, 0", "Implant, PI-based, 0", "Implant, No ARV, 0",
    "MEFP, Efavirenz, 0"   , "MEFP, Nevirapine, 0"   , "MEFP, PI-based, 0"   , "MEFP, No ARV, 0",
    "Depo, Efavirenz, 0"   , "Depo, Nevirapine, 0"   , "Depo, PI-based, 0"   , "Depo, No ARV, 0",
    "No FP, Efavirenz, 0"  , "No FP, Nevirapine, 0"  , "No FP, PI-based, 0"  , "No FP, No ARV, 0"
  )
  
  del_records <- c(grep('LEFP', data_all$combination), grep('Oral', data_all$combination),
                   grep('NRTIs', data_all$combination), grep('Other ARV', data_all$combination))
  keep_records <- !(1:nrow(data_all) %in% del_records)
  data_all  <- data_all %>% mutate(index = paste(combination, pregnant, sep=', ')) %>%
    filter(!(is.na(combination) | combination=='' | is.na(pregnant)) & keep_records) %>% group_by(person_id)
  data_all_temp <- matrix(0, nrow=length(unique(data_all$person_id)), ncol=33) %>% as.data.frame(.)
  colnames(data_all_temp) <- c('person_id', ord)
  
  un_ids <- unique(data_all$person_id)
  for (j in 1:length(un_ids)){
    exp_col <- which(colnames(data_all_temp) %in% data_all[data_all$person_id == un_ids[j],]$index)
    data_all_temp[j, c(1, exp_col)] <- c(unique(data_all$person_id)[j], rep(1, length(exp_col)))
    if (j %% 5000 == 0)
      cat('Calc Prob. Person = ', j, ' out of ', length(un_ids), '\n')
  }
  
  data_all_temp_ph2 <- data_all_temp %>% filter(person_id %in% unique(data_samp$person_id))
  
  dataMatch1 <- read.table('Match1.csv', sep=',', header=TRUE)
  dataMatch2 <- read.table('Match2.csv', sep=',', header=TRUE) %>% dplyr::select(record_id, id)
  dataMatch <- rbind(dataMatch1, dataMatch2)
  data_Rain2 <- data_Rain %>% left_join(dataMatch, 'id')# %>% left_join(dataStudy, 'id') %>% filter(study == 'FACES ')

  data_samp  <- data_samp %>% mutate(index = paste(exposure, preg_status, sep=', ')) %>%
    filter(!(is.na(exposure) | is.na(preg_status))) %>% group_by(person_id) %>% slice(1)
  data_samp2 <- data_samp %>% filter(record_id %in% data_Rain2$record_id)
  data_samp_temp <- matrix(0, nrow=nrow(data_samp2), ncol=33) %>% as.data.frame(.)
  colnames(data_samp_temp) <- c('person_id', ord)
  
  ## subset dos selecionados para fase 2
  ### mat2, dimensions [n, 32] ## index (reason why she was selected)
  for (j in 1:length(unique(data_samp2$person_id))){
    data_samp2_temp <- data_samp2[j,]
    exp_col <- which(colnames(data_samp_temp) == data_samp2[j,]$index)
    data_samp_temp[j, c(1,exp_col)] <- c(data_samp2$person_id[j], 1)
    
    row_all <- which(data_all_temp$person_id == data_samp2$person_id[j])
    pos_all <- which(colnames(data_all_temp[row_all, ]) == as.character(data_samp2_temp$index))
    data_all_temp[row_all, pos_all] <- 2
  }
  
  total <- apply(data_all_temp, 2, FUN = function(x) sum(x>0))
  
  for (ii in 1:ncol(data_all_temp)){
    data_all_temp[,ii] <- as.numeric(as.character(data_all_temp[,ii]))
  }
  
  den <- rep(0, 32); den[1] <- total[2]
  data_all_ttemp <- data_all_temp[,-1]
  for (i in 2:32){
    data_all_tempv2 <- as.matrix(data_all_ttemp[data_all_ttemp[,i] > 0, 1:(i-1)])
    nn <- sum(data_all_tempv2 > 1)
    if (i <= 8)
      nn <- sum(apply(data_all_tempv2, 1, FUN = function(x) any(x > 0)))
    den[i] <- total[i+1] - nn	
  }
  p <- apply(data_samp_temp[,-1], 2, sum)/den
  
  for (ii in 1:ncol(data_all_ttemp)){
    data_all_ttemp[,ii] <- as.numeric(as.character(data_all_ttemp[,ii]))
  }
  
  Prob <- matrix(NA, nrow(data_all_temp), 2)
  colnames(Prob) <- c('id', 'Probs')
  for (j in 1:nrow(data_all_temp)){
    exp_id <- which(data_all_temp[j,-1] > 0)
    if (any(exp_id %in% c(1:8))){
      exp_id <- min(exp_id[which(exp_id %in% c(1:8))])
    }
    Prob[j,] <- c(data_all_temp$person_id[j], probs(p, cat=exp_id))
  }
  Prob <- Prob %>% as.data.frame(.) %>% mutate(wgt = 1/Probs)
  return(Prob)
}


rak_func <- function(data_emr, data_chart, data_phone, formula) {
  
  IF_func <- function(fit){
    MM <- model.matrix(fit)
    Score.P   <- fit$y*MM - MM*exp(fit$linear.predictors)
    InfoMat.P <- t(MM*sqrt(exp(fit$linear.predictors))) %*% (MM*sqrt(exp(fit$linear.predictors)))/nrow(MM)
    inffunc <- (Score.P) %*% solve(InfoMat.P)
    inffunc
  }

  if (!isTRUE(all.equal(levels(data_emr$fptype), levels(data_chart$fptype))))
  {
    print('levels on FP differ')
    data_emr$fptype  <- factor(data_emr$fptype, levels=levels(data_chart$fptype))
  }
  if (!isTRUE(all.equal(levels(data_emr$arttype), levels(data_chart$arttype))))
  {
    print('levels on ART differ')
    data_emr$arttype  <- factor(data_emr$arttype, levels=levels(data_chart$arttype))
  } 
  
  mod1 <- glm(formula, data=data_emr, family=poisson(link=log))
  inffun_EMR <- IF_func(mod1)
  colnames(inffun_EMR) <- paste("if", 1:ncol(inffun_EMR), sep="")
  inffun_EMRb <- data.frame(id=data_emr$id, inffun_EMR)
  inffun_EMR_sum  <- aggregate(. ~ id, sum, data=inffun_EMRb)
  totals_EMR <- c(`(Intercept)`=nrow(data_emr), colSums(inffun_EMR))

  mod2 <- glm(formula, data=data_chart, family=poisson(link=log))
  if_chart <- IF_func(mod2)
  if_chartb <- data.frame(id=data_chart$id, if_chart)
  if_chart2 <- aggregate(. ~ id, sum, data=if_chartb)
  colnames(if_chart)  <- paste("if", 1:ncol(if_chart), sep="")
  colnames(if_chart2) <- c('id', paste("if", 1:ncol(if_chart), sep=""))
  totals_chart_record <- c(`(Intercept)`=nrow(data_chart), colSums(if_chart))

  data_chartb_ifs <- data_chart %>% left_join(inffun_EMR_sum, by='id')
  data_phoneb_ifs <- data_phone %>% left_join(if_chart2, by='id') 

  #############################
  #Calibrate phase-2 on phase-1 
  #############################
  if_designR <- svydesign(id = ~id, weights = ~ wgt, data = data_chart)
  index <- match(data_chart$id, inffun_EMR_sum$id)
  if_designR1 <- update(if_designR,
                        if1=inffun_EMR_sum[index,2],  if9=inffun_EMR_sum[index,10], if17=inffun_EMR_sum[index,18], if25=inffun_EMR_sum[index,26], if33=inffun_EMR_sum[index,34], if41=inffun_EMR_sum[index,42],
                        if2=inffun_EMR_sum[index,3], if10=inffun_EMR_sum[index,11], if18=inffun_EMR_sum[index,19], if26=inffun_EMR_sum[index,27], if34=inffun_EMR_sum[index,35], if42=inffun_EMR_sum[index,43],
                        if3=inffun_EMR_sum[index,4], if11=inffun_EMR_sum[index,12], if19=inffun_EMR_sum[index,20], if27=inffun_EMR_sum[index,28], if35=inffun_EMR_sum[index,36],
                        if4=inffun_EMR_sum[index,5], if12=inffun_EMR_sum[index,13], if20=inffun_EMR_sum[index,21], if28=inffun_EMR_sum[index,29], if36=inffun_EMR_sum[index,37],
                        if5=inffun_EMR_sum[index,6], if13=inffun_EMR_sum[index,14], if21=inffun_EMR_sum[index,22], if29=inffun_EMR_sum[index,30], if37=inffun_EMR_sum[index,38],
                        if6=inffun_EMR_sum[index,7], if14=inffun_EMR_sum[index,15], if22=inffun_EMR_sum[index,23], if30=inffun_EMR_sum[index,31], if38=inffun_EMR_sum[index,39],
                        if7=inffun_EMR_sum[index,8], if15=inffun_EMR_sum[index,16], if23=inffun_EMR_sum[index,24], if31=inffun_EMR_sum[index,32], if39=inffun_EMR_sum[index,40],
                        if8=inffun_EMR_sum[index,9], if16=inffun_EMR_sum[index,17], if24=inffun_EMR_sum[index,25], if32=inffun_EMR_sum[index,33], if40=inffun_EMR_sum[index,41])
  keep_cols_efv1 <- union(union(1:9, grep('arttypeEFV', rownames(summary(mod1)$coef))), grep('fptypeDepo', rownames(summary(mod1)$coef)))
  keep_cols_efv2 <- c(1, keep_cols_efv1+1)
  
  keep_cols_efv1 <- keep_cols_efv1[keep_cols_efv1<42]
  keep_cols_efv2 <- keep_cols_efv2[keep_cols_efv2<=42]
  keep_cols_efv1 <- c(keep_cols_efv1[keep_cols_efv1<10], keep_cols_efv1[keep_cols_efv1<42 & keep_cols_efv1>36])
  keep_cols_efv2 <- c(keep_cols_efv2[keep_cols_efv2<=10], keep_cols_efv2[keep_cols_efv2<=42 & keep_cols_efv2>36])
  
  if_formula1 <- as.formula(paste('~ ', paste(colnames(if_chart[,keep_cols_efv1]), collapse=' + ')))
  if_cal2ph <- survey::calibrate(design=if_designR1, calfun='raking', formula=if_formula1, population=totals_EMR[keep_cols_efv2])
  resCAL_2ph <- svyglm(formula, design=if_cal2ph, family=poisson(link=log))

  #############################
  #Calibrate phase-3 on phase-1 
  #############################
  if_design <- svydesign(id = ~id, weights = ~wgt, data = data_phone)
  index <- match(data_phone$id, inffun_EMR_sum$id)
  if_design1 <- update(if_design,
                       if1=inffun_EMR_sum[index,2],  if9=inffun_EMR_sum[index,10], if17=inffun_EMR_sum[index,18], if25=inffun_EMR_sum[index,26], if33=inffun_EMR_sum[index,34], if41=inffun_EMR_sum[index,42],
                       if2=inffun_EMR_sum[index,3], if10=inffun_EMR_sum[index,11], if18=inffun_EMR_sum[index,19], if26=inffun_EMR_sum[index,27], if34=inffun_EMR_sum[index,35], if42=inffun_EMR_sum[index,43],
                       if3=inffun_EMR_sum[index,4], if11=inffun_EMR_sum[index,12], if19=inffun_EMR_sum[index,20], if27=inffun_EMR_sum[index,28], if35=inffun_EMR_sum[index,36],
                       if4=inffun_EMR_sum[index,5], if12=inffun_EMR_sum[index,13], if20=inffun_EMR_sum[index,21], if28=inffun_EMR_sum[index,29], if36=inffun_EMR_sum[index,37],
                       if5=inffun_EMR_sum[index,6], if13=inffun_EMR_sum[index,14], if21=inffun_EMR_sum[index,22], if29=inffun_EMR_sum[index,30], if37=inffun_EMR_sum[index,38],
                       if6=inffun_EMR_sum[index,7], if14=inffun_EMR_sum[index,15], if22=inffun_EMR_sum[index,23], if30=inffun_EMR_sum[index,31], if38=inffun_EMR_sum[index,39],
                       if7=inffun_EMR_sum[index,8], if15=inffun_EMR_sum[index,16], if23=inffun_EMR_sum[index,24], if31=inffun_EMR_sum[index,32], if39=inffun_EMR_sum[index,40],
                       if8=inffun_EMR_sum[index,9], if16=inffun_EMR_sum[index,17], if24=inffun_EMR_sum[index,25], if32=inffun_EMR_sum[index,33], if40=inffun_EMR_sum[index,41])
  if_cal_3p  <- survey::calibrate(design=if_design1, calfun='raking', formula=if_formula1, population=totals_EMR[keep_cols_efv2], aggregate.index=~id)
  resCAL_3ph <- svyglm(formula, design=if_cal_3p, family=poisson(link=log))

  #############################
  #Calibrate phase-3 on phase-2 
  #############################
  data_phoneb <- data_phone %>% filter(id %in% if_chart2$id)
  if_design <- svydesign(id = ~id, weights = ~wgt2, data = data_phoneb)
  index <- match(data_phoneb$id, if_chart2$id)
  if_design1 <- update(if_design,
                       if1=if_chart2[index,2],  if9=if_chart2[index,10], if17=if_chart2[index,18], if25=if_chart2[index,26], if33=if_chart2[index,34], if41=if_chart2[index,42],
                       if2=if_chart2[index,3], if10=if_chart2[index,11], if18=if_chart2[index,19], if26=if_chart2[index,27], if34=if_chart2[index,35], if42=if_chart2[index,43],
                       if3=if_chart2[index,4], if11=if_chart2[index,12], if19=if_chart2[index,20], if27=if_chart2[index,28], if35=if_chart2[index,36],
                       if4=if_chart2[index,5], if12=if_chart2[index,13], if20=if_chart2[index,21], if28=if_chart2[index,29], if36=if_chart2[index,37],
                       if5=if_chart2[index,6], if13=if_chart2[index,14], if21=if_chart2[index,22], if29=if_chart2[index,30], if37=if_chart2[index,38],
                       if6=if_chart2[index,7], if14=if_chart2[index,15], if22=if_chart2[index,23], if30=if_chart2[index,31], if38=if_chart2[index,39],
                       if7=if_chart2[index,8], if15=if_chart2[index,16], if23=if_chart2[index,24], if31=if_chart2[index,32], if39=if_chart2[index,40],
                       if8=if_chart2[index,9], if16=if_chart2[index,17], if24=if_chart2[index,25], if32=if_chart2[index,33], if40=if_chart2[index,41])
  if_cal_23p  <- survey::calibrate(design=if_design1, calfun='raking', formula=if_formula1, population=totals_chart_record[keep_cols_efv2], aggregate.index=~id)
  resCAL_23ph <- svyglm(formula, design=if_cal_23p, family=poisson(link=log))

  ################################################
  # Raking, 3-phase, product of calibrated weights
  ################################################
  if_data_phase2 <- data.frame(id=if_cal2ph$cluster$id, p1=if_cal2ph$prob) %>% filter(id %in% if_cal_23p$cluster$id) %>%
    group_by(id) %>% slice(1)
  if_data_phase3 <- data.frame(id=if_cal_23p$cluster$id, p2=if_cal_23p$prob) %>% group_by(id) %>% slice(1)
  pi2_new <- if_data_phase3 %>% left_join(if_data_phase2, by='id') %>% mutate(pprod_new = p1*p2)
  data_phoneb <- data_phoneb %>% left_join(pi2_new, by='id')
  
  data_chart_ifs2  <- data_chart  %>% left_join(inffun_EMR_sum, by='id')
  data_phoneb_ifs2 <- data_phoneb %>% left_join(if_chart2, by='id')
  raking_new  <- glm(formula, weights = 1/pprod_new, data=data_phoneb_ifs2, family=poisson(link=log))
  
  rr <- residuals(raking_new)
  XX <- model.matrix(raking_new) * weights(raking_new, "working")
  r <- residuals(raking_new, "working") * sqrt(weights(raking_new, "working"))
  X <- model.matrix(raking_new) * sqrt(weights(raking_new, "working"))
  I0 <- X*residuals(raking_new)
  
  m12 <- glm(formula, weights=1/Probs, family=poisson(link=log), data=data_phoneb_ifs2)
  ifd1 <- survey::svydesign(id = ~id, probs = ~ Probs, data = data_phoneb_ifs2)
  m12  <- svyglm(formula, family=poisson(link=log), design=ifd1)  
  r12 <- residuals(m12, "response") * sqrt(weights(m12, "prior"))
  x12 <- model.matrix(m12) * sqrt(weights(m12, "prior"))
  ifs12 <- as.matrix(cbind(1, data_phoneb_ifs2[, paste0('if', keep_cols_efv2)]))
  t12 <- t(ifs12 * sqrt(weights(m12, "prior"))) %*% (ifs12 * sqrt(weights(m12, "prior")))
  ifs12 <- as.matrix(cbind(1, data_phoneb_ifs2[, paste0('if', keep_cols_efv2)]))
  ifs12b <- as.matrix(cbind(1, data_phoneb_ifs2[, paste0('if', keep_cols_efv2)]))
  yy12 <- t(ifs12) %*% (r12 * x12)
  r12T <- (ifs12b * weights(m12, "prior")) %*% solve(t12) %*% yy12
  
  r12_temp <- data.frame(id = data_phoneb_ifs2$id, r12T) %>% group_by(id) %>% slice(1)
  r23_temp <- data.frame(id=data_phoneb_ifs2$id) %>% left_join(r12_temp, by='id') %>% dplyr::select(-id)
  r12T <- r23_temp
  
  m23 <- glm(formula, weights=1/p2, family=poisson(link=log), data=data_phoneb_ifs2)
  ifd <- survey::svydesign(id = ~id, probs = ~ p2, data = data_phoneb_ifs2)
  m23 <- svyglm(formula, family=poisson(link=log), design=ifd)
  r23 <- residuals(m23, "response") * sqrt(weights(m23, "prior"))
  X23 <- model.matrix(m23) * sqrt(weights(m23, "prior"))
  ifs23 <- as.matrix(cbind(1, data_phoneb_ifs2[, paste0('if', keep_cols_efv2)]))
  t23 <- t(ifs23 * sqrt(weights(m23, "prior"))) %*% (ifs23 * sqrt(weights(m23, "prior")))
  yy23  <- t(ifs23) %*% (r23 * X23)
  r23 <- (ifs23 * weights(m23, "prior")) %*% solve(t23) %*% yy23
  
  varsRak_new <- solve(t(X) %*% X) %*% (t(as.matrix(I0 - r12T - r23 + r12T*r23)) %*%
                                 as.matrix(I0 - r12T - r23 + r12T*r23)) %*% solve(t(X) %*% X)

  resCAL2ph <- resCAL_2ph
  resCAL3ph <- resCAL_3ph
  raking_23ph  <- raking_new
  varsRak_23ph <- varsRak_new
  
  return(list(rescal2ph=resCAL2ph, rescal3ph=resCAL3ph,
              raking3phb=raking_23ph, varsRak3phb=varsRak_23ph))
}


rm(list = ls())
library(tidyverse)
library(parallel)
library(survey)
source("zmisc.R")
set.seed(02138)

earningstypes <- c("log_1000")
var <- "chdearn2"

sample <- list()
sample[[1]] <- readRDS(file = paste0("Samples/", var, "drop_full.RDS"))
names(sample) <- c("drop_full")

for (earn in earningstypes) {
  
  nlsy97_full <- sample[[1]] %>%
    map( ~ mutate(.x, 
                  hs = if_else(educAB2_grad >= 2, 1, 0),
                  college = if_else(educAB2_grad >= 4, 1, 0),
                  ba = if_else(educAB2_grad >= 5, 1, 0),
                  gschool = if_else(educAB2_grad >= 6, 1, 0)))
  
  if (earn == "no_log") {
    nlsy97_full <- nlsy97_full %>%
      map( ~ mutate(.x, earnings = !!sym(var)))
  } else if (earn == "log_10") {
    nlsy97_full <- nlsy97_full %>%
      map( ~ mutate(.x, earnings = log(!!sym(var) + 10)))
  } else if (earn == "log_100") {
    nlsy97_full <- nlsy97_full %>%
      map( ~ mutate(.x, earnings = log(!!sym(var) + 100)))
  } else if (earn == "log_1000") {
    nlsy97_full <- nlsy97_full %>%
      map( ~ mutate(.x, earnings = log(!!sym(var) + 1000)))
  }
  
  
  ##########################################################
  # Formulas for treatment and outcome models
  ##########################################################
  
  a <- "hs"
  m1 <- "college"
  m2 <- "ba"
  m3 <- "gschool"
  y <- "earnings"
  w <- "weight"
  id <- "id"
  
  x <- nlsy97_full[[1]] %>%
    dplyr::select(female, black, hispan, parinc, paredu, parast, livebio, father, #other
                  afqt3, gpa, children18, substind, delinqind, rural97, south97,
                  frcoll_75, frcoll_90, schstolen, schthreat, schfight) %>%
    names()
  
  z <- nlsy97_full[[1]] %>%
    dplyr::select(stem, college_gpa) %>%
    names()
  
  y_form_ax <- as.formula(paste(y, " ~ (", paste(c(x), collapse= "+"), ")* (", a, ")"))
  y_form_axm1 <- as.formula(paste(y, " ~ (", paste(c(x), collapse= "+"), ")* (", a, "+", m1, ")"))
  y_form_axzm1m2 <- as.formula(paste(y, " ~  (", paste(c(x), collapse= "+"), ")* (", a, "+", m1, "+", m2, ") + ", paste(c(z), collapse= "+"), "+", m2))
  y_form_axzm1m2m3 <- as.formula(paste(y, " ~  (", paste(c(x), collapse= "+"), ")* (", a, "+", m1, "+", m2,  "+", m3,") + ", paste(c(z), collapse= "+"), "+ (", m2, "+", m3, ")"))

  m1_form <- as.formula(paste(m1, " ~ ", paste(c(x), collapse= "+")))
  m2_form <- as.formula(paste(m2, " ~ ", paste(c(x), collapse= "+"))) # M1 -> m2
  m3_form <- as.formula(paste(m3, " ~ (", paste(c(x), collapse= "+"), ")+ (", paste(c(z), collapse= "+"), ")"))
  
  z1_form <- as.formula(paste("college_gpa  ~ ", paste(c(x), collapse= "+")))
  z2_form <- as.formula(paste("stem  ~ ", paste(c(x), collapse= "+")))

  ##########################################################
  # Create bootstrap function to parallelise
  ##########################################################
  
  boot <- function(d) {
    s <- sample(1:nrow(d), replace = TRUE)
    df_0 <- d[s,] # df_0 <- nlsy97_full[[1]]
    
    design_z <- svydesign(ids = ~ 1, data = filter(df_0, college == 1), weights = ~ weight)
    
    # Residualise intermediate confounders
    z1.mod <- svyglm(z1_form, design = design_z)
    z2.mod <- svyglm(z2_form, design = design_z)
    
    df <- df_0 %>%
      mutate_at(vars(x), demean) %>%
      mutate(college_gpa = if_else(college == 1, residualize2(z1.mod, college_gpa, .), 0),
             stem = if_else(college == 1, residualize2(z2.mod, stem, .), 0))
    
    design_rwr <- svydesign(ids = ~ 1, data = df, weights = ~ weight)
    design_rwr_m1 <- svydesign(ids = ~ 1, data = filter(df, hs == 1), weights = ~ weight) # ie M1 ~ X | A=1 (hs)
    design_rwr_m2 <- svydesign(ids = ~ 1, data = filter(df, college == 1), weights = ~ weight)# ie M2 ~ X  | M1=1 (college)
    design_rwr_m3 <- svydesign(ids = ~ 1, data = filter(df, ba == 1), weights = ~ weight)
    
    rwr_y_ax <- svyglm(y_form_ax, design = design_rwr)
    rwr_y_xm1 <- svyglm(y_form_axm1, design = design_rwr)
    rwr_y_axzm1m2 <- svyglm(y_form_axzm1m2, design = design_rwr)
    rwr_y_axzm1m2m3 <- svyglm(y_form_axzm1m2m3, design = design_rwr)
    
    rwr_m1 <- svyglm(m1_form, design = design_rwr_m1)
    rwr_m2 <- svyglm(m2_form, design = design_rwr_m2)
    rwr_m3 <- svyglm(m3_form, design = design_rwr_m3)
    
    boot_out <- tibble(
      # Delta terms
      delta0_rwr = rwr_y_xm1$coefficients[["hs"]],  # E[Y(1,0)-Y(0,0)]
      delta1_rwr = rwr_y_axzm1m2$coefficients[["college"]], # E[Y(1,1,0)-Y(1,0,0)]
      delta2_rwr = rwr_y_axzm1m2m3$coefficients[["ba"]],  # E[Y(1,1,1,0)-Y(1,1,0,0)]
      delta3_rwr = rwr_y_axzm1m2m3$coefficients[["gschool"]], # E[Y(1,1,1,1)-Y(1,1,1,0)] 
      
      # Tau terms
      tau0_rwr = rwr_y_ax$coefficients[["hs"]],  # ATE
      tau1_rwr = rwr_y_xm1$coefficients[["college"]], # E[Y(1,1)-Y(1,0)]
      tau2_rwr = rwr_y_axzm1m2$coefficients[["ba"]],  # E[Y(1,1,1)-Y(1,1,0)]
      tau3_rwr = delta3_rwr,
      
      pi1_rwr = rwr_m1$coefficients[["(Intercept)"]],
      pi2_rwr = rwr_m2$coefficients[["(Intercept)"]],
      pi3_rwr = rwr_m3$coefficients[["(Intercept)"]],
      
      # Eta terms
      eta1_rwr = tau0_rwr - delta0_rwr - pi1_rwr*tau1_rwr,
      eta2_rwr = tau1_rwr - delta1_rwr - pi2_rwr*tau2_rwr,
      eta3_rwr = tau2_rwr - delta2_rwr - pi3_rwr*tau3_rwr,
      
      # Direct and continuation effects (theta)
      theta0_rwr = delta0_rwr,
      theta1_rwr = pi1_rwr*delta1_rwr + eta1_rwr,
      theta2_rwr = pi1_rwr*pi2_rwr*delta2_rwr + pi1_rwr*eta2_rwr,
      theta3_rwr = pi1_rwr*pi2_rwr*pi3_rwr*delta3_rwr + pi1_rwr*pi2_rwr*eta3_rwr) %>% t() #%>% as.data.frame()
    
    
    return(boot_out)
  }
  ##########################################################
  # Main analyses
  ##########################################################
  no_of_cores <- detectCores()
  my_cores_to_use <- no_of_cores - 2
  cl <- makeCluster(my_cores_to_use)
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(survey))
  clusterSetRNGStream(cl, iseed = 02138)
  
  I <- length(nlsy97_full)
  
  out <- vector(mode = "list", I)
  overall_out <-  vector(mode = "list", I)
  reps <- 1000
  reps <- 250
  
  for(i in 1:I){  
    
    cat("imputed sample ", i, "\n")
    
    df_i <- nlsy97_full[[i]]
    clusterExport(cl, c("boot", "df_i", "residualize2", "y_form_ax", "y_form_axm1", "y_form_axzm1m2", "y_form_axzm1m2m3", 'm1_form', "m2_form", "m3_form", 
                        "z1_form", "z2_form", "x", "demean"))
    out <- parLapply(cl, 1:reps, function(i) boot(df_i))
    out2 <- do.call(cbind, out)
    
    output <- data.frame(estimand = names(out2[,1]),
                         est = apply(out2, MARGIN = 1, mean),
                         se = apply(out2, MARGIN = 1, sd))
    
    overall_out[[i]] <- output 
    
  }
  
  stopCluster(cl)
  
  ############################
  # Output
  ############################
  
  overall_df_rwr <- overall_out %>%
    imap(., ~ mutate(.x, index = .y)) %>%
    reduce(., bind_rows) %>%
    group_by(estimand) %>%
    summarise(coef_est = mean(est),
              coef_se = sqrt(mean(se^2) + (1 + 1/I)*var(est)),
    ) %>%
    ungroup() %>%
    pivot_longer(names_to = "measure", coef_est:coef_se) %>%
    separate(measure, c("type", "measure")) %>%
    pivot_wider(names_from = "measure", values_from = "value") %>%
    select(-type) %>% separate(estimand, sep = "_", into = c("estimand", "estimator")) %>% dplyr::select(-estimator)
  

  name <- paste0("02_rwr_", earn, ".RDS")
  saveRDS(overall_df_rwr, file = name)
  
  
}
  
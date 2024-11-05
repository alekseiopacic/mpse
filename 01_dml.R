rm(list = ls())
library(tidyverse)
library(survey)
library(rlang)
library(caret)
library(glmnet)
library(ranger)
library(SuperLearner)
source("zmisc.R")
set.seed(02138)

vars <- c("chdearn2") 
earningstypes <- c("log_1000")
for (var in vars) {
  
  sample <- readRDS(file = paste0("Samples/", var, "drop_full.RDS"))
  
  for (earn in earningstypes) {
    
    nlsy97_full <- sample %>%
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
      dplyr::select(female, black, hispan, parinc, paredu, parast, livebio, father, re!
                    afqt3, gpa, children18, substind, delinqind, rural97, south97,
                    frcoll_75, frcoll_90, schstolen, schthreat, schfight) %>%
      names()
    
    z <- nlsy97_full[[1]] %>%
      dplyr::select(stem, college_gpa) %>%
      names()
    
    
    # Remember to drop constant variables for model matrices
    form_x <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2, m3), collapse= "+")))
    form_xz <- as.formula(paste(y, " ~ ", paste(c(x, z,  a, m1, m2, m3), collapse= "+")))
    
    form_x.w <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2, m3, w), collapse= "+")))
    form_xz.w <- as.formula(paste(y, " ~ ", paste(c(x, z,  a, m1, m2, m3, w), collapse= "+")))
    
    form_xz.id <- as.formula(paste(y, " ~ ", paste(c(x, z,  a, m1, m2, m3, id), collapse= "+")))
    
    
    ##########################################################
    # Main analyses
    ##########################################################
    
    I <- length(nlsy97_full)
    
    K <- 5
    
    out <- vector(mode = "list", I)
    overall_out <-  vector(mode = "list", I)
    
    for(i in 1:I){
      
      cat("imputed sample ", i, "\n")
      
      df_0 <- nlsy97_full[[i]] 
      
      rownames(df_0) <- df_0$id
    
      
      # create cross-fitting split
      cf_fold <- createFolds(df_0$earnings, K)  
      
      main_list <- vector(mode = "list", K)
      
      # k = 1
      for(k in 1:K){
        
        cat(" cross-fitting fold ", k, "\n")
        
        df <- df_0
        
        df_x <- model.matrix(form_x, data = df)[, -1] %>% as_tibble()
        df_x.w <- model.matrix(form_x.w, data = df)[, -1] %>% as_tibble()
        
        df_xz <- model.matrix(form_xz, data = df)[, -1] %>% as_tibble()
        df_xz.w <- model.matrix(form_xz.w, data = df)[, -1] %>% as_tibble()
        
        #################################################
        # Design matrices : ensure filter before taking fold  
        #################################################
        
        aux <- df[-cf_fold[[k]], ]
        
        aux_x <- df_x[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_x.w <- df_x.w[-cf_fold[[k]], ]
        
        aux_a <- dplyr::filter(df, hs == 1)[-cf_fold[[k]], ]
        aux_an <- dplyr::filter(df, hs == 0)[-cf_fold[[k]], ]
        
        aux_x_a <- dplyr::filter(df_x, hs == 1)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_x_a.w <- dplyr::filter(df_x.w, hs == 1)[-cf_fold[[k]], ] 
        
        aux_x_an <- dplyr::filter(df_x, hs == 0)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_x_an.w <- dplyr::filter(df_x.w, hs == 0)[-cf_fold[[k]], ] 
        
        aux_a_m1 <- dplyr::filter(df, hs == 1, college == 1)[-cf_fold[[k]], ]
        aux_x_a_m1 <- dplyr::filter(df_x, hs == 1, college == 1)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_x_a_m1.w <- dplyr::filter(df_x.w, hs == 1, college == 1)[-cf_fold[[k]], ] 
        
        aux_a_m1n <- dplyr::filter(df, hs == 1, college == 0)[-cf_fold[[k]], ]
        aux_x_a_m1n <- dplyr::filter(df_x, hs == 1, college == 0)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_x_a_m1n.w <- dplyr::filter(df_x.w, hs == 1, college == 0)[-cf_fold[[k]], ] 
        
        
        aux_xz_a_m1 <- dplyr::filter(df_xz, hs == 1, college == 1)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_xz_a_m1.w <- dplyr::filter(df_xz.w, hs == 1, college == 1)[-cf_fold[[k]], ] 
        
        aux_a_m1_m2 <- dplyr::filter(df, hs == 1, college == 1, ba == 1)[-cf_fold[[k]], ]
        aux_xz_a_m1_m2 <- dplyr::filter(df_xz, hs == 1, college == 1, ba == 1)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_xz_a_m1_m2.w <- dplyr::filter(df_xz.w, hs == 1, college == 1, ba == 1)[-cf_fold[[k]], ] 
        
        aux_a_m1_m2n <- dplyr::filter(df, hs == 1, college == 1, ba == 0)[-cf_fold[[k]], ]
        aux_xz_a_m1_m2n <- dplyr::filter(df_xz, hs == 1, college == 1, ba == 0)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_xz_a_m1_m2n.w <- dplyr::filter(df_xz.w, hs == 1, college == 1, ba == 0)[-cf_fold[[k]], ] 
        
        
        aux_x_a_m1_m2 <- dplyr::filter(df_x, hs == 1, college == 1, ba == 1)[-cf_fold[[k]], ] %>%  dplyr::select(-c(hs, college, ba, gschool))
        aux_x_a_m1_m2.w <- dplyr::filter(df_x.w, hs == 1, college == 1, ba == 1)[-cf_fold[[k]], ] %>%  dplyr::select(-c(hs, college, ba, gschool))
        
        aux_x_a_m1_m2n <- dplyr::filter(df_x, hs == 1, college == 1, ba == 0)[-cf_fold[[k]], ] %>%  dplyr::select(-c(hs, college, ba, gschool))
        aux_x_a_m1_m2n.w <- dplyr::filter(df_x.w, hs == 1, college == 1, ba == 0)[-cf_fold[[k]], ] %>%  dplyr::select(-c(hs, college, ba, gschool))
        
        
        aux_a_m1_m2_m3 <- dplyr::filter(df, hs == 1, college == 1, ba == 1, gschool == 1)[-cf_fold[[k]], ]
        aux_xz_a_m1_m2_m3 <- dplyr::filter(df_xz, hs == 1, college == 1, ba == 1, gschool == 1)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_xz_a_m1_m2_m3.w <- dplyr::filter(df_xz.w, hs == 1, college == 1, ba == 1, gschool == 1)[-cf_fold[[k]], ] 
        
        aux_a_m1_m2_m3n <- dplyr::filter(df, hs == 1, college == 1, ba == 1, gschool == 0)[-cf_fold[[k]], ]
        aux_xz_a_m1_m2_m3n <- dplyr::filter(df_xz, hs == 1, college == 1, ba == 1, gschool == 0)[-cf_fold[[k]], ] %>% dplyr::select(-c(hs, college, ba, gschool))
        aux_xz_a_m1_m2_m3n.w <- dplyr::filter(df_xz.w, hs == 1, college == 1, ba == 1, gschool == 0)[-cf_fold[[k]], ] 
        
        ## Dfs for prediction
        df_x_pred <- model.matrix(form_x, data = df)[, -1] %>% as_tibble() %>% dplyr::select(-c(hs, college, ba, gschool))
        
        df_xz_pred_a_m1 <- model.matrix(form_xz.id, data = df)[, -1] %>% as.data.frame()  
        df_xz_pred_a_m1 <- df_xz_pred_a_m1 %>%
          `rownames<-`(df_xz_pred_a_m1$id) %>% 
          filter(hs == 1, college == 1) %>% dplyr::select(-c(hs, college, ba, gschool, id))
        
        df_xz_pred_a_m1_m2 <- model.matrix(form_xz.id, data = df)[, -1] %>% as.data.frame() 
        df_xz_pred_a_m1_m2 <- df_xz_pred_a_m1_m2 %>%
          `rownames<-`(df_xz_pred_a_m1_m2$id) %>%
          filter(hs == 1, college == 1, ba == 1) %>% dplyr::select(-c(hs, college, ba, gschool, id))
        
        df_xz_pred_a_m1_m2_m3 <- model.matrix(form_xz.id, data = df)[, -1] %>% as.data.frame() 
        df_xz_pred_a_m1_m2_m3 <- df_xz_pred_a_m1_m2_m3 %>%
          `rownames<-`(df_xz_pred_a_m1_m2_m3$id) %>%
          filter(hs == 1, college == 1, ba == 1, gschool == 1) %>% dplyr::select(-c(hs, college, ba, gschool, id))
        
  
        #################################################
        # Treatment and Mediator Models
        #################################################
        
        a_sl <- SuperLearner(
          Y          = aux$hs,
          X          = aux_x,
          newX       = df_x_pred,
          family     = binomial(),
          obsWeights = aux$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
          cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
        )
        
        df$a_fit = predict.SuperLearner(a_sl, newdata = df_x_pred)$pred
        
        m1_sl <- SuperLearner(
          Y          = aux_a$college,
          X          = aux_x_a,
          newX       = df_x_pred,
          family     = binomial(),
          obsWeights = aux_x_a.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
          cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
        )
        
        df$m1_fit = predict.SuperLearner(m1_sl, newdata = df_x_pred)$pred
        
        ## M2 as outcome model (no Z)
        m2_sl_noz <- SuperLearner(
          Y          = aux_a_m1$ba,
          X          = aux_x_a_m1,
          newX       = df_x_pred,
          family     = binomial(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
          cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
        )
        
        df$m2_fit_noz = predict.SuperLearner(m2_sl_noz, newdata = df_x_pred)$pred
        
        ## Models with  predictions only required of a subset of the data ##
        
        # M2 for IPW  (with Z)
        m2_sl <- SuperLearner(
          Y          = aux_a_m1$ba,
          X          = aux_xz_a_m1,
          newX       = df_xz_pred_a_m1,
          family     = binomial(),
          obsWeights = aux_xz_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
          cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
        )
        
        df_xz_pred_to.merge <- df_xz_pred_a_m1
        df_xz_pred_to.merge$m2_fit <- predict.SuperLearner(m2_sl, newdata = df_xz_pred_a_m1)$pred
        df_xz_pred_to.merge <- df_xz_pred_to.merge %>% tibble::rownames_to_column(var = "id") %>% dplyr::select(id, m2_fit)
        df <- df %>% merge(df_xz_pred_to.merge, by = "id", all.x=TRUE)
        
        ## M3 outcome model: fitted to and imputed for units with hs == 1, college == 1, ba == 1
        
        m3_sl <- SuperLearner(
          Y          = aux_a_m1_m2$gschool,
          X          = aux_xz_a_m1_m2,
          newX       = df_xz_pred_a_m1, # predict among all with M1 and higher
          family     = binomial(),
          obsWeights = aux_xz_a_m1_m2.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
          cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
        )
        
        # also need to create an m3fit var in the full dataframe 
        
        df_xz_pred_to.merge <- df_xz_pred_a_m1
        df_xz_pred_to.merge$m3_fit <- predict.SuperLearner(m3_sl, newdata = df_xz_pred_a_m1)$pred
        df_xz_pred_to.merge2 <- df_xz_pred_to.merge %>% tibble::rownames_to_column(var = "id") %>% dplyr::select(id, m3_fit)
        df <- df %>% merge(df_xz_pred_to.merge2, by = "id", all.x=TRUE)
        
        m3_sl.x <- SuperLearner(
          Y          = df_xz_pred_to.merge$m3_fit[-cf_fold[[k]]],
          X          = aux_x_a_m1,
          newX       = df_x_pred, #technically shoudn't be all units but indicator in IPW will make unwanted values 0 
          family     = gaussian(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$m3_fit.x <- predict.SuperLearner(m3_sl.x, newdata = df_x_pred)$pred
        
        
        #################################################
        # Outcome models
        #################################################
        
        #### Y(1,1,1,m3)
        
        # Y(1,1,1,1)
        y1111_sl <- SuperLearner(
          Y          = aux_a_m1_m2_m3$earnings,
          X          = aux_xz_a_m1_m2_m3,
          newX       = df_xz_pred_a_m1, # predict among all with M1 and higher
          family     = gaussian(),
          obsWeights = aux_xz_a_m1_m2_m3.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df_xz_pred_to.merge <- df_xz_pred_a_m1
        df_xz_pred_to.merge$y1111_fit <- predict.SuperLearner(y1111_sl, newdata = df_xz_pred_a_m1)$pred
        df_xz_pred_to.merge2 <- df_xz_pred_to.merge %>% tibble::rownames_to_column(var = "id") %>% dplyr::select(id, y1111_fit)
        df <- df %>% merge(df_xz_pred_to.merge2, by = "id", all.x=TRUE)
        
        y1111_sl.x <- SuperLearner(
          Y          = df_xz_pred_to.merge$y1111_fit[-cf_fold[[k]]],
          X          = aux_x_a_m1,
          newX       = df_x_pred,  #technically shoudn't be all units but indicator in IPW will make unwanted values 0 
          family     = gaussian(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y1111_fit.x <- predict.SuperLearner(y1111_sl.x, newdata = df_x_pred)$pred
        
        ###
        
        # Y(1,1,1,0)
        y1110_sl <- SuperLearner(
          Y          = aux_a_m1_m2_m3n$earnings,
          X          = aux_xz_a_m1_m2_m3n,
          newX       = df_xz_pred_a_m1, # predict among all with M1 and higher
          family     = gaussian(),
          obsWeights = aux_xz_a_m1_m2_m3n.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df_xz_pred_to.merge <- df_xz_pred_a_m1
        df_xz_pred_to.merge$y1110_fit <- predict.SuperLearner(y1110_sl, newdata = df_xz_pred_a_m1)$pred
        df_xz_pred_to.merge2 <- df_xz_pred_to.merge %>% tibble::rownames_to_column(var = "id") %>% dplyr::select(id, y1110_fit)
        df <- df %>% merge(df_xz_pred_to.merge2, by = "id", all.x=TRUE)
        
        y1110_sl.x <- SuperLearner(
          Y          = df_xz_pred_to.merge$y1110_fit[-cf_fold[[k]]],
          X          = aux_x_a_m1,
          newX       = df_x_pred,  #technically shoudn't be all units but indicator in IPW will make unwanted values 0 
          family     = gaussian(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y1110_fit.x <- predict.SuperLearner(y1111_sl.x, newdata = df_x_pred)$pred
        
        
        #################################################
        #### Y(a), Y(1,m1), Y(1,1,m2)
        #################################################
        
        # Y(1)
        y1_sl <- SuperLearner(
          Y          = aux_a$earnings,
          X          = aux_x_a,
          family     = gaussian(),
          obsWeights = aux_x_a.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y1_fit <- predict.SuperLearner(y1_sl, newdata = df_x_pred)$pred
        
        # Y(0)
        y0_sl <- SuperLearner(
          Y          = aux_an$earnings,
          X          = aux_x_an,
          family     = gaussian(),
          obsWeights = aux_x_an.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y0_fit <- predict.SuperLearner(y0_sl, newdata = df_x_pred)$pred
        
        
        # Y(1,1)
        y11_sl <- SuperLearner(
          Y          = aux_a_m1$earnings,
          X          = aux_x_a_m1,
          newX       = df_x_pred,
          family     = gaussian(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y11_fit <- predict.SuperLearner(y11_sl, newdata = df_x_pred)$pred
        
        # Y(1,0)
        y10_sl <- SuperLearner(
          Y          = aux_a_m1n$earnings,
          X          = aux_x_a_m1n,
          family     = gaussian(),
          obsWeights = aux_x_a_m1n.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y10_fit <- predict.SuperLearner(y10_sl, newdata = df_x_pred)$pred
        
        
        # Y(1,1,1) and fitted model
        
        y111_sl <- SuperLearner(
          Y          = aux_a_m1_m2$earnings,
          X          = aux_xz_a_m1_m2,
          newX = df_xz_pred_a_m1, # predict among all with M1 and higher - i believe this is correct
          family     = gaussian(),
          obsWeights = aux_xz_a_m1_m2.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        
        # also need to create an predition var in the full dataframe 
        
        df_xz_pred_to.merge <- df_xz_pred_a_m1
        df_xz_pred_to.merge$y111_fit <- predict.SuperLearner(y111_sl, newdata = df_xz_pred_a_m1)$pred
        df_xz_pred_to.merge2 <- df_xz_pred_to.merge %>% tibble::rownames_to_column(var = "id") %>% dplyr::select(id, y111_fit)
        df <- df %>% merge(df_xz_pred_to.merge2, by = "id", all.x=TRUE)
        
        y111_sl.x <- SuperLearner(
          Y          = df_xz_pred_to.merge$y111_fit[-cf_fold[[k]]],
          X          = aux_x_a_m1,
          newX       = df_x_pred, #technically shoudn't be all units but indicator in IPW will make unwanted values 0 
          family     = gaussian(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y111_fit.x <- predict.SuperLearner(y111_sl.x, newdata = df_x_pred)$pred
        
        
        
        # Y(1,1,0)
        y110_sl <- SuperLearner(
          Y          = aux_a_m1_m2n$earnings,
          X          = aux_xz_a_m1_m2n,
          newX = df_xz_pred_a_m1, # predict among all with M1 and higher - i believe this is correct
          family     = gaussian(),
          obsWeights = aux_xz_a_m1_m2n.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        
        # also need to create an predition var in the full dataframe 
        
        df_xz_pred_to.merge <- df_xz_pred_a_m1
        df_xz_pred_to.merge$y110_fit <- predict.SuperLearner(y110_sl, newdata = df_xz_pred_a_m1)$pred
        df_xz_pred_to.merge2 <- df_xz_pred_to.merge %>% tibble::rownames_to_column(var = "id") %>% dplyr::select(id, y110_fit)
        df <- df %>% merge(df_xz_pred_to.merge2, by = "id", all.x=TRUE)
        
        y110_sl.x <- SuperLearner(
          Y          = df_xz_pred_to.merge$y110_fit[-cf_fold[[k]]],
          X          = aux_x_a_m1,
          newX       = df_x_pred, #technically shoudn't be all units but indicator in IPW will make unwanted values 0 
          family     = gaussian(),
          obsWeights = aux_x_a_m1.w$weight,
          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
          control    = list(saveFitLibrary = TRUE),
          cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
        )
        
        df$y110_fit.x <- predict.SuperLearner(y110_sl.x, newdata = df_x_pred)$pred
        
        
        main_list[[k]] <- df[cf_fold[[k]], ]
      }
      
      main_df <- reduce(main_list, bind_rows)
      
      #############################
      # Construct rEIFs
      #############################
      # main_df <- df
      main_df[is.na(main_df)] <- 0
      
      main_df <- main_df %>% mutate(
        w0_y = if_else(hs == 1, 1/a_fit, 0) %>% trimQ(0.01, 0.99),
        w0_n = if_else(hs == 0, 1/(1 - a_fit), 0) %>% trimQ(0.01, 0.99),
        w1_yy = if_else(college == 1, w0_y/m1_fit, 0) %>% trimQ(0.01, 0.99),
        w1_yn = if_else(college == 0, w0_y/(1-m1_fit), 0) %>% trimQ(0.01, 0.99),
        w2_yyy = if_else(ba == 1, w1_yy/m2_fit, 0) %>% trimQ(0.01, 0.99),
        w2_yyn = if_else(ba == 0, w1_yy/(1 - m2_fit), 0) %>% trimQ(0.01, 0.99),
        w3_yyyn = if_else(gschool == 0, w2_yyy/(1 - m3_fit), 0) %>% trimQ(0.01, 0.99),
        w3_yyyy = if_else(gschool == 1, w2_yyy/m3_fit, 0) %>% trimQ(0.01, 0.99),
      )
      
      
      main_df <- main_df %>%
        mutate(
          
          # M1(1), M2(1,1), M3(1,1,1)
          m1.ay_eif = m1_fit +  w0_y * (college - m1_fit),
          m1.aym2y_eif = m2_fit_noz + w1_yy * (ba - m2_fit_noz),
          m1.aym2ym3y_eif = m3_fit.x + w1_yy * (m3_fit - m3_fit.x) + w2_yyy * (gschool - m3_fit),
          
          # Y(a)
          y1_eif = y1_fit +  w0_y * (earnings - y1_fit),
          y0_eif = y0_fit +  w0_n * (earnings - y0_fit),
          
          # Y(a,m1)
          y11_eif = y11_fit +  w1_yy * (earnings - y11_fit),
          y10_eif = y10_fit +  w1_yn * (earnings - y10_fit),
          
          # Y(1,1,m2)
          y111_eif = y111_fit.x +  w1_yy * (y111_fit - y111_fit.x) + w2_yyy * (earnings - y111_fit),
          y110_eif = y110_fit.x +  w1_yy * (y110_fit - y110_fit.x) + w2_yyn * (earnings - y110_fit),
          
          # Y(a,1,1,m3) 
          y1111_eif = y1111_fit.x +  w1_yy * (y1111_fit - y1111_fit.x) + w3_yyyy * (earnings - y1111_fit),
          y1110_eif = y1110_fit.x +  w1_yy * (y1110_fit - y1110_fit.x) + w3_yyyn * (earnings - y1110_fit),
          
          # Delta terms
          delta0_eif = y10_eif - y0_eif,
          delta1_eif = y110_eif - y10_eif,
          delta2_eif = y1110_eif- y110_eif,
          delta3_eif = y1111_eif - y1110_eif,
          
          # Tau terms
          tau0_eif = y1_eif - y0_eif,
          tau1_eif = y11_eif - y10_eif,
          tau2_eif = y111_eif - y110_eif,
          tau3_eif = delta3_eif,
          
          # Pi terms
          pi1_eif = m1.ay_eif,
          pi2_eif = m1.aym2y_eif,
          pi3_eif = m1.aym2ym3y_eif,
          
          # Eta terms
          eta1_eif = tau0_eif - delta0_eif - pi1_eif*weighted.mean(tau1_eif, weight) - tau1_eif*weighted.mean(pi1_eif, weight) + weighted.mean(tau1_eif, weight)*weighted.mean(pi1_eif, weight),
          eta2_eif = tau1_eif - delta1_eif - pi2_eif*weighted.mean(tau2_eif, weight) - tau2_eif*weighted.mean(pi2_eif, weight) + weighted.mean(tau2_eif, weight)*weighted.mean(pi2_eif, weight),
          eta3_eif = tau2_eif - delta2_eif - pi3_eif*weighted.mean(tau3_eif, weight) - tau3_eif*weighted.mean(pi3_eif, weight) + weighted.mean(tau3_eif, weight)*weighted.mean(pi3_eif, weight),
          
          
          # Direct and continuation effects (theta)
          theta0_eif = delta0_eif,
          theta1_eif = weighted.mean(delta1_eif, weight)*pi1_eif + weighted.mean(pi1_eif, weight)*delta1_eif + eta1_eif - weighted.mean(delta1_eif, weight)*weighted.mean(pi1_eif, weight),
          
          theta2_eif = weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*delta2_eif +
            weighted.mean(pi1_eif, weight)*weighted.mean(delta2_eif, weight)*pi2_eif +
            weighted.mean(pi2_eif, weight)*weighted.mean(delta2_eif, weight)*pi1_eif +
            weighted.mean(pi1_eif, weight)*eta2_eif + weighted.mean(eta2_eif, weight)*pi1_eif -
            2*(weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*weighted.mean(delta2_eif, weight)) - weighted.mean(eta2_eif, weight)*weighted.mean(pi1_eif, weight),
          
          theta3_eif = weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*weighted.mean(pi3_eif, weight)*delta3_eif +
            weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*weighted.mean(delta3_eif, weight)*pi3_eif +
            weighted.mean(pi3_eif, weight)*weighted.mean(pi1_eif, weight)*weighted.mean(delta3_eif, weight)*pi2_eif +
            weighted.mean(pi3_eif, weight)*weighted.mean(pi2_eif, weight)*weighted.mean(delta3_eif, weight)*pi1_eif +
            
            weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*eta3_eif +
            weighted.mean(eta3_eif, weight)*weighted.mean(pi2_eif, weight)*pi1_eif +
            weighted.mean(eta3_eif, weight)*weighted.mean(pi1_eif, weight)*pi2_eif -
            3*(weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*weighted.mean(pi3_eif, weight)*weighted.mean(delta3_eif, weight)) - 
            2*(weighted.mean(pi1_eif, weight)*weighted.mean(pi2_eif, weight)*weighted.mean(eta3_eif, weight)),
        )
      # 1.15 +  0.02977888 +  0.1905948 +  0.002171962 = 1.38   
      # -0.00288*0.424 +   0.0310 =  0.02977888
      # 0.424*0.268*0.737 + 0.424*0.252  =  0.1905948
      # 0.424*0.268*0.286*0.249 +  0.424*0.268*-0.0521 = 0.002171962
      
      out[[i]] <- main_df
      
      
      #################################################
      # rEIF regressions
      #################################################
      
      design <- svydesign(ids = ~ 1, weights = ~ weight, data = main_df)
      
      depvars <- select(main_df, starts_with("theta"), starts_with("pi"), starts_with("tau"), starts_with("delta"), starts_with("eta")) %>% names()
      
      overall_mods <- map(depvars, ~ svyglm(as.formula(expr(!!sym(.x) ~ 1)), design = design))
      
      ############################
      # Output
      ############################
      
      overall_out[[i]] <- overall_mods %>%
        set_names(str_sub(depvars, end = -5)) %>%
        enframe(name = "estimand", value = "model") %>%
        mutate(est = map(model, coef),
               vcov = map(model, vcov),
               se = map(vcov, ~ sqrt(diag(.x)))) %>%
        select(estimand, est, se) %>%
        unnest(c(est, se))
      
    }
    
    overall_df <- overall_out %>%
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
      select(-type)
    
    name <- paste0("01_dml_", earn, ".RDS")
    saveRDS(overall_df, file = name)
    
  }
}



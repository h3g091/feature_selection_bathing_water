source("read_and_lag_preprocessed_data.R")
data_kbw<-preprocess_river_data_from_path(data_path_kleine_badewiese, "kleine badewiese", kleine_BW = T)
data_isar<-preprocess_river_data_from_path(data_path_isar, "Isar")
data_ilz<-preprocess_river_data_from_path(data_path_ilz, "Ilz")
data_rhein<-preprocess_river_data_from_path(data_path_rhein, "Rhein")
data_mosel<-preprocess_river_data_from_path(data_path_mosel, "Mosel")
data_ruhr<-preprocess_river_data_from_path(data_path_ruhr, "Ruhr")

build_models<-function(data_path, river_name, kleine_BW = F){
  
  #data_path<-data_path_ilz
  #river_name<- "Ilz"
  #kleine_BW = F
  source("read_and_lag_preprocessed_data.R")
  river_data<-preprocess_river_data_from_path(data_path, river_name, kleine_BW)
  #rm(preprocess_river_data_from_path)
  
  seed <- 123
  
  
  #split data for fs
  
  set.seed(seed)
  training.samples <- river_data$log_e.coli %>%
    createDataPartition(p = 0.8, list = FALSE)
  train_data  <- river_data[training.samples, ]
  test_data <- river_data[-training.samples, ]
  
  
  
  #scaling right before every fs 
  ######## not neccessary
  train_data_full_scaled<-scale(train_data)
  scale_attributes_full <- attributes(train_data_full_scaled)
  train_data_full_scaled<-as.data.frame(train_data_full_scaled)
  
  ### in ressource function
  scale_test_data_with_attributes<-function(data,attributes_from_scale){
    #  attributes_from_scale<-scale_attributes_full
    #data<-train_data_full
    a<-sweep(data,2, attributes_from_scale$`scaled:center`,FUN = "-" )
    data_scaled<-sweep(a,2, attributes_from_scale$`scaled:scale`,FUN = "/" )
    return(data_scaled)
  }
  rescale_data_with_attributes<-function(data,attributes_from_scale){
    #  attributes_from_scale<-scale_attributes_full
    #data<-train_data_full
    a<-sweep(data,2, attributes_from_scale$`scaled:scale`,FUN = "*" )
    data_scaled<-sweep(a,2, attributes_from_scale$`scaled:center`,FUN = "+" )
    return(data_scaled)
  }
  
  #holdout test_set scaled
  test_data_full_scaled<-scale_test_data_with_attributes(test_data,scale_attributes_full)
  
  
  #regularized regression model building 
  
  #build 5-fold-cross with foldid
  n = nrow(train_data_full_scaled)
  foldid=sample(rep(seq(5),length=n))
  
  # Predictor variables
  train_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, train_data_full_scaled)[,-1]
  # Outcome variable
  log_e.coli_true <- train_data$log_e.coli
  #enable parallel computing
  require(doMC)
  registerDoMC(cores=4)
  
  cv_lasso <- cv.glmnet(train_data_matrix_with_interactions, log_e.coli_true, alpha = 1, parallel = T, foldid = foldid)
  
  # Fit the final model on the training data
  #for prediction R2 etc.
  lasso_min_model <- glmnet(train_data_matrix_with_interactions, log_e.coli_true,
                            alpha = 1, lambda = cv_lasso$lambda.min)
  lasso_1se_model <- glmnet(train_data_matrix_with_interactions, log_e.coli_true,
                            alpha = 1, lambda = cv_lasso$lambda.1se)
  
  extract_coefficients_from_cv.glmnet_model<- function(cv.glmnet.fit, lambda_value = "lambda.min"){
    tmp_coeffs <- coef(cv.glmnet.fit, s = lambda_value)
    return(data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)%>% arrange(desc(abs(coefficient))))
    
  }
  
  features_lasso_lambda_min<- extract_coefficients_from_cv.glmnet_model(cv_lasso)
  features_lasso_lambda_1se<- extract_coefficients_from_cv.glmnet_model(cv_lasso, "lambda.1se")
  
  
  
  #stepwise model building
  source("FS_algo_functions.R")
  step_5_aic_models_full_dataset<- get_list_of_step_aic_models(train_data_full_scaled, state_of_the_art = T)
  step_aic_models_full_dataset<- get_list_of_step_aic_models(train_data_full_scaled, state_of_the_art = F)
  
  step_5_bic_models_full_dataset <- get_list_of_step_bic_models(train_data_full_scaled, state_of_the_art = T)
  step_bic_models_full_dataset <- get_list_of_step_bic_models(train_data_full_scaled, state_of_the_art = F)
  
  # rm(get_list_of_step_bic_models)
  #  rm(get_list_of_step_aic_models)
  
  #######for identification  not neccessary right now
  
  
  
  #split data for fs
  
  
  
  features_step_5_aic_full<-as.character(formula(step_5_aic_models_full_dataset))[[3]]
  features_step_aic_full<-as.character(formula(step_aic_models_full_dataset))[[3]]
  features_step_5_bic_full<-as.character(formula(step_5_bic_models_full_dataset))[[3]]
  features_step_bic_full<-as.character(formula(step_bic_models_full_dataset))[[3]]
  
  return(list(step_5_aic_models_full_dataset,step_aic_models_full_dataset, step_5_bic_models_full_dataset,step_bic_models_full_dataset ,features_lasso_lambda_min, features_lasso_lambda_1se))
}

models_kbw<-build_models(data_path_kleine_badewiese, "kleine badewiese", kleine_BW = T)
models_isar<-build_models(data_path_isar, "Isar")
models_ilz<-build_models(data_path_ilz, "Ilz")
models_rhein<-build_models(data_path_rhein, "Rhein")
models_mosel<-build_models(data_path_mosel, "Mosel")
models_ruhr<-build_models(data_path_ruhr, "Ruhr")
##retrain with whole dataset to check significance
extract_features_from_models <- function(models_river){
  #models_river<- models_isar
  features_temp<- lapply(models_river[1:4], coef)  
  features_temp<-lapply(features_temp, names)
  
  features_temp<-lapply(features_temp, function(x) paste(x[-1], collapse = " +"))
  features_temp<- append(features_temp,paste(models_river[[5]]$name[-1],collapse = " + " ))
  features_temp<- append(features_temp,paste(models_river[[6]]$name[-1],collapse = " + " ))
  return(features_temp)
}

features_kbw <-  extract_features_from_models(models_kbw)
features_isar <- extract_features_from_models(models_isar)
features_ilz <- extract_features_from_models(models_ilz)
features_mosel <- extract_features_from_models(models_mosel)
features_rhein <-  extract_features_from_models(models_rhein)
features_ruhr <-  extract_features_from_models(models_ruhr)


build_final_models_with_full_dataset<- function(features_river, data_full){
  
  #data_full<- data_kbw
  #features_river <- features_kbw[[1]]
  return(lm(paste("log_e.coli ~ ",features_river), as.data.frame(scale(as.data.frame(data_full)))))
  
}

list_final_models_kbw<-lapply(features_kbw,build_final_models_with_full_dataset, data_kbw )
build_coefficient_oveview <- function(trained_final_model){
  #trained_final_model<-list_final_models_kbw[[1]]
  df_temp<-as.data.frame(coef(trained_final_model))
  df_temp<-cbind(names(coef(trained_final_model)),df_temp)
  df_temp<-cbind(df_temp,summary(trained_final_model)$coef[,4] <= .05)
  names(df_temp) <- c("features", "coef", "significant")
  return(df_temp)
}
coefficients_final_model<-lapply(list_final_models_kbw,build_coefficient_oveview)
coefficients_final_model_df<-as.data.frame( do.call(rbind,coefficients_final_model))

#write.csv(coefficients_final_model_df, file = "/Users/heiko.langer/Masterarbeit_lokal/Masterprojekt/automated_method/last_try_to_fit_timos_expectation/results/kbw_details.csv")

writexl::write_xlsx(coefficients_final_model_df, "/Users/heiko.langer/Masterarbeit_lokal/Masterprojekt/automated_method/last_try_to_fit_timos_expectation/results/kbw_details.xlsx")

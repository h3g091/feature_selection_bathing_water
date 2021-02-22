
library(roxygen2)
library(glmnet)
library(fhpredict)
library(tidyverse)
#important
library(kwb.flusshygiene)
library(dplyr)
library(Matrix)
library(caret)

#automated method
# add path of data. Data  
#data_path_isar<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Bayern/Isar/DATA_preprocessed_csv")
data_path_isar<-list(havel = "/Data/Daten_Bayern/Isar/DATA_preprocessed_csv")
data_path_ilz<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Bayern/Ilz/DATA_preprocessed_csv")
data_path_mosel<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Rhein_Mosel_Lahn/Mosel/DATA_preprocessed_csv")
data_path_rhein<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Rhein_Mosel_Lahn/Rhein/DATA_preprocessed_csv")
data_path_ruhr <- list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Ruhr/Ruhr/DATA_preprocessed_csv")
data_path_kleine_badewiese<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/kleine_Badewiese/Havel/DATA_preprocessed_csv")





#' Builds tibble that contains the results of the validation and accuracy of the build models for a specific dataset.
#' 
#' Preprocesses data from river_path and subsequently builds two lasso and four stepwise models and uses percentage-coverage method to valdiate them
#' 
#' @param path string in with the data_preprocesses_csv files directory
#' @param river_name name of the bathing site or river name as it is shown in tibble
#' @param kleine_BW TRUE/FALSE value if Data from kleine Badewiese is used. Some inconsistencies in the data haave led to extra needed steps to treat the data the same
#' @return return tibble of results
build_river_tibble<-function(data_path, river_name, kleine_BW = F){
  
  #data_path<-data_path_mosel
  #river_name<- "Isar"
  #kleine_BW = F
  ##preprocesses river data and 
  source("read_and_lag_preprocessed_data.R")
  river_data<-preprocess_river_data_from_path(data_path, river_name, kleine_BW)
  rm(preprocess_river_data_from_path)
  # fix pseudo random seed to get always same results
  seed <- 123
  
  #split data for feature selection in 90% Train-set and 20% Test-set
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
  
  ###for preprocessing data "whitening of features" 
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
  
  #scale test-set with same attributes than Train-set 
  test_data_full_scaled<-scale_test_data_with_attributes(test_data,scale_attributes_full)
  
  #regularized regression model building 
  
  #build split of five-fold-cross validation with foldid
  n = nrow(train_data_full_scaled)
  foldid=sample(rep(seq(5),length=n))
  
  # Predictor variables with all interactions possible
  train_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, train_data_full_scaled)[,-1]
  # Outcome variable
  log_e.coli_true <- train_data$log_e.coli
  #enable parallel computing
  require(doMC)
  registerDoMC(cores=4)
  # build lasso models with fixed folds for feature selection to make the results reproducable
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
  
  features_step_5_aic_full<-as.character(formula(step_5_aic_models_full_dataset))[[3]]
  features_step_aic_full<-as.character(formula(step_aic_models_full_dataset))[[3]]
  features_step_5_bic_full<-as.character(formula(step_5_bic_models_full_dataset))[[3]]
  features_step_bic_full<-as.character(formula(step_bic_models_full_dataset))[[3]]
  

  #########################
  #predict on train to calculate sd
  train_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, train_data_full_scaled)[,-1]
  
  train_predictions_lasso_min <- lasso_min_model %>% predict(train_data_matrix_with_interactions) %>% as.vector()
  train_predictions_lasso_1se <- lasso_1se_model %>% predict(train_data_matrix_with_interactions) %>% as.vector()
  
  train_prediction_step_5_aic<-predict(step_5_aic_models_full_dataset, train_data_full_scaled)%>% as.vector()
  train_prediction_step_aic<-predict(step_aic_models_full_dataset, train_data_full_scaled)%>% as.vector()
  train_prediction_step_5_bic<-predict(step_5_bic_models_full_dataset, train_data_full_scaled)%>% as.vector()
  train_prediction_step_bic<-predict(step_bic_models_full_dataset, train_data_full_scaled)%>% as.vector()
  
  
  
  #predict on test
  test_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, test_data_full_scaled)[,-1]
  
  predictions_lasso_min <- lasso_min_model %>% predict(test_data_matrix_with_interactions) %>% as.vector()
  predictions_lasso_1se <- lasso_1se_model %>% predict(test_data_matrix_with_interactions) %>% as.vector()
  
  prediction_step_5_aic<-predict(step_5_aic_models_full_dataset, test_data_full_scaled)%>% as.vector()
  prediction_step_aic<-predict(step_aic_models_full_dataset, test_data_full_scaled)%>% as.vector()
  prediction_step_5_bic<-predict(step_5_bic_models_full_dataset, test_data_full_scaled)%>% as.vector()
  prediction_step_bic<-predict(step_bic_models_full_dataset, test_data_full_scaled)%>% as.vector()
  
  #' calculate model performance metric R2
  #' @param prediction A vector of predicted values
  #' @param true A vector of the corresponding true values to prediction
  #' @param return returns R2, A number 
  get_R2_model_performance_metrics <- function(prediction, true){
    cor(prediction, true)^2
  }
  
  
  Train_set_R2<-c(get_R2_model_performance_metrics(train_prediction_step_5_aic, train_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(train_prediction_step_aic, train_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(train_prediction_step_5_bic, train_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(train_prediction_step_bic, train_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(train_predictions_lasso_min, train_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(train_predictions_lasso_1se, train_data_full_scaled$log_e.coli))
  
  Test_set_R2<-c(get_R2_model_performance_metrics(prediction_step_5_aic, test_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(prediction_step_aic, test_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(prediction_step_5_bic, test_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(prediction_step_bic, test_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(predictions_lasso_min, test_data_full_scaled$log_e.coli),
                 get_R2_model_performance_metrics(predictions_lasso_1se, test_data_full_scaled$log_e.coli))
  #combine all linear models in one list for validation
  
  #validation of data
  
  
  
  
  #split data for validation
  source("get_split_for_validation.R")
  data_split_for_validation<-get_split_for_validation(river_data, foldid) 
  #rm(get_split_for_validation)
  
  #' Extracts features from cv.glmnet object at given lambda value
  #' 
  #' @param lambda_value A number or "lambda.1se" or "lambda.min" regurlarization value that determines the model and therefore the features in the model
  #' @param cv.glmnet_model A cv.glmnet model 
  #' @param return returns A string, reutrns the features as formula with '+' between each feature
  extract_coefficients_as_formula_from_glmnet<-function(cv.glmnet_model, lambda_value= "lambda.min"){
    coefficients_temp<-coefplot::extract.coef(cv.glmnet_model, lambda = lambda_value)
    return(paste(coefficients_temp$Coefficient[-1], collapse = " + "))
    
  }
  lasso.min_formula<-extract_coefficients_as_formula_from_glmnet(cv_lasso, lambda_value = "lambda.min")
  lasso.1se_formula<-extract_coefficients_as_formula_from_glmnet( cv_lasso, lambda_value = "lambda.1se")
  
  features <- c(features_step_5_aic_full, 
                features_step_aic_full,
                features_step_5_bic_full, 
                features_step_bic_full,
                lasso.min_formula,
                lasso.1se_formula
                )
  source("validation_methods.R")
  ##calculate validations for every lm and count how often it could be validated
  list_linear_models <- list()
  #rm(list_linear_models)
  
  validation_df_step_aic <- validation_5_fold(features_step_aic_full, data_split_for_validation[[1]], data_split_for_validation[[2]])
  validation_df_step_5_aic <- validation_5_fold(features_step_5_aic_full, data_split_for_validation[[1]], data_split_for_validation[[2]])
  
  validation_df_step_bic <- validation_5_fold(features_step_bic_full, data_split_for_validation[[1]], data_split_for_validation[[2]])
  validation_df_step_5_bic <- validation_5_fold(features_step_5_bic_full, data_split_for_validation[[1]], data_split_for_validation[[2]])
  
  validation_df_lasso.min<-validation_5_fold(lasso.min_formula, data_split_for_validation[[1]], data_split_for_validation[[2]])
  validation_df_lasso.1se<-validation_5_fold(lasso.1se_formula, data_split_for_validation[[1]], data_split_for_validation[[2]])
  
  check_if_5_fold_validation_pass<-function(validation_df){
    return(validation_df$in95 ==5 & validation_df$below95 ==5 & validation_df$below90 ==5 & validation_df$in50 ==5)
  }
  

  five_fold_validated<-c(
    validation_df_step_5_aic$passed_folds,
    validation_df_step_aic$passed_folds,
    validation_df_step_5_bic$passed_folds,
    validation_df_step_bic$passed_folds,
    validation_df_lasso.min$passed_folds,
    validation_df_lasso.1se$passed_folds)
  validation_on_test<-c(
    validation_on_test_stepwise(model = step_5_aic_models_full_dataset, train_data_full_scaled, test_data_full_scaled),
    validation_on_test_stepwise(model = step_aic_models_full_dataset, train_data_full_scaled, test_data_full_scaled),
    validation_on_test_stepwise(model = step_5_bic_models_full_dataset, train_data_full_scaled, test_data_full_scaled),
    validation_on_test_stepwise(model = step_bic_models_full_dataset, train_data_full_scaled, test_data_full_scaled),
    validation_on_test_glmnet(model = lasso_min_model, train_data_matrix_with_interactions, test_data_matrix_with_interactions, train_data_full_scaled,test_data_full_scaled),
    validation_on_test_glmnet(model = lasso_1se_model, train_data_matrix_with_interactions, test_data_matrix_with_interactions, train_data_full_scaled,test_data_full_scaled)
  )
  
 # rm(validation_5_fold)
 # rm(test_beta)
#  rm(get_validation_per_unique_formula_over_all_validion_iterations)
  
  algo_names <- c("step_5_aic","step_aic","step_5_bic", "step_bic","lasso_min","lasso_1se")
  
  
  get_min_criterion_value_of_step_wise_model <- function(stepwise_model){
    return(stepwise_model$anova$AIC[nrow(stepwise_model$anova)])
  }
  
  
  
  regularization_parameter <- c( "AIC","AIC", "BIC", "BIC","Lambda", "Lambda")
  
  regularization_parameter_value <- c(get_min_criterion_value_of_step_wise_model(step_5_aic_models_full_dataset),
                                get_min_criterion_value_of_step_wise_model(step_aic_models_full_dataset),
                                get_min_criterion_value_of_step_wise_model(step_5_bic_models_full_dataset),
                                get_min_criterion_value_of_step_wise_model(step_bic_models_full_dataset),
                                cv_lasso$lambda.min, 
                                cv_lasso$lambda.1se
  )
  
  
  
  
  river_tibble<-tibble(rep(river_name, 6),algo_names, regularization_parameter, regularization_parameter_value,  Train_set_R2, Test_set_R2,five_fold_validated, validation_on_test,features)
  return(river_tibble)
}


# build tibble for each dataset
tibble_kleine_Badewiese<-build_river_tibble(data_path_kleine_badewiese, "kleine badewiese", kleine_BW = T)
tibble_isar<-build_river_tibble(data_path_isar, "Isar")
tibble_ilz<-build_river_tibble(data_path_ilz, "Ilz")
tibble_rhein<-build_river_tibble(data_path_rhein, "Rhein")
tibble_mosel<-build_river_tibble(data_path_mosel, "Mosel")
tibble_ruhr<-build_river_tibble(data_path_ruhr, "Ruhr")

big_tibble<- bind_rows(tibble_kleine_Badewiese, 
                    tibble_ilz, 
                    tibble_isar, 
                    tibble_mosel,
                    tibble_rhein,
                    tibble_ruhr)
names(big_tibble)<- c("River","Method","para","para value", "Train R2", "Test R2", "5-fold valid","test valid", "features"       )
big_tibble$`Train R2`<-round(big_tibble$`Train R2`,2)
big_tibble$`Test R2`<-round(big_tibble$`Test R2`,2)

big_tibble$n_features <- unlist(lapply(big_tibble$features, str_count, "\\+"))+1
big_tibble<-big_tibble[,c(1:(ncol(big_tibble)-2),10,9)]
#big_tibble$n_features
print(big_tibble, n=40)
overwrite_big_tibble <- F
if(overwrite_big_tibble){
  
saveRDS(big_tibble,"results/big_tibble/big_tibble.rds")
}


writexl::write_xlsx(big_tibble,"/Users/heiko.langer/Masterarbeit_lokal/Masterprojekt/automated_method/last_try_to_fit_timos_expectation/results/new_results.xlsx")
#write.csv(big_tibble, file = "/Users/heiko.langer/Masterarbeit_lokal/Masterprojekt/automated_method/last_try_to_fit_timos_expectation/results/new_results.csv")
#take only that are at least 10 biggest valdiations!!



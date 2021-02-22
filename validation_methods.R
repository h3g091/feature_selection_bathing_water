#' Performs two sided beta test between list of TRUE and FALSE values and corresponding percentile that shall be tested
#' Beta test confidence niveau of 5%
#' @param is_true A list of T and F values. 
#' @param percentile A number/double, represents the percentile to test. p.ex. 0.95 to test below95 and within95 uncertainty interval
#' @param return returns logical, T or F if test was passed or not
{
  test_beta <- function(is_true, percentile)  {
    u<- F
    if(stats::pbeta(
      q = percentile,
      shape1 = sum(is_true) + 1,
      shape2 = sum(! is_true) + 1
    ) > 0.05 & stats::pbeta(
      q = percentile,
      shape1 = sum(is_true) + 1,
      shape2 = sum(! is_true) + 1
    ) < 0.95  ){
      u<- T
    }
    return(u)    
    
  }
}

#' calculate model performance metric R2
#' @param prediction A vector of predicted values
#' @param true A vector of the corresponding true values to prediction
#' @param return returns R2, A number 
rsq <- function (x, y) cor(x, y) ^ 2


#' checks on how many folds a model with specific formula passes percentage coverage criteria in five fold cross validation on train_set
#' @param formel A string with feature names pased together with "+". "log_E.coli ~ " is missing in formel and is pasted to it
#' @param data_train A list of dataframes with length five. Each dataframe contains the training datapoints for the specific fold 
#' @param data_test A list of dataframes with length five. Each dataframe contains the test datapoints for the corresponding fold of data_train
#' @param return returns A dataframe, dataframe contains different insights. Most important, the number of folds, the modell passed the percentage coverage in the five-fold cross valdiation. 

{     
  validation_5_fold <- function(formel,data_train,data_test){      
    #formel <- features_step_aic_full
    #data_train <- data_split_for_validation[[1]]
    #data_test <- data_split_for_validation[[2]]
    #data_train <-validation_train_data_random_split[[1]]
    #data_test <-validation_test_data_random_split[[1]]
    
    dummy_df <- data.frame(0,0,0,0,0,0,0,0,0)
    names(dummy_df) <- c("mean_MSE", "mean_R2", "n_obs", "mean_sd_pred_train_to_true","in95", "below95", "below90", "in50", "passed_folds")
    
    n_features<-str_count(string = formel, pattern = "\\+")+1
    
    list_test_error <- list()
    list_sd_pred_train_to_true <- list()
    list_R2_iteration <- list()
    
    ##### each of the folds 1 to 5 
    for(fold_indx in 1:5){
      
      data_train_fold <- data_train[[fold_indx]]
      data_test_fold  <- data_test[[fold_indx]]
      
      train_bacteria  <- data_train_fold$log_e.coli
      test_bacteria <- data_test_fold$log_e.coli
      full_with_interaction <- formula(paste("log_e.coli ~ ",formel,sep = ""))
      linear_model<- lm(formula = full_with_interaction, data = data_train_fold )
      
      
      n_obs<- nrow(data_train_fold)
      
      #prediction
      pred_train <- predict(linear_model,newdata = data_train_fold)
      pred_test <- predict(linear_model, newdata =  data_test_fold)
      
      
      
      r2<- rsq( train_bacteria, pred_train)
      list_R2_iteration <-append(list_R2_iteration,list(r2))
      sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/(length(train_bacteria)-1))
      
      
      percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
      
      
      below95 = test_bacteria < percentile_train_95
      
      below90 = test_bacteria < percentile_train_90
      
      within95 = test_bacteria < percentile_train_97_5 & test_bacteria > percentile_train_2_5
      
      within50 = test_bacteria < percentile_train_75 & test_bacteria > percentile_train_25
      
      
      within95_fold_result<-test_beta(is_true = within95, percentile = .95 )
      within50_fold_result<-test_beta(is_true = within50, percentile = .5)
      below95_fold_result<-test_beta(is_true = below95, percentile = .95 )
      below90_fold_result<-test_beta(is_true = below90, percentile = .90 )
      
      
      dummy_df$in95 <- dummy_df$in95 + within95_fold_result
      dummy_df$below95 <- dummy_df$below95 + below95_fold_result
      dummy_df$below90 <- dummy_df$below90 + below90_fold_result
      dummy_df$in50 <- dummy_df$in50+within50_fold_result
      
      if(within95_fold_result & below95_fold_result & below90_fold_result& within50_fold_result){
        dummy_df$passed_folds <- dummy_df$passed_folds+ 1
      }
      
      list_test_error <- append(list_test_error,list(mean((pred_test- data_test_fold$log_e.coli)^2)))
      list_sd_pred_train_to_true <-append(list_sd_pred_train_to_true, list(sd_pred_train_to_true))        
      
    }
    
    dummy_df$n_obs<-nrow(data_train_fold)+nrow(data_test_fold)
    dummy_df$mean_MSE<-mean(unlist(list_test_error))
    dummy_df$mean_sd_pred_train_to_true<-mean(unlist(list_sd_pred_train_to_true))
    dummy_df$mean_R2 <-mean(unlist(list_R2_iteration))
    
    return(dummy_df)
  }
}



#' Tests if a already trained glmnet model passes the validation on the Test-set
#' @param model A cv.glmnet object with fix lambda value
#' @param train_data_matrix_with_interactions A matrix, can be sparse. Contains true Log_e.coli value in train-set and all features with interactions
#' @param test_data_matrix_with_interactions A matrix, can be sparse. Contains true Log_e.coli value in test-set and all features with interactions
#' @param return A logical, if the model passes the percentage-coverage criteria on test-set
{     
  validation_on_test_glmnet <- function(model, train_data_matrix_with_interactions, test_data_matrix_with_interactions, data_train, data_test){      
   
    train_bacteria  <- data_train$log_e.coli
    test_bacteria <- data_test$log_e.coli
    
    #prediction
    
    pred_train <- model %>% predict(train_data_matrix_with_interactions) %>% as.vector()
    pred_test <- model %>% predict(test_data_matrix_with_interactions) %>% as.vector()


    #get sd from difference between train_bacteria_true to the predicted! mean is the prediction
    sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/length(train_bacteria))
    
    #dummy_df$sd_pred_train_to_true <- sd_pred_train_to_true
    
    percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
    
    below95 = test_bacteria < percentile_train_95
    below90 = test_bacteria < percentile_train_90
    within95 = test_bacteria < percentile_train_97_5 & test_bacteria > percentile_train_2_5
    within50 = test_bacteria < percentile_train_75& test_bacteria > percentile_train_25
    
    
    
    
 
    return(test_beta(is_true = within95, percentile = .95 ) &
             test_beta(is_true = below95, percentile = .95 )&
             test_beta(is_true = below90, percentile = .90 )&
             test_beta(is_true = within50, percentile = .5)
    )
  }
}
#' Tests if a already trained glmnet model passes the validation on the Test-set
#' @param model A linear model object already trained on full Train-set
#' @param data_train A dataframe Contains true Log_e.coli value in train-set and all features with interactions
#' @param data_test A dataframe. Contains true Log_e.coli value in test-set and all features with interactions
#' @param return A logical, if the model passes the percentage-coverage criteria on test-set
{     
  validation_on_test_stepwise <- function(model, data_train, data_test){      
     #model <- step_5_aic_models_full_dataset 
     # data_train<- train_data_full_scaled
    #  data_test<- test_data_full_scaled
    train_bacteria  <- data_train$log_e.coli
    test_bacteria <- data_test$log_e.coli
    
    #prediction
  
    pred_train <- model %>% predict(data_train) %>% as.vector()
    pred_test <- model %>% predict(data_test) %>% as.vector()
    
    
    #get sd from difference between train_bacteria_true to the predicted! mean is the prediction
    sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/length(train_bacteria))
    
    percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
    
  
    below95 = test_bacteria < percentile_train_95
    below90 = test_bacteria < percentile_train_90
    within95 = test_bacteria < percentile_train_97_5 & test_bacteria > percentile_train_2_5
    within50 = test_bacteria < percentile_train_75& test_bacteria > percentile_train_25

    
    return(test_beta(is_true = within95, percentile = .95 ) &
             test_beta(is_true = below95, percentile = .95 )&
             test_beta(is_true = below90, percentile = .90 )&
             test_beta(is_true = within50, percentile = .5)
    )
  }
}
#' Checks if all percentage coverage criteria == 5 in the valdiation dataframe of a specific model and it passes the five-fold cross validation
#' @param validation_df A datafram. Contains the number of folds a model passes each criteria in five-fold cross validation
#' @param return A logical, if all percentage coverage criteria are 5 --> T
check_if_all_folds_validated <-function(validation_df){
  
  if(all(validation_df$in95==5,validation_df$below95==5,validation_df$below90==5,validation_df$in50==5)){
    return(T)
  }else if(!all(validation_df$in95==5,validation_df$below95==5,validation_df$below90==5,validation_df$in50==5)){
    return(F)
  }
}

#' Calculates the number of times a model passes the percentage-coverage criteria in the five-fold cross validation if used on different splits of the same data.
#' can take several models as formulas in a list and calculates validation df for each of the formulas in every of the different splits
#' @param unique_river_algo_formulas A list of formulas, each formula contains only feature names pasted together with "+"
#' @param all_train_iterations_validation_scaled_river A list  of lsits of dataframes. Prescalled!!  Outer list contains the different splits of the same data with n datapoints. each of these lists contains five different dataframes for each of the five-folds each containing the datapoints for training in a specific split
#' @param all_test_iterations_validation_scaled_river A list  of lists of dataframes. all_train_iterations_validation_scaled_river
#' @param return A list object with two lists of length m, the number of different iteratoions. - p.ex. 50 if 50 different splits are used. The first list contains the number of iterations a model passed all five folds in the five fold cross validation. The second one contains each of the validation dataframes for each of the iterations
get_validation_per_unique_formula_over_all_validion_iterations<-function(unique_river_algo_formulas, all_train_iterations_validation_scaled_river, all_test_iterations_validation_scaled_river){
 
  list_number_validated<- list()
  list_validation_dfs_lasso_havel <- list()
  for(formula_indx in 1: length(unique_river_algo_formulas)){
    #for(formula_indx in 1: 10){
    #formula_indx<-1
    validation_dfs_lasso_havel <- list()
    for(iteration_indx in 1:length(all_test_iterations_validation_scaled_river)){
      #iteration_indx<- 1
      
      validation_dfs_lasso_havel <- append(validation_dfs_lasso_havel, list(
        validation(unique_river_algo_formulas[formula_indx],
                   all_train_iterations_validation_scaled_river[[iteration_indx]], 
                   all_test_iterations_validation_scaled_river[[iteration_indx]])))
    }
    list_validation_dfs_lasso_havel <- append(list_validation_dfs_lasso_havel, list(validation_dfs_lasso_havel))
    number_validated_per_50<-sum(unlist((lapply(validation_dfs_lasso_havel,check_if_all_folds_validated))))
    list_number_validated <- append(list_number_validated, list(number_validated_per_50))
    
    
  }
  return_object <- c(list(list_number_validated), list(list_validation_dfs_lasso_havel))
  
  return(return_object)
}

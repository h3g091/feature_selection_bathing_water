
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

rsq <- function (x, y) cor(x, y) ^ 2
{     
  validation_5_fold <- function(formel,data_train,data_test){      
    
    #formel<-as.character(formula(step_5_aic_models_full_dataset))[3]
    
    
    #data_train <-validation_train_data_random_split[[1]]
    #data_test <-validation_test_data_random_split[[1]]
    
    
    
    dummy_df <- data.frame(0,0,0,0,0,0,0,0,0)
    names(dummy_df) <- c("mean_MSE", "mean_R2", "n_obs", "mean_sd_pred_train_to_true","in95", "below95", "below90", "in50", "Formula")
    dummy_df$Formula <- formel
    n_features<-str_count(string = formel, pattern = "\\+")+1
    
    fold_indx <- 1
    
    #river_stat_tests$in95 <- river_stat_tests$below95 <-river_stat_tests$below90 <- river_stat_tests$in50 <- river_stat_tests$MSE <- 0
    data_train_fold <- data_train[[fold_indx]]
    data_test_fold  <- data_test[[fold_indx]]
    n_obs<- nrow(data_train_fold)
    train_bacteria  <- data_train_fold$log_e.coli
    test_bacteria <- data_test_fold$log_e.coli
    full_with_interaction <- formula(paste("log_e.coli ~ ",formel,sep = ""))
    linear_model<- lm(formula = full_with_interaction, data = data_train_fold )
    #list_linear_models <<- append(list_linear_models, list(linear_model))
    #list_R2_iteration <- summary(linear_model)[["r.squared"]]
    #list_R2_iteration <- summary(linear_model)[["adj.r.squared"]]
    
    #prediction
    pred_train <- predict(linear_model,newdata = data_train_fold)
    pred_test <- predict(linear_model, newdata =  data_test_fold)
    
    r2<- rsq( train_bacteria, pred_train)
    #r2<- 1 - sum((data_train_fold$log_e.coli-pred_train)^2)/sum((data_train_fold$log_e.coli-mean(data_train_fold$log_e.coli))^2)
    
    #rsq <- function (x, y) cor(x, y) ^ 2
    #rsq( train_bacteria, pred_train)
    #e<-cbind(as.data.frame(train_bacteria), as.data.frame(pred_train))
    #e$error <- sqrt((e$train_bacteria - e$pred_train)^2)
    
    
    
    #adj_R2 <- 1-((1-r2)*(n_obs-1))/(n_obs-n_features-1)
    #list_R2_iteration <-adj_R2
    list_R2_iteration <-r2
    list_test_error <- mean((pred_test- data_test_fold$log_e.coli)^2)
    
    #get sd from difference between train_bacteria_true to the predicted! mean is the prediction
    sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/length(train_bacteria))
    list_sd_pred_train_to_true <- sd_pred_train_to_true
    #dummy_df$sd_pred_train_to_true <- sd_pred_train_to_true
    
    percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
    
    #evaluating ther model has to be classified correctly with every single test train split
    #--> here 5 different splits, if all validations correct than everywhere ==5
    
    below95 = test_bacteria < percentile_train_95
    below90 = test_bacteria < percentile_train_90
    within95 = test_bacteria < percentile_train_97_5 & test_bacteria > percentile_train_2_5
    within50 = test_bacteria < percentile_train_75& test_bacteria > percentile_train_25
    
    
    
    
    
    #dummy_df$in95 <- dummy_df$in95 + test_beta(true = sum(within95), false = sum(!within95), percentile = .95 )
    #dummy_df$below95 <- dummy_df$below95+test_beta(true = sum(below95), false = sum(!below95), percentile = .95 )
    #dummy_df$below90 <- dummy_df$below90+test_beta(true = sum(below90), false = sum(!below90), percentile = .90 )
    #dummy_df$in50 <- dummy_df$in50+test_beta(true = sum(within50), false = sum(!within50), percentile = .5)
    
    
    dummy_df$in95 <- dummy_df$in95 + test_beta(is_true = within95, percentile = .95 )
    dummy_df$below95 <- dummy_df$below95 + test_beta(is_true = below95, percentile = .95 )
    dummy_df$below90 <- dummy_df$below90 + test_beta(is_true = below90, percentile = .90 )
    dummy_df$in50 <- dummy_df$in50+test_beta(is_true = within50, percentile = .5)
    
    
    
    for(fold_indx in 2:5){
      #first_fold
      #fold_indx <- 2
      data_train_fold <- data_train[[fold_indx]]
      data_test_fold  <- data_test[[fold_indx]]
      
      train_bacteria  <- data_train_fold$log_e.coli
      test_bacteria <- data_test_fold$log_e.coli
      full_with_interaction <- formula(paste("log_e.coli ~ ",formel,sep = ""))
      linear_model<- lm(formula = full_with_interaction, data = data_train_fold )
      
      
      n_obs<- nrow(data_train_fold)
      
      
      
      
      #r2<- 1 - sum((data_train_fold$log_e.coli-pred_train)^2)/sum((data_train_fold$log_e.coli-mean(data_train_fold$log_e.coli))^2)
      #prediction
      pred_train <- predict(linear_model,newdata = data_train_fold)
      pred_test <- predict(linear_model, newdata =  data_test_fold)
      
      
      #r2<- 1 - sum((data_train_fold$log_e.coli-pred_train)^2)/sum((data_train_fold$log_e.coli-mean(data_train_fold$log_e.coli))^2)
      r2<- rsq( train_bacteria, pred_train)
      
      #adj_R2 <- 1-((1-r2)*(n_obs-1))/(n_obs-n_features-1)
      #list_R2_iteration <-append(list_R2_iteration,list(adj_R2))
      list_R2_iteration <-append(list_R2_iteration,list(r2))
      #get sd from difference between train_bacteria_true to the predicted! mean is the prediction
      sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/(length(train_bacteria)-1))
      
      
      percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
      percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
      
      
      #evaluating ther model has to be classified correctly with every single test train split
      #--> here 5 different splits, if all validations correct than everywhere ==5
      
      below95 = test_bacteria < percentile_train_95
      
      below90 = test_bacteria < percentile_train_90
      
      within95 = test_bacteria < percentile_train_97_5 & test_bacteria > percentile_train_2_5
      
      within50 = test_bacteria < percentile_train_75 & test_bacteria > percentile_train_25
      
      # test_beta <- function(true, false, percentile)
      #  true= 58
      #false = 0  
      #  { if( pbeta(q = percentile, shape1 = true + 1, shape2 = false + 1) > 0.05)
      
      #  {TRUE}
      
      #   else{FALSE}
      
      #  }
      
      
      
      #dummy_df$in95 <- dummy_df$in95 + test_beta(true = sum(within95), false = sum(!within95), percentile = .95 )
      #dummy_df$below95 <- dummy_df$below95 + test_beta(true = sum(below95), false = sum(!below95), percentile = .95 )
      #dummy_df$below90 <- dummy_df$below90 + test_beta(true = sum(below90), false = sum(!below90), percentile = .90 )
      #dummy_df$in50 <- dummy_df$in50 + test_beta(true = sum(within50), false = sum(!within50), percentile = .5)
      
      dummy_df$in95 <- dummy_df$in95 + test_beta(is_true = within95, percentile = .95 )
      dummy_df$below95 <- dummy_df$below95 + test_beta(is_true = below95, percentile = .95 )
      dummy_df$below90 <- dummy_df$below90 + test_beta(is_true = below90, percentile = .90 )
      dummy_df$in50 <- dummy_df$in50 + test_beta(is_true = within50, percentile = .5)
      
      
      
      
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

#validation_on_test_glmnet(model = lasso_1se_model, train_data_matrix_with_interactions, test_data_matrix_with_interactions, train_data_full_scaled,test_data_full_scaled)
{     
  validation_on_test_glmnet <- function(model, train_data_matrix_with_interactions, test_data_matrix_with_interactions, data_train, data_test){      
   # model <- lasso_1se_model 
  #  data_train<- train_data_full_scaled
  #  data_test<- test_data_full_scaled
    
    #train_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, data_train)[,-1]
    #test_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, data_test)[,-1]
    
   
    
    #formel<-as.character(formula(step_5_aic_models_full_dataset))[3]
    
    
    #data_train <-validation_train_data_random_split[[1]]
    #data_test <-validation_test_data_random_split[[1]]
    
    
    
    #dummy_df <- data.frame(0,0,0,0,0,0,0,0,0)
    #names(dummy_df) <- c("mean_MSE", "mean_R2", "n_obs", "mean_sd_pred_train_to_true","in95", "below95", "below90", "in50")
    #fold_indx <- 1
    
    #river_stat_tests$in95 <- river_stat_tests$below95 <-river_stat_tests$below90 <- river_stat_tests$in50 <- river_stat_tests$MSE <- 0
    #data_train_fold <- data_train[[fold_indx]]
    #data_test_fold  <- data_test[[fold_indx]]
    train_bacteria  <- data_train$log_e.coli
    test_bacteria <- data_test$log_e.coli
    
    #prediction
    #pred_train <- predict(linear_model,newdata = data_train_fold)
    #pred_test <- predict(linear_model, newdata =  data_test_fold)
    pred_train <- model %>% predict(train_data_matrix_with_interactions) %>% as.vector()
    pred_test <- model %>% predict(test_data_matrix_with_interactions) %>% as.vector()
    
    
    #r2<- 1 - sum((data_train_fold$log_e.coli-pred_train)^2)/sum((data_train_fold$log_e.coli-mean(data_train_fold$log_e.coli))^2)
    
    #rsq <- function (x, y) cor(x, y) ^ 2
    #rsq( train_bacteria, pred_train)
    #e<-cbind(as.data.frame(train_bacteria), as.data.frame(pred_train))
    #e$error <- sqrt((e$train_bacteria - e$pred_train)^2)
    
    
    
    #adj_R2 <- 1-((1-r2)*(n_obs-1))/(n_obs-n_features-1)
    #list_R2_iteration <-adj_R2
    
    #get sd from difference between train_bacteria_true to the predicted! mean is the prediction
    sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/length(train_bacteria))
    
    #dummy_df$sd_pred_train_to_true <- sd_pred_train_to_true
    
    percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
    
    #evaluating ther model has to be classified correctly with every single test train split
    #--> here 5 different splits, if all validations correct than everywhere ==5
    
    below95 = test_bacteria < percentile_train_95
    below90 = test_bacteria < percentile_train_90
    within95 = test_bacteria < percentile_train_97_5 & test_bacteria > percentile_train_2_5
    within50 = test_bacteria < percentile_train_75& test_bacteria > percentile_train_25
    
    
    
    
    
    #dummy_df$in95 <- dummy_df$in95 + test_beta(true = sum(within95), false = sum(!within95), percentile = .95 )
    #dummy_df$below95 <- dummy_df$below95+test_beta(true = sum(below95), false = sum(!below95), percentile = .95 )
    #dummy_df$below90 <- dummy_df$below90+test_beta(true = sum(below90), false = sum(!below90), percentile = .90 )
    #dummy_df$in50 <- dummy_df$in50+test_beta(true = sum(within50), false = sum(!within50), percentile = .5)
    
    
    return(test_beta(is_true = within95, percentile = .95 ) &
             test_beta(is_true = below95, percentile = .95 )&
             test_beta(is_true = below90, percentile = .90 )&
             test_beta(is_true = within50, percentile = .5)
    )
  }
}
{     
  validation_on_test_stepwise <- function(model, data_train, data_test){      
     #model <- step_5_aic_models_full_dataset 
     # data_train<- train_data_full_scaled
    #  data_test<- test_data_full_scaled
    
    #train_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, data_train)[,-1]
    #test_data_matrix_with_interactions <- sparse.model.matrix(log_e.coli~(.)^2, data_test)[,-1]
    
    train_bacteria  <- data_train$log_e.coli
    test_bacteria <- data_test$log_e.coli
    
    #prediction
    #pred_train <- predict(linear_model,newdata = data_train_fold)
    #pred_test <- predict(linear_model, newdata =  data_test_fold)
    pred_train <- model %>% predict(data_train) %>% as.vector()
    pred_test <- model %>% predict(data_test) %>% as.vector()
    
    
    #r2<- 1 - sum((data_train_fold$log_e.coli-pred_train)^2)/sum((data_train_fold$log_e.coli-mean(data_train_fold$log_e.coli))^2)
    
    #rsq <- function (x, y) cor(x, y) ^ 2
    #rsq( train_bacteria, pred_train)
    #e<-cbind(as.data.frame(train_bacteria), as.data.frame(pred_train))
    #e$error <- sqrt((e$train_bacteria - e$pred_train)^2)
    
    
    
    #adj_R2 <- 1-((1-r2)*(n_obs-1))/(n_obs-n_features-1)
    #list_R2_iteration <-adj_R2
    
    #get sd from difference between train_bacteria_true to the predicted! mean is the prediction
    sd_pred_train_to_true<-sqrt((sum((train_bacteria - pred_train)^ 2))/length(train_bacteria))
    
    #dummy_df$sd_pred_train_to_true <- sd_pred_train_to_true
    
    percentile_train_2_5<- qnorm(0.025, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_25<- qnorm(0.25, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_75<- qnorm(0.75, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_90<- qnorm(0.90, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_95<- qnorm(0.95, mean = pred_test, sd = sd_pred_train_to_true)
    percentile_train_97_5<- qnorm(0.975, mean = pred_test, sd = sd_pred_train_to_true)
    
    #evaluating ther model has to be classified correctly with every single test train split
    #--> here 5 different splits, if all validations correct than everywhere ==5
    
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

check_if_all_folds_validated <-function(validation_df){
  
  if(all(validation_df$in95==5,validation_df$below95==5,validation_df$below90==5,validation_df$in50==5)){
    return(T)
  }else if(!all(validation_df$in95==5,validation_df$below95==5,validation_df$below90==5,validation_df$in50==5)){
    return(F)
  }
}
get_validation_per_unique_formula_over_all_validion_iterations<-function(unique_river_algo_formulas, all_train_iterations_validation_scaled_river, all_test_iterations_validation_scaled_river){
  #unique_river_algo_formulas<-formula(step_5_aic_models_full_dataset)
  #unique_river_algo_formulas<-unique_formulas_river
  #all_train_iterations_validation_scaled_river<-validation_train_data_random_split
  #all_test_iterations_validation_scaled_river<-validation_test_data_random_split
  #<-all_test_iterations_validation_scaled_river[[river_number]]
  #all_test_iterations_validation_scaled_river<-all_test_iterations_validation_scaled_only_validation_data[[1]]
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

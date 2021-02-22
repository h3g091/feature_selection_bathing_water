#### validate on fifty validation iterations


library(glmnet)
#library(fhpredict)
library(tidyverse)
#important
library(kwb.flusshygiene)
library(dplyr)
library(Matrix)
library(caret)

### in source function
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
#automated method
#add path for data here
data_path_isar<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Bayern/Isar/DATA_preprocessed_csv")
data_path_ilz<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Bayern/Ilz/DATA_preprocessed_csv")
data_path_mosel<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Rhein_Mosel_Lahn/Mosel/DATA_preprocessed_csv")
data_path_rhein<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Rhein_Mosel_Lahn/Rhein/DATA_preprocessed_csv")
data_path_ruhr <- list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/Daten_Ruhr/Ruhr/DATA_preprocessed_csv")
data_path_kleine_badewiese<-list(havel = "/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/kleine_Badewiese/Havel/DATA_preprocessed_csv")


#river_name <- "kleine badewiese"
#data_path_river <- data_path_kleine_badewiese
#kleine_BW<- T

big_tibble<-readRDS("results/big_tibble/big_tibble.rds")

get_dataset_for_5_fold<-function(data_path, river_name, kleine_BW){
    #data_path<- data_path_isar
    #river_name <- "Isar"
    #kleine_BW <- F
    #build outer train test split --> use train and build 50 foldid iterations
    source("read_and_lag_preprocessed_data.R")
    river_data<-preprocess_river_data_from_path(data_path, river_name, kleine_BW)
    #rm(preprocess_river_data_from_path)
    
    seed <- 123
    
    
    #split data for fs
    
    set.seed(seed)
    training.samples <- river_data$log_e.coli %>%
      createDataPartition(p = 0.8, list = FALSE)
    train_data  <- river_data[training.samples, ]
    #test_data <- river_data[-training.samples, ]
    return(train_data)
}
build_multiple_validation_tibble <- function(data_path_river, river_name, kleine_BW, number_of_iterations=100) {
  list_method_features <- big_tibble%>% filter(River == river_name)%>% select(Method,features )
  data_set_for_five_fold<- get_dataset_for_5_fold(data_path_river, river_name, kleine_BW )
  
  #build 50 - 5-fold-cross with foldid
  n = nrow(data_set_for_five_fold)
  list_fold_ids <- list()
  iterations<- number_of_iterations
  fixed_seeds <- seq(1:iterations)
  for(seed in fixed_seeds){
    set.seed(seed)
    list_fold_ids<- append(list_fold_ids,list(sample(rep(seq(5),length=n))))
    
  }
  
  
  #split data for validation
  source("get_split_for_validation.R")
  
  list_train_test_split_data<-list()
  for(foldids in list_fold_ids){
    list_train_test_split_data<- append(list_train_test_split_data, list(get_split_for_validation(data_set_for_five_fold, foldids) ))
  }
  
  
  source("validation_methods.R")
  ##calculate validations for every lm and count how often it could be validated
  list_linear_models <- list()
  #rm(list_linear_models)
  
  #validation_df_step_aic <- validation_5_fold(features_step_aic_full, data_split_for_validation[[1]], data_split_for_validation[[2]])
  list_passed_folds_per_method <- list()
  for(method_features in list_method_features$features){
    
  
  list_number_of_passed_folds<- list()
    for(iteration_indx in 1:length(list_train_test_split_data)){
      list_number_of_passed_folds<- append(list_number_of_passed_folds, validation_5_fold(method_features, list_train_test_split_data[[iteration_indx]][[1]], list_train_test_split_data[[iteration_indx]][[2]])$passed_folds)
      
    }
  list_passed_folds_per_method <- append(list_passed_folds_per_method, list(unlist(list_number_of_passed_folds)))
  }  
  
  names(list_passed_folds_per_method)<-  list_method_features$Method
  river_validation_tibble <- tibble()
  for(method_indx in 1:6){
  #method_indx <- 3
    method_name<- list_method_features$Method[[method_indx]]
    five_fold_validation_river<- as.data.frame(list_passed_folds_per_method)
    
    #frequency_method<-count(five_fold_validation_river, step_5_aic)
    frequency_method<-count(five_fold_validation_river, !!sym(list_method_features$Method[method_indx]))
    names(frequency_method)<- c("N_Folds", "Frequency" )
    frequency_method$Method <- method_name
    river_validation_tibble<- rbind(river_validation_tibble, frequency_method)
  }
  
  river_validation_tibble$River <- river_name
  return(river_validation_tibble)
}

vallidation_occurence_kleine_badewiese<-build_multiple_validation_tibble(data_path_kleine_badewiese,"kleine badewiese", kleine_BW = T)
vallidation_occurence_isar<-build_multiple_validation_tibble(data_path_isar,"Isar", kleine_BW = F)
vallidation_occurence_ilz<-build_multiple_validation_tibble(data_path_ilz,"Ilz", kleine_BW = F)
vallidation_occurence_mosel<-build_multiple_validation_tibble(data_path_mosel,"Mosel", kleine_BW = F)
vallidation_occurence_rhein<-build_multiple_validation_tibble(data_path_rhein,"Rhein", kleine_BW = F)
vallidation_occurence_Ruhr<-build_multiple_validation_tibble(data_path_ruhr,"Ruhr", kleine_BW = F)

validation_occurence_big_table<-do.call("rbind", list(vallidation_occurence_kleine_badewiese,
                      vallidation_occurence_isar,
                      vallidation_occurence_ilz,
                      vallidation_occurence_mosel,
                      vallidation_occurence_rhein,
                      vallidation_occurence_Ruhr))
#five_fold_validation_river_tible<-as_tibble(rep(river_name, 6))
#five_fold_validation_river_tible<-cbind(five_fold_validation_river_tible, unique(big_tibble$Method))


#features_step_5_aic<-as.character(big_tibble%>% filter(Method== "step_5_aic", River == river_name)%>% select(features ))


river_validation_tibble<-validation_occurence_big_table %>% pivot_wider(names_from = N_Folds, values_from = Frequency)

###replace NA wit 0
long_temp<-river_validation_tibble%>% replace_na(list(`0`= 0,`1`= 0, `2`= 0,`3`= 0,`4`= 0,`5`= 0))

long_temp<-long_temp%>% gather(`Number of folds all criteria are passed`, `Occurence in 100 Iterations`, - Method, - River)
long_temp$River[long_temp$River== "kleine badewiese"]<- "Havel"

long_temp$River <- factor(long_temp$River, levels = c("Havel","Ilz","Isar", "Mosel", "Rhein", "Ruhr" ))
long_temp$Method <- factor(long_temp$Method, levels = c("step_5_aic","step_aic","step_5_bic", "step_bic", "lasso_min", "lasso_1se" ))

#long_temp$`Number of folds all criteria are passed`<-as.integer(long_temp$`Number of folds all criteria are passed`)
long_temp$validated <- F
long_temp$validated[long_temp$`Number of folds all criteria are passed` ==5] <- T
long_temp$validated<- as.factor(long_temp$validated)
#long_temp$Occurence<-long_temp$Occurence/iterations

p<-ggplot(data = long_temp, aes(x= `Number of folds all criteria are passed`, fill = validated))+
  geom_bar( aes(y = `Occurence in 100 Iterations`), stat = "identity", show.legend = F)+
  geom_text( aes(y =`Occurence in 100 Iterations`,label=`Occurence in 100 Iterations`),vjust=-0.4,size =3.5)+
  ylim(c(0,105))+
  theme(text = element_text(size = 16))       
  
  

p + facet_grid( River  ~ Method )
#p + facet_grid( Method ~ River  )

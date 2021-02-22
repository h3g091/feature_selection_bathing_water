get_split_for_validation <-function(data, foldid){
  #set.seed(random_seed)
  #number_of_folds <- number_of_folds_for_validation  
  number_of_folds<-max(foldid)
  #validation_samples<-sample(1:1000, size = 1)
  datapoints<-seq(1, nrow(data))
  #multiple train_test_splits for dsame dataframe
  
  #foldid=sample(rep(seq(number_of_folds),length=nrow(river_data)))
  train_rows<- list()
  test_rows<- list()
  for (fold in 1:number_of_folds) {
    train_rows_temp<-datapoints[foldid!=fold] 
    test_rows_temp<-datapoints[foldid==fold] 
    train_rows <- append(train_rows,list(data[train_rows_temp,]))
    test_rows <- append(test_rows,list(data[test_rows_temp,]))
  }
  
  return(list(train_rows, test_rows))
}

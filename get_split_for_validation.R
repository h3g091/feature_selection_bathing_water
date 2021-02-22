#' Splits up data with given foldid
#' @param data A table  of data with n rows
#' @param foldid A list of length n with values between 1 and m, m representing mthe maximum number of folds
#' @param return A list, that contains two lists, both with length m. Each train_rows element contains the datapoints of fold x for training and each corresponding test_rows all datapoints of fold x for testing
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

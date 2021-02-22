
###75 as 15*5 iterations!!
get_all_lm_features_from_list_as_string <- function(model_list_with_coef){
  return(as.character(formula(model_list_with_coef)))
  #return(unlist(lapply(1:75, function(x) lapply(lapply(model_list_with_coef, formula), as.character)[[x]][3])))
}



get_10_biggest_iterations_from_list_of_top_unique_validations<-function(unique_validations){
  unique_validations<-unique(sort (unlist(unique_validations), decreasing = T ))
  if(length(unique_validations[unique_validations!=0]) >10){
    max_10_iteration_numbers<-unique_validations[unique_validations!=0]
    
    return(max_10_iteration_numbers[1:10])
  }else if(length(unique_validations[unique_validations!=0]) <=10){
    max_10_iteration_numbers<-unique_validations[unique_validations!=0]
    return(max_10_iteration_numbers)
  }
}

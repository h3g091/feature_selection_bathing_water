#FS_algorithms

get_list_of_step_aic_models <- function(river_data_temp, state_of_the_art= F){
  #river_data_temp<- all_data_scaled_df[[2]]
  
  null <- lm(log_e.coli ~ 1, data = river_data_temp)
  full <- lm(log_e.coli ~ .^2, data = river_data_temp)
  if(state_of_the_art==F){
    selection_aic <- step(null, data = river_data_temp ,
                          
                          direction = "forward",
                          
                          list(lower=null, upper=full), k = 2)
    
  }else if(state_of_the_art!=F){
    selection_aic <- step(null, data = river_data_temp ,
                          
                          direction = "forward",
                          
                          list(lower=null, upper=full), k = 2, steps = 5)
    
  }
  return(selection_aic)
  
}
get_list_of_step_bic_models <- function(river_data_temp,  state_of_the_art= F){
  #river_data_temp<- all_data[[1]]
  n<-nrow(river_data_temp)
  null <- lm(log_e.coli ~ 1, data = river_data_temp)
  full <- lm(log_e.coli ~ .^2, data = river_data_temp)
  if(state_of_the_art==F){
    selection_bic <- step(null, data = river_data_temp ,
                          
                          direction = "forward",
                          
                          list(lower=null, upper=full), k = log(n))
    
  }else if(state_of_the_art!=F){
    selection_bic <- step(null, data = river_data_temp ,
                          
                          direction = "forward",
                          
                          list(lower=null, upper=full), k = log(n), steps = 5)
    
  }
  return(selection_bic)
  
}

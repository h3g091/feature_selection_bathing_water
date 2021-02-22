preprocess_river_data_from_path <- function(data_path, river_name, kleine_BW=F){
  #data_path<-"/Users/heiko.langer/Masterarbeit_lokal/Data_preprocess/kleine_Badewiese/Havel/DATA_preprocessed_csv"
  #river_name<- "kleine badewiese"
  river_paths<-data_path
  river_list_data<-list()
  
  
  
  {
    calc_t <- function (datalist=river_data$havel, onlysummer) {
      #heiko
      #datalist<- river_data1$havel
      phy_data <- datalist[-1] # Entfernung der Hygienedaten
      
      if(onlysummer==T){
        hyg_df <- subset(datalist[[1]],
                         subset = lubridate::month(datum) %in% 5:9) # 
        
        data_summer <- lapply(phy_data, function(df){
          
          df <- subset(df, subset = lubridate::month(datum) %in% 4:9) 
        }
        )
      }  
      
      
      # z_standardize <- function (x) {
      
      #   y = (x - mean(x, na.rm=T))/sd(x, na.rm=T)
      
      # }
      
      log_transorm_rain <- function(df) { #log transforming rain data
        
        for (site in names(df)[-1]) { # every col gets treatment
          
          df2 <- subset(df, select = c("datum", site))
          
          if (grepl("^r_.*",site)) { # rain gets log-transformed and 1/sigma2
            
            df2[[site]] <- log(df2[[site]]+1)
            
            # df2[[site]] <- df2[[site]]/sd(df2[[site]], na.rm=T)
            
          } #else {
          
          #   df[[site]] <- z_standardize(df2[[site]]) # standardize
          
          # }
          
          df[[site]] <- df2[[site]]
          
        }
        
        return(df)
        
      }
      
      data_t <- lapply(data_summer, log_transorm_rain)
      
      result <- append(list(hyg_df), data_t)
      
      names(result) <- names(datalist)
      
      return(result)
      
    }
  }
  
  river_data <- lapply(river_paths, kwb.flusshygiene::import_riverdata)
  names(river_data) <- river_name
  #renaming of columns
  if(kleine_BW ==T){
  colnames(river_data$`kleine badewiese`$r_radolan_2350)[2] <- "r"
    
  }
  
  
  #river with lag days
  {
    river_data_ts <- lapply(river_data, function(river_list){
      
      river_ts <- calc_t(river_list, onlysummer = T) # use function
      
      add_meancol <- function (df) { # for rain and i #edit: + ka #2ndedit: + q
        
        prefix <- unique(sub("([a-z])_.*","\\1",names(df)[-1]))
        
        for (pre in prefix) {
          
          df2 <- dplyr::select(df, dplyr::starts_with(pre))
          
          df[,paste0(pre,"_mean")] <- rowMeans(df2, na.rm=T)
          
        }
        
        
        
        return(df)
        
      }
      
      add_sumcol <- function (df) { # originally for ka, but not used
        
        prefix <- unique(sub("([a-z])_.*","\\1",names(df)[-1]))
        
        if (length(df) > 2)
          
          df[,paste0(prefix,"_sum")] <- rowSums(df[,-1], na.rm=T)
        
        return(df)
        
      }
      
      
      
      q_pos <- grep("^q", names(river_ts)[-1])+1
      
      
      if (length(q_pos) == 1)
        
        river_ts[[q_pos]] <- add_meancol(river_ts[[q_pos]])
      
      ka_pos <- grep("^ka", names(river_ts)[-1])+1
      
      if (length(ka_pos) == 1)
        
        river_ts[[ka_pos]] <- add_meancol(river_ts[[ka_pos]])
      
      i_pos <- grep("^i", names(river_ts)[-1])+1
      
      if (length(i_pos) == 1)
        
        river_ts[[i_pos]] <- add_meancol(river_ts[[i_pos]])
      
      r_pos <- grep("^r", names(river_ts)[-1])+1
      
      river_ts[[r_pos]] <- add_meancol(river_ts[[r_pos]])
      
      return(river_ts)
      
    })  
  }
  
  
  pattern = "(i_mean|q_mean|r_mean|ka_mean)"
  
  riverdata <- river_data_ts
  riverdata <- river_data_ts[[river_name]]
  
  # prepare variables out of all cominations (given by pattern)
  # variables for interaction get replaced by q_new (remove q_old)
  vars1 <- (riverdata[-1] %>% unroll_physical_data() %>%
              lapply(names) %>% unlist() %>% unique())[-1]
  
  vars2 <- vars1[stringr::str_detect(vars1, pattern)]
  # prepare formulas
  
  data <- process_model_riverdata(riverdata, c("log_e.coli", vars2)) %>%  dplyr::select(-datum) 
  #get rid of tiefwerder variable that is double
  data<-(data %>% select(- contains("tiefwerder")))
  data <- na.omit(data)
  
  return(data)   
}





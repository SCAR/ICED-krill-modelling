getTiles <- function(df, xTrans=-0.5+c(0,1,1), yTrans=-0.5+c(0,0,1)){
  require(dplyr)
  
  # function to get the points of a triangle
  getPoints <- function(vector, trans){
    for(i in 1:length(vector)){
      temp <- vector[i]+trans
      if(i==1){
        out <- temp
      }else{
        out <- c(out, temp)
      }
    }
    return(out)
  }
  
  df <- df %>%
    arrange(Latitude, Depth)%>%
    mutate(id = paste0("surv.", row_number()))
  x <- as.numeric(df$Depth)
  y <- as.numeric(df$Latitude)
  df <- rbind(df,df,df) 
  df <- df %>% arrange(Latitude, Depth)
  df$x <- getPoints(x, xTrans)
  df$y <- getPoints(y, yTrans)
  
  return(df)
}
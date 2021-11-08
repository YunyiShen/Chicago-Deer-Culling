extract_things_we_want <- function(scheme_list, what_we_want){
  list_stuff_we_want <- lapply(scheme_list,`[[`, what_we_want[[1]])
  #print(class(list_stuff_we_want[[1]]))
  dimensions <- dim(list_stuff_we_want[[1]])
  #print(dimensions)
  res <- matrix(NA,dimensions[1], dimensions[2])
  rownames(res) <- c(paste0("F",1:8),paste0("M",1:3))
  colnames(res) <- 1991+(1:dimensions[2])
  
  for(i in 1:dimensions[1]){
    for (j in 1:dimensions[2]) {
      res_vec <- sapply(list_stuff_we_want,`[`,i,j)
      #print(res_vec)
      if(all(is.na(res_vec))){
        res[i,j] <- NA
      }
      else{
        res[i,j] <- quantile(res_vec, probs = what_we_want[[2]], na.rm = TRUE)
      }
    }
  }
  return(res)
}
# --- functions --- #
#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param cv_group (numeric) cv_group
#' @param k_fold (numeric) k_fold cv
#' @return (numeric) errors

error_compute=function(data_X,data_y,cv_group,k_fold){
  sum_error=0
  if(ncol(data_X%>%as.matrix)>1){
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i,]
      data_X_test=data_X[!cv_group==i,]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      Xi= matlib::Ginv(data_X_train) %*% data_y_train
      sum_error=sum_error+sum(data_y_test-data_X_test%*%Xi)^2
      # mod=lm(data_y_train~data_X_train-1)
      # sum_error=sum_error+sum(data_y_test-data_X_test%*%(mod$coefficients))^2
    }
  }
  else{
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i]
      data_X_test=data_X[!cv_group==i]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      Xi= cov(data_X_train,data_y_train)/cov(data_X_train,data_X_train)
      sum_error=sum_error+sum(data_y_test-data_X_test*Xi)^2
      # mod=lm(data_y_train~data_X_train-1)
      # sum_error=sum_error+sum(data_y_test-data_X_test*(mod$coefficients))^2
    }
  }
  return(sum_error)
}



#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param cv_group (numeric) cv_group
#' @param k_fold (numeric) k_fold cv
#' @return (numeric) errors

error_compute_MASS=function(data_X,data_y,cv_group,k_fold){
  sum_error=0
  if(ncol(data_X%>%as.matrix)>1){
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i,]
      data_X_test=data_X[!cv_group==i,]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      Xi= MASS::ginv(data_X_train,0) %*% data_y_train
      sum_error=sum_error+sum(data_y_test-data_X_test%*%Xi)^2
      # mod=lm(data_y_train~data_X_train-1)
      # sum_error=sum_error+sum(data_y_test-data_X_test%*%(mod$coefficients))^2
    }
  }
  else{
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i]
      data_X_test=data_X[!cv_group==i]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      Xi= cov(data_X_train,data_y_train)/cov(data_X_train,data_X_train)
      sum_error=sum_error+sum(data_y_test-data_X_test*Xi)^2
      # mod=lm(data_y_train~data_X_train-1)
      # sum_error=sum_error+sum(data_y_test-data_X_test*(mod$coefficients))^2
    }
  }
  return(sum_error)
}




#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param cv_group (numeric) cv_group
#' @param k_fold (numeric) k_fold cv
#' @return (numeric) errors

error_compute_lm=function(data_X,data_y,cv_group,k_fold){
  sum_error=0
  if(ncol(data_X%>%as.matrix)>1){
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i,]
      data_X_test=data_X[!cv_group==i,]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      mod=lm(data_y_train~data_X_train-1)
      sum_error=sum_error+sum(data_y_test-data_X_test%*%(mod$coefficients))^2
    }
  }
  else{
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i]
      data_X_test=data_X[!cv_group==i]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      mod=lm(data_y_train~data_X_train-1)
      sum_error=sum_error+sum(data_y_test-data_X_test*(mod$coefficients))^2
    }
  }
  return(sum_error)
}


#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @return (numeric) How many time a variable is chosen

order_coef=function(data_X,data_y){
  n=ncol(data_X)
  order_col=rep(n,n)
  coef_dot=rep(TRUE,n)
  Xi=rep(0,n)
  
  for(i in 1:(n-1)){
    Xi[coef_dot]= MASS::ginv(data_X[,coef_dot],0) %*% data_y
    order_col[coef_dot][which.min(abs(Xi[coef_dot]))]=i
    coef_dot[coef_dot][which.min(abs(Xi[coef_dot]))]=FALSE
  }
  return(order_col)
}





#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param k_fold (numeric) k_fold cv
#' @return (list) list of some results

sSINDy=function(data_X,data_y,k_fold){
  library(magrittr)

  n_row=nrow(data_X)
  n_col=ncol(data_X)
  cv_group=sample(rep(1:k_fold,ceiling(n_row/k_fold)),n_row)
  errors=rep(0,n_col)
  coef_order=order_coef(data_X,data_y)
  coef_final=rep(0,n_col)

  errors=sapply((1:n_col),function(i)
    {error_compute_lm(data_X[,coef_order>=(n_col-i+1)],data_y,cv_group,k_fold)})
  # for(i in 1:n_col)
  # {
  #   errors[i]=error_compute(data_X[,coef_order>=(n_col-i+1)],data_y,cv_group,k_fold)
  # }
  min_coef=which.min(errors)
  coef_final[coef_order>(n_col-min_coef)]=
    lm(data_y~data_X[,coef_order>(n_col-min_coef)]-1)%>%coef()
  
  return(list("coef_order"=coef_order,"cv_error"=errors,
              "min#"=min_coef,"coef"=coef_final))
}



#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param k_fold (numeric) k_fold cv
#' @return (list) list of some results

sSINDy_lm=function(data_X,data_y,k_fold){
  library(magrittr)
  
  n_row=nrow(data_X)
  n_col=ncol(data_X)
  cv_group=sample(rep(1:k_fold,ceiling(n_row/k_fold)),n_row)
  errors=rep(0,n_col)
  coef_order=order_coef(data_X,data_y)
  coef_final=rep(0,n_col)
  
  
  for(i in 1:n_col)
  {
    errors[i]=error_compute_lm(data_X[,coef_order>=(n_col-i+1)],data_y,cv_group,k_fold)
  }
  min_coef=which.min(errors)
  coef_final[coef_order>(n_col-min_coef)]=
    lm(data_y~data_X[,coef_order>(n_col-min_coef)]-1)%>%coef()
  
  return(list("coef_order"=coef_order,"cv_error"=errors,
              "min#"=min_coef,"coef"=coef_final))
}


#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param k_fold (numeric) k_fold cv
#' @return (list) list of some results

sSINDy_lm_vectorized=function(data_X,data_y,k_fold){
  library(magrittr)
  
  n_row=nrow(data_X)
  n_col=ncol(data_X)
  cv_group=sample(rep(1:k_fold,ceiling(n_row/k_fold)),n_row)
  errors=rep(0,n_col)
  coef_order=order_coef(data_X,data_y)
  coef_final=rep(0,n_col)
  
  
  errors=sapply((1:n_col),function(i){error_compute_MASS(data_X[,coef_order>=(n_col-i+1)],data_y,cv_group,k_fold)})
  
  min_coef=which.min(errors)
  coef_final[coef_order>(n_col-min_coef)]=
    lm(data_y~data_X[,coef_order>(n_col-min_coef)]-1)%>%coef()
  
  return(list("coef_order"=coef_order,"cv_error"=errors,
              "min#"=min_coef,"coef"=coef_final))
}



###############################################################################
###This function calculates Mean recurrence time of a markov chain

ergodic_projector <- function(P, n){
  n <- n
  A <- array(1, dim=c(dim(P)[1], dim(P)[2] ,n))
  for (i in 0:n){
    
    result <- P %^% i

    A[,,i] <- result
    
  }

  #print(X)
  output <- apply(A, c(1,2), mean, na.rm = TRUE)
  return(output)

}



deviation_matrix <- function(P,n){
  ep <- ergodic_projector(P,n)
  c <- ncol(P)
  D <- inv(diag(c) - P + ep) - ep
  return(D)
}


MRT <- function(P,n){
  c <- ncol(P)
  D <- (deviation_matrix(P,n))
  E <- diag(diag(ergodic_projector(P,n)), c,c)
  R <- ((diag(c)) - D + (as.matrix(rep(1,c))) %*% t(as.matrix(rep(1,c))* diag(D))) %*% inv(E)
  
  return(R)
}


MRT_distance <- function(P,n){
  M <- MRT(P,n)
  average <- matrix(0, ncol(M), ncol(M))
  for (i in 1:ncol(M)){
    for (j in 1:ncol(M)){
      
    average[i,j] <- mean(as.vector(c(M[i,j], M[j,i])))
      
      
    }
  
  }

  return(average) 
  
}



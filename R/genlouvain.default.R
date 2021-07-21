metanetwork <- function(J,S)
{
  PP <- sparseMatrix(i=1:length(S),j=S,x=1)
  t(PP) %*% J %*% PP 
}

tidyconfig <- function(S)
{
  T <- rep(0,length(S))
  for (i in 1:length(S))
  {
    if (T[i]==0) T[S==S[i]] <- max(T)+1
  }
  T
}

genlouvain.default<-function(B){
  eps <- 2e-16
  if (sum(B-t(B))>0) B <- (B+t(B))/2  #force symmetric matrix
  M <- B
  
  n <- nrow(B)
  S <- t(1:n)
  
  dtot <- 0
  
  S2 <- 1:n
  Sb <- NULL
  
  n.outer <- 0
  
  while(!identical(Sb,S2))
  {
    n.outer <- n.outer+1
    if (n.outer>50) {
      print("Reached greater than 50 outer iterations.")
      return(0)
    }
    
    y <- unique(S2)
    y <- y[order(y,decreasing=F)]
    Sb <- S2
    print(paste(c("Merging",length(y),"communities"),collapse=" "))
    
    yb <- NULL
    
    G <- sparseMatrix(i=1:length(y),j=y,x=1)
    dstep <- 1
    nsteps <- 0
    
    while((!identical(yb,y)) && (dstep/dtot>2*eps))
    {
      yb <- y
      dstep <- 0
      nsteps <- nsteps+1
      print(nsteps)
      if (nsteps>50) {
        print("Reached greater than 50 inner iterations.")
        return(0)
      }
      
      #ord.i <- sample(1:nrow(M),nrow(M),replace=F)
      ord.i <- 1:nrow(M)
      
      for (i in ord.i)  
      {
        u <- unique(c(y[i],y[M[,i]>0]))
        u <- u[order(u,decreasing=F)]
        dH <- t(M[,i]) %*% G[,u]
        
        yi <- which(u==y[i])
        dH[yi] <- dH[yi] - M[i,i]
        k <- max.col(dH)
        #if (length(k)>1) k <- sample(k,1)
        
        if (dH[k]>dH[yi]){
          dtot <- dtot+dH[k]-dH[yi]
          dstep <- dstep+dH[k]-dH[yi]
          G[i,y[i]] <- 0
          G[i,u[k]] <- 1
          y[i] <- u[k]
        }
      }
    }
    
    y <- tidyconfig(y)

    
    for (i in 1:length(y))
    {
      S[S==i] <- y[i]
      S2[S2==i] <- y[i]
    }
    
    if (identical(Sb,S2))
    {
      return(S)
    }
    
    M <- metanetwork(B,S2)		
  }
}


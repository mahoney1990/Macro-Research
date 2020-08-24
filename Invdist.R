invdist<-function(prob){
  tol = max((1e-6)/length(prob[,1]),1e-11)
  maxiter = max(800,200*length(prob[,1]))
  startDist = (1/length(prob[,1]))*matrix(1,length(prob[,1]),length(prob[,1]))
  
  nrows=c(as.numeric(length(prob[1,])))
  ncol=c(as.numeric(length(prob[,1])))
  
  if(nrows != ncol){
    print("You have fucked up")  
    quit()
  }
  
  if(colSums(t(prob)) != matrix(1,1,nrows)){
    print("You have fucked up")  
    quit()}
  
  print("No Errors")
  
  invdist = startDist;        
  test = 1.0                  
  niter = 1 
  diff=10
  while(diff>tol){
    if(niter>maxiter){
      print("Max Iterations Reached!")
      quit()
    }
    probdist=t(prob) %*% invdist
    diff=max(abs(probdist-invdist))
    invdist=probdist
    niter=niter+1
  }
  
  return(t(invdist))}
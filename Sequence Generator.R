

generate_seq <- function(n)
{
  x    <- numeric(n)
  x[1] <- sample(4, 1)
  had_a_four <- FALSE
  
  for(i in seq(n - 1)) {
    if(!had_a_four)
    {
      x[i + 1]  <- switch(x[i], sample(1:2, 1, prob = c(3, 1)), 
                          sample(2:3, 1, prob = c(3, 1)), 
                          sample(3:4, 1, prob = c(3, 1)), 
                          sample(4, 1))
    }
    else
    {
      x[i + 1]  <- switch(x[i], sample(1:2, 1, prob = c(3, 1)), 
                          sample(2:3, 1, prob = c(3, 1)), 
                          sample(3:4, 1, prob = c(3, 1)),  
                          4)
    }
    if(x[i + 1] == 4 & !all(x[1:(i+1)] == 4)) had_a_four <- TRUE
  }
  x
}


X<-matrix(replicate(2000000, generate_seq(12)), ncol = 2000000)
X<-X[, !duplicated(t(X))]
write.csv(X,"Valid_Seq_4dim.csv")



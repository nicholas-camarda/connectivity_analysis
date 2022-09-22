library(reprex)
my_function <- function(){
  res <- rbinom(n = 8, size = 1, prob = c(1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8))
  return(sum(res)/length(res) == 5/8)
}
reps <- replicate(n = 10000, expr = my_function(), simplify = TRUE)
sum(reps)/length(reps)

reprex()
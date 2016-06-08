#Set initial parameters
#Set an initial resource level
K <- 1000
# Lambda is growth rate defined as (N(t)/N(t-1))-1
LAMBDA<-1.9
#N is the initial population
N = 1000

#Get empty vector to hold populations
POPULATIONS<-NULL

#Append initial value
POPULATIONS<-c(POPULATIONS,N)

for (i in 1:999){
  N<-N + (LAMBDA*N*(1-(N/K)))
  POPULATIONS<-c(POPULATIONS,N)
  #Vary carrying capacity
  K<-rnorm(1,1000,100)
}

plot.ts(POPULATIONS,ylim=c(0,1500))
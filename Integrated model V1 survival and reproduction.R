#Sets up a timer
a<-proc.time()

#So here we want to build a model of differential reproductive success for male strategies.
#We'll use the bones of the simple drift model to figure out how these should fit together.

#Set number of generations
GENERATIONS=1000
ITERATIONS=100

#Set up data collection vectors over the iterations
DATA_EXTINCTION<-NULL
DATA_FIXED<-NULL

for (j in 1:ITERATIONS){
  
  #Set initial parameter values for the allele and sex frequencies: 
  POPULATION = 1000
  ALLELETERR = 0.5
  ALLELESNEAK = 1-ALLELETERR
  SEXRATIOF = 0.5
  SEXRATIOM = 1-SEXRATIOF
  MALE_TERR_SUCCESS = 0.5
  MALE_SNEAK_SUCCESS = 0.5
  
  #Ecological survival parameters based on K1 and K2 (carrying capacity for terr. males and females+sneakers)
  K1_mean = 500
  K1_STD = 10
  K2_mean = 625
  K2_STD = 15
  
  #Get number of individuals for each sex and genotype
  #Set initial values to 0 outside loop
  NUMBERS_FEMALE<-rmultinom(1,POPULATION*SEXRATIOF,c(ALLELETERR^2,ALLELETERR*2*ALLELESNEAK,ALLELESNEAK^2))
  NUMBERS_MALE<-rmultinom(1,POPULATION*SEXRATIOM,c(ALLELETERR^2,ALLELETERR*ALLELESNEAK*2,ALLELESNEAK^2))
  DATA_MALE<-NULL
  DATA_FEMALE<-NULL
  DATA_S<-NULL
  DATA_T<-NULL
  DATA_POPULATION<-NULL
  DATA_K1<-NULL
  DATA_K2<-NULL
  
  
  for (i in 1:GENERATIONS){
    
    #New generation
    #Probabilities for different genotypes of offpring
    
    #Here's the first place where our code changes. Now we want to first find out what proportion of male
    #even survive to breed in the first place, plus payouts for males that do survive. 
    #So how do we want to model this? 
    
    #Creating these things for the sake of convenience
    TOTAL_FEMALE<-sum(NUMBERS_FEMALE)
    TOTAL_TERR<-sum(NUMBERS_MALE[1:2])
    TOTAL_SNEAK<-NUMBERS_MALE[3]
    #Do survival. FEMSNEAK is ALL the females and sneakers--lumped together--not sneaker females
    K1<-max(0,rnorm(1,K1_mean,K1_STD))
    K2<-max(0,rnorm(1,K2_mean,K2_STD))
    
    if(K1<TOTAL_TERR&TOTAL_TERR>0){
      NUMBERS_MALE[1:2]<-rmultinom(1,K1,c(NUMBERS_MALE[1]/TOTAL_TERR,NUMBERS_MALE[2]/TOTAL_TERR))
    } 
    if(K2<TOTAL_SNEAK+TOTAL_FEMALE){
      NUMBERS_FEMSNEAK<-rmultinom(1,min(K2,TOTAL_SNEAK+TOTAL_FEMALE),c((TOTAL_SNEAK/(TOTAL_FEMALE+TOTAL_SNEAK)),(NUMBERS_FEMALE[1]/(TOTAL_FEMALE+TOTAL_SNEAK)),(NUMBERS_FEMALE[2]/(TOTAL_FEMALE+TOTAL_SNEAK)),(NUMBERS_FEMALE[3]/(TOTAL_FEMALE+TOTAL_SNEAK))))
      NUMBERS_MALE[3]<-NUMBERS_FEMSNEAK[1]
      NUMBERS_FEMALE<-NUMBERS_FEMSNEAK[2:4]
    } 
    
    
    #Now we look at females. Each female has a different likelihood of meeting with males of different strategies.
    #So here, we're working in the mating success with females that males of each strategy has.
    #First, get number of mating males of all three types
    NUMBER_MALE_SS_MATING<-ceiling(NUMBERS_MALE[3]*MALE_SNEAK_SUCCESS)
    NUMBER_MALE_TT_MATING<-ceiling(NUMBERS_MALE[1]*MALE_TERR_SUCCESS)
    NUMBER_MALE_ST_MATING<-ceiling(NUMBERS_MALE[2]*MALE_TERR_SUCCESS)
    TOTAL_MALE_MATING<-NUMBER_MALE_SS_MATING+NUMBER_MALE_ST_MATING+NUMBER_MALE_TT_MATING
    #Now get probability of each offsping genotype for each female genotype
    SS_TO_SS<-(NUMBER_MALE_SS_MATING/TOTAL_MALE_MATING)+(0.5*NUMBER_MALE_ST_MATING/TOTAL_MALE_MATING)
    SS_TO_ST<-(NUMBER_MALE_TT_MATING/TOTAL_MALE_MATING)+(0.5*NUMBER_MALE_ST_MATING/TOTAL_MALE_MATING)
    TT_TO_TT<-(NUMBER_MALE_TT_MATING/TOTAL_MALE_MATING)+(0.5*NUMBER_MALE_ST_MATING/TOTAL_MALE_MATING)
    TT_TO_ST<-(NUMBER_MALE_SS_MATING/TOTAL_MALE_MATING)+(0.5*NUMBER_MALE_ST_MATING/TOTAL_MALE_MATING)
    ST_TO_TT<-(0.5*NUMBER_MALE_TT_MATING/TOTAL_MALE_MATING)+(0.25*NUMBER_MALE_ST_MATING/TOTAL_MALE_MATING)
    ST_TO_SS<-(0.5*NUMBER_MALE_SS_MATING/TOTAL_MALE_MATING)+(0.25*NUMBER_MALE_ST_MATING/TOTAL_MALE_MATING)
    ST_TO_ST<-(1-ST_TO_TT-ST_TO_SS)
    
    #Now calculating offspring numbers from females of each genotype...
    OFFSPRING_OF_SS<-rmultinom(1,NUMBERS_FEMALE[3]*2,c(SS_TO_SS,SS_TO_ST))
    OFFSPRING_OF_TT<-rmultinom(1,NUMBERS_FEMALE[1]*2,c(TT_TO_ST,TT_TO_TT))
    OFFSPRING_OF_ST<-rmultinom(1,NUMBERS_FEMALE[2]*2,c(ST_TO_TT,ST_TO_ST,ST_TO_SS))
    
    #For next generation, let's find out how many offspring of each genotype...
    POPULATION<-sum(OFFSPRING_OF_TT,OFFSPRING_OF_ST,OFFSPRING_OF_SS)
    TOTAL_SS<-sum(OFFSPRING_OF_SS[1],OFFSPRING_OF_ST[3])
    TOTAL_ST<-sum(OFFSPRING_OF_SS[2],OFFSPRING_OF_ST[2],OFFSPRING_OF_TT[1])
    TOTAL_TT<-sum(OFFSPRING_OF_ST[1],OFFSPRING_OF_TT[2])
    FREQ_T<-(2*TOTAL_TT + TOTAL_ST) / (POPULATION*2)
    FREQ_S<-(2*TOTAL_SS + TOTAL_ST) / (POPULATION*2)
    
    #Calculating our numbers of males and females for the next gen...
    SEXRATIOF=rbeta(1,11,11)
    SEXRATIOM = 1-SEXRATIOF
    NUMBERS_FEMALE<-c(ceiling(SEXRATIOF*TOTAL_TT),ceiling(SEXRATIOF*TOTAL_ST),ceiling(SEXRATIOF*TOTAL_SS))
    NUMBERS_MALE<-c(floor(SEXRATIOM*TOTAL_TT),floor(SEXRATIOM*TOTAL_ST),floor(SEXRATIOM*TOTAL_SS))
    DATA_MALE<-rbind(DATA_MALE,NUMBERS_MALE)
    DATA_FEMALE<-rbind(DATA_FEMALE,NUMBERS_FEMALE)
    DATA_S<-c(DATA_S,FREQ_S)
    DATA_T<-c(DATA_T,FREQ_T)
    DATA_POPULATION<-c(DATA_POPULATION,POPULATION)
    DATA_K1<-c(DATA_K1,K1)
    DATA_K2<-c(DATA_K2,K2)
  }

X<-seq(GENERATIONS)
#par(mfrow=c(3,1))
plot(cbind(X,DATA_T),ylim=c(0,1),col='red')
points(cbind(X,DATA_S))
#plot.ts(DATA_T)
#plot.ts(DATA_MALE[,3])

  #Did the population go extinct?
  if(DATA_POPULATION[GENERATIONS]== 0){
    EXTINCTION = 'YES'
  } else {EXTINCTION='NO'}
  DATA_EXTINCTION<-c(DATA_EXTINCTION,EXTINCTION)
  
  #Did something fix?
  if (DATA_S[GENERATIONS]==0){
    FIXED='T'
  } else if (DATA_S[GENERATIONS]==1) {
    FIXED = 'S'
  } else {FIXED='N'}
  DATA_FIXED<-c(DATA_FIXED,FIXED)
}

sum(DATA_EXTINCTION=='YES')
sum(DATA_FIXED=='NO')
#Returns runtime at the end
proc.time()-a
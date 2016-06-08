#So here we want to build a model of differential reproductive success for male strategies.
#We'll use the bones of the simple drift model to figure out how these should fit together.

#Set initial parameter values for the allele and sex frequencies: 
POPULATION = 1000
ALLELETERR = 0.2
ALLELESNEAK = 0.6
SEXRATIOF = 0.5
SEXRATIOM = 1-SEXRATIOF
MALE_TERR_SUCCESS = 0.8
MALE_SNEAK_SUCCESS = 0.4

#Ecological survival parameters based on K1 and K2
K1_mean = 375
K1_STD = 100
K2_mean = 625
K2_STD = 150

#Get number of individuals for each sex and genotype
#Use a for loop through number of individuals plus if statement
#There has to be a better way to optimize this
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


for (i in 1:1000){
  
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
  K1<-rnorm(1,K1_mean,K1_STD)
  K2<-rnorm(1,K2_mean,K2_STD)
  NUMBERS_MALE[1:2]<-rmultinom(1,min(K1,TOTAL_TERR),c(NUMBERS_MALE[1]/TOTAL_TERR,NUMBERS_MALE[2]/TOTAL_TERR))
  NUMBERS_FEMSNEAK<-rmultinom(1,min(K2,TOTAL_SNEAK+TOTAL_FEMALE),c((TOTAL_SNEAK/(TOTAL_FEMALE+TOTAL_SNEAK)),(NUMBERS_FEMALE[1]/(TOTAL_FEMALE+TOTAL_SNEAK)),(NUMBERS_FEMALE[2]/(TOTAL_FEMALE+TOTAL_SNEAK)),(NUMBERS_FEMALE[3]/(TOTAL_FEMALE+TOTAL_SNEAK))))
  NUMBERS_MALE[3]<-NUMBERS_FEMSNEAK[1]
  NUMBERS_FEMALE<-NUMBERS_FEMSNEAK[2:4]
  
  #Now we look at females. Each female has a different likelihood of meeting with males of different strategies.
  #So here, we're working in the mating success with females that males of each strategy has.
  SS_TO_SS<-(MALE_SNEAK_SUCCESS*NUMBERS_MALE[3]/sum(NUMBERS_MALE))+(MALE_TERR_SUCCESS*0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  SS_TO_ST<-(MALE_TERR_SUCCESS*NUMBERS_MALE[1]/sum(NUMBERS_MALE))+(MALE_TERR_SUCCESS*0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  TT_TO_TT<-(MALE_TERR_SUCCESS*NUMBERS_MALE[1]/sum(NUMBERS_MALE))+(MALE_TERR_SUCCESS*0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  TT_TO_ST<-(MALE_SNEAK_SUCCESS*NUMBERS_MALE[3]/sum(NUMBERS_MALE))+(MALE_TERR_SUCCESS*0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  ST_TO_TT<-((MALE_TERR_SUCCESS*0.5*(NUMBERS_MALE[1]/sum(NUMBERS_MALE)))+(0.25*(MALE_TERR_SUCCESS*NUMBERS_MALE[2]/sum(NUMBERS_MALE))))
  ST_TO_SS<-((0.5*(MALE_SNEAK_SUCCESS*NUMBERS_MALE[3]/sum(NUMBERS_MALE)))+(0.25*(MALE_TERR_SUCCESS*NUMBERS_MALE[2]/sum(NUMBERS_MALE))))
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

X<-seq(1000)
#par(mfrow=c(3,1))
plot(cbind(X,DATA_T),ylim=c(0,1),col='red')
points(cbind(X,DATA_S))
plot.ts(DATA_POPULATION)
plot.ts(DATA_MALE[,3])

#Time to spit out the record-keeping numbers for each run:
if (DATA_POPULATION[1000]==0){
  EXTINCTION="YES"
} else {EXTINCTION="NO"}
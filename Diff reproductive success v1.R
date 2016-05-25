#So here we want to build a model of differential reproductive success for male strategies.
#We'll use the bones of the simple drift model to figure out how these should fit together.

#Set initial parameter values for the allele and sex frequencies: 
POPULATION = 1000
ALLELETERR = 0.5
ALLELESNEAK = 0.5
SEXRATIOF = 0.5
SEXRATIOM = 1-SEXRATIOF
MALE_TERR_SURVIVAL = 0.3
MALE_SNEAK_SURVIVAL= 0.8
FEMALE_SURVIVAL = 1
MALE_TERR_SUCCESS = 0.7
MALE_SNEAK_SUCCESS = 0.1

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

for (i in 1:1000){
  
  #New generation
  #Probabilities for different genotypes of offpring
  
  #Here's the first place where our code changes. Now we want to first find out what proportion of male
  #even survive to breed in the first place, plus payouts for males that do survive. 
  #So how do we want to model this? 
  
  #First we need to account for male phenotype as being separate from genotype.... or do we? 
  #No, we can just lump NUMBERS_MALE numbers. Okay. So now we want to incorporate male survival before reproduction.
  NUMBERS_MALE[1]<-NUMBERS_MALE[1]*MALE_TERR_SURVIVAL
  NUMBERS_MALE[2]<-NUMBERS_MALE[2]*MALE_TERR_SURVIVAL
  NUMBERS_MALE[3]<-NUMBERS_MALE[3]*MALE_SNEAK_SURVIVAL
  #Now we look at females. Each female has a different likelihood of meeting with males of different strategies.
  #So here, we're working in the mating success with females that males of each strategy has.
  #Okay. So let's say that here, 80% of 
 
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
  NUMBERS_FEMALE<-c(ceiling(0.5*TOTAL_TT),ceiling(0.5*TOTAL_ST),ceiling(0.5*TOTAL_SS))
  NUMBERS_MALE<-c(floor(0.5*TOTAL_TT),floor(0.5*TOTAL_ST),floor(0.5*TOTAL_SS))
  DATA_MALE<-rbind(DATA_MALE,NUMBERS_MALE)
  DATA_FEMALE<-rbind(DATA_FEMALE,NUMBERS_FEMALE)
  DATA_S<-c(DATA_S,FREQ_S)
  DATA_T<-c(DATA_T,FREQ_T)
}

X<-seq(1000)
plot(cbind(X,DATA_T),ylim=c(0,1),col='red')
points(cbind(X,DATA_S))
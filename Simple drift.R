#Set initial parameter values - population size and allele frequency
POPULATION = 1000
ALLELETERR = 0.5
ALLELESNEAK = 0.5
SEXRATIOF = 0.5
SEXRATIOM = 1-SEXRATIOF

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
  #HS - homozygote sneaker, HT - homozygote territorial, HET - heterozygote
  #Notation: e.g. HSTOHS - probability of HS female having HS offspring
  SS_TO_SS<-(NUMBERS_MALE[3]/sum(NUMBERS_MALE))+(0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  SS_TO_ST<-(NUMBERS_MALE[1]/sum(NUMBERS_MALE))+(0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  TT_TO_TT<-(NUMBERS_MALE[1]/sum(NUMBERS_MALE))+(0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  TT_TO_ST<-(NUMBERS_MALE[3]/sum(NUMBERS_MALE))+(0.5*NUMBERS_MALE[2]/sum(NUMBERS_MALE))
  ST_TO_TT<-((0.5*(NUMBERS_MALE[1]/sum(NUMBERS_MALE)))+(0.25*(NUMBERS_MALE[2]/sum(NUMBERS_MALE))))
  ST_TO_SS<-((0.5*(NUMBERS_MALE[3]/sum(NUMBERS_MALE)))+(0.25*(NUMBERS_MALE[2]/sum(NUMBERS_MALE))))
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
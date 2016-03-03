#Set initial parameter values - population size and allele frequency
POPULATION = 1000
ALLELETERR = 0.9
ALLELESNEAK = 0.1
SEXRATIOF = 0.5
SEXRATIOM = 1-SEXRATIOF

#Get number of individuals for each sex and genotype
#Use a for loop through number of individuals plus if statement
#There has to be a better way to optimize this
#Set initial values to 0 outside loop
NUMBERSFEMALE<-rmultinom(1,POPULATION*SEXRATIOF,c(ALLELETERR^2,ALLELETERR*2*ALLELESNEAK,ALLELESNEAK^2))
NUMBERSMALE<-rmultinom(1,POPULATION*SEXRATIOM,c(ALLELETERR^2,ALLELETERR*ALLELESNEAK*2,ALLELESNEAK^2))

#New generation
#Probabilities for different genotypes of offpring
#HS - homozygote sneaker, HT - homozygote territorial, HET - heterozygote
#Notation: e.g. HSTOHS - probability of HS female having HS offspring
HSTOHS<-(NUMBERSMALE[3]/sum(NUMBERSMALE))+(0.5*NUMBERSMALE[2]/sum(NUMBERSMALE))
HSTOHET<-(NUMBERSMALE[1]/sum(NUMBERSMALE))+(0.5*NUMBERSMALE[2]/sum(NUMBERSMALE))

OFFSPRINGHOMOSNEAK<-rmultinom(1,NUMBERSFEMALE[3]*2,c(HSTOHS,HSTOHET))
OFFSPRINGHET<-rmultinom(1,NUMBERSFEMALE[2]*2,c((NUMBERSMALE[3]/sum(NUMBERSMALE))+(0.5*NUMBERSMALE[2]/sum(NUMBERSMALE)),(NUMBERSMALE[1]/sum(NUMBERSMALE))+(0.5*NUMBERSMALE[2]/sum(NUMBERSMALE))))

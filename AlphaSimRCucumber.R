# install.packages(pkgs = "AlphaSimR")


library(package = "AlphaSimR")


#here I am adding detail that is not fully used in the simplified version
Nchromosomes<-7
ChrLengths<-rep(100,Nchromosomes)
SnpsDensityInCM<-1
Nsnps<-sum(ChrLengths)/SnpsDensityInCM
NfounderIndividuals<-100

  
  #use maCS
  ChrLength<-mean(ChrLengths)
  NsnpsPerChr<-Nsnps/Nchromosomes
  
  FounderPop<-runMacs(
    nInd=NfounderIndividuals,
    nChr = Nchromosomes,
    segSites = NsnpsPerChr,
    inbred = T,
    species = "GENERIC",
    split = NULL,
    ploidy = 2L,
    manualCommand = NULL,
    manualGenLen = NULL,
    nThreads = NULL
  )

#alternative function:
# runMacs2(
#   nInd,
#   nChr = 1,
#   segSites = NULL,
#   Ne = 100,
#   bp = 1e+08,
#   genLen = 1,
#   mutRate = 2.5e-08,
#   histNe = c(500, 1500, 6000, 12000, 1e+05),
#   histGen = c(100, 1000, 10000, 1e+05, 1e+06),
#   inbred = FALSE,
#   split = NULL,
#   ploidy = 2L,
#   returnCommand = FALSE,
#   nThreads = NULL
# )

FounderPop

#set simulation parameters
SP <- SimParam$new(FounderPop)
SP

#hybrid plot size 10
#10QTLs
Nqtl<-2#Per chromosome!!!

SP$nTraits
SP$traits
SP$varA
SP$varE


# Add trait (can include additive, dominance, epistasis and GxE) (can give vectors for multiple traits)
SP$addTraitA(nQtlPerChr=Nqtl ,mean = 0, var = 1, corA = NULL, gamma = FALSE, shape = 1, force = FALSE)
SP$varA
SP$varE
SP$nTraits
SP$traits
SP$traits[[1]]

getQtlMap()


SP$tsepMapSP$traitNames
SP$traitNames<-"Digibitism"


# Set trial error variance
SP$setVarE(h2=0.7, H2 = NULL, varE = NULL, corE = NULL)
SP$varE


#manTraits#make Pop object
FounderPop<-newPop(rawPop=FounderPop, simParam = SP)
FounderPop

mean(pheno(FounderPop))
var(pheno(FounderPop))

mean(gv(FounderPop))
var(gv(FounderPop))

#select breeding crosses
NPopCrosses<-100
CandidateCrosses<-t(combn(rownames(pullSegSiteGeno(FounderPop)),2))#might have to switch to explicitly puting ind ids in mappop
PopCrosses<-CandidateCrosses[sample(1:nrow(CandidateCrosses),NPopCrosses),]#random selection
head(PopCrosses)

#make DH donor plants
F1DHdonors<-makeCross(pop=FounderPop, crossPlan=PopCrosses, nProgeny = 1, simParam = SP)
F1DHdonors

#make DH families
DHpop<-makeDH(F1DHdonors, nDH = 20, useFemale = TRUE, keepParents = TRUE, simParam = NULL)
DHpop

#make hybrids
#NHybCrosses<-100 #not needed if there is no selection of hybrids
CandidateP1<-rownames(pullSegSiteGeno(DHpop))
Np1<-20
P1<-sample(CandidateP1,Np1)#random selection
CandidateP2<-rownames(pullSegSiteGeno(DHpop))
CandidateP2<-CandidateP2[!CandidateP2%in%P1]
Np2<-20
P2<-sample(CandidateP2,Np2)#random selection
HybCrosses<-as.matrix(expand.grid(P1,P2))#full factorial
head(HybCrosses)

#option1
HybridPop<-makeCross(pop=DHpop,crossPlan = HybCrosses,nProgeny = 10,simParam = SP)#Be careful with number of simulated plants. It should agree with how h2 is defined.
#option2
makeCross2(females=DHpop[P1], males=DHpop[P2], crossPlan=HybCrosses, nProgeny = 10, simParam = SP)
#option3
hybridCross(females=DHpop[P1], males=DHpop[P2], crossPlan="testcross", returnHybridPop = F, simParam = SP)

hybridCross(females=DHpop[P1], males=DHpop[P2], crossPlan="testcross", returnHybridPop = T, simParam = SP)


#see phenotypes
pheno(HybridPop)
dim(pheno(HybridPop))
dimnames(pheno(HybridPop))
mean(pheno(HybridPop))
var(pheno(HybridPop))

getPed(HybridPop)

View(genParam(HybridPop))
mean(gv(HybridPop))
var(gv(HybridPop))


?selectCross #wrapper around randcross
?randCross#wrapper around makeCross

?selectInd
?selectInd
?selIndex
?smithHazel
?usefulness

#################
?newMapPop #to import founder genotypic data
?importHaplo #wrapper around newMapPop
?importInbredGeno #012 geno matrix

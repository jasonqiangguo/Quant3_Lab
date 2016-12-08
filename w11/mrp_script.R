rm(list=ls(all=TRUE))

library("arm")
library("foreign")

#read in megapoll and attach
marriage.data <- read.dta("gay_marriage_megapoll.dta", convert.underscore = TRUE) 

#read in state-level dataset

Statelevel <- read.dta("state_level_update.dta",convert.underscore = TRUE)
Statelevel <- Statelevel[order(Statelevel$sstate.initnum),]

#read in Census data
Census <- read.dta("poststratification 2000.dta",convert.underscore = TRUE)
Census <- Census[order(Census$cstate),]
Census$cstate.initnum <-  match(Census$cstate, Statelevel$sstate)

#Create index variables

  #At level of megapoll

marriage.data$race.female <- (marriage.data$female *3) + marriage.data$race.wbh# from 1 for white males to 6 for hispanic females
marriage.data$age.edu.cat <- 4 * (marriage.data$age.cat -1) + marriage.data$edu.cat# from 1 for 18-29 with low edu to 16 for 65+ with high edu
marriage.data$p.evang.full <- Statelevel$p.evang[marriage.data$state.initnum]# proportion of evangelicals in respondent's state
marriage.data$p.mormon.full <-Statelevel$p.mormon[marriage.data$state.initnum]# proportion of mormon's in respondent's state
marriage.data$p.relig.full <- marriage.data$p.evang.full + marriage.data$p.mormon.full# combined evangelical + mormom proportions
marriage.data$p.kerry.full <- Statelevel$kerry.04[marriage.data$state.initnum]# kerry's % of 2-party vote in respondent's state in 2004

  #At census level (same coding as above for all variables)

Census$crace.female <- (Census$cfemale *3) + Census$crace.WBH 
Census$cage.edu.cat <- 4 * (Census$cage.cat -1) + Census$cedu.cat 
Census$cp.evang.full<-  Statelevel$p.evang[Census$cstate.initnum]
Census$cp.mormon.full <- Statelevel$p.mormon[Census$cstate.initnum]
Census$cp.relig.full <- Census$cp.evang.full + Census$cp.mormon.full
Census$cp.kerry.full <-  Statelevel$kerry.04[Census$cstate.initnum]


#run individual-level opinion model

individual.model <- glmer(formula = yes.of.all ~ (1|race.female) + (1|age.cat) 
  + (1|edu.cat) + (1|age.edu.cat) + (1|state) + (1|region) + (1|poll) + p.relig.full 
        + p.kerry.full,data=marriage.data, family=binomial(link="logit"))
display(individual.model)

#examine random effects and standard errors for race-female
ranef(individual.model)$race.female
se.ranef(individual.model)$race.female

#create vector of state ranefs and then fill in missing ones
state.ranefs <- array(NA,c(51,1))
dimnames(state.ranefs) <- list(c(Statelevel$sstate),"effect")
for(i in Statelevel$sstate){
    state.ranefs[i,1] <- ranef(individual.model)$state[i,1]
}
state.ranefs[,1][is.na(state.ranefs[,1])] <- 0 #set states with missing REs (b/c not in data) to zero


#create a prediction for each cell in Census data
cellpred <- invlogit(fixef(individual.model)["(Intercept)"]
            +ranef(individual.model)$race.female[Census$crace.female,1]
            +ranef(individual.model)$age.cat[Census$cage.cat,1]
            +ranef(individual.model)$edu.cat[Census$cedu.cat,1]
            +ranef(individual.model)$age.edu.cat[Census$cage.edu.cat,1]
            +state.ranefs[Census$cstate,1]
            +ranef(individual.model)$region[Census$cregion,1]   
            +(fixef(individual.model)["p.relig.full"] *Census$cp.relig.full)
            +(fixef(individual.model)["p.kerry.full"] *Census$cp.kerry.full)
               )

#weights the prediction by the freq of cell                                       
cellpredweighted <- cellpred * Census$cpercent.state

#calculates the percent within each state (weighted average of responses)
statepred <- 100* as.vector(tapply(cellpredweighted,Census$cstate,sum))
statepred
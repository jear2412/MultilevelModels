

#Linear Mixed Effects Models


library(lme4)
library(splines)
library(lattice)
library(ggplot2)
library(brms) # Bayesian Mixed Models

#Mixed effects models

#response variable vs covariates 

#ME: at least one of the covariates
#is a categorical variable representing
#experimental or observational units.

#units: subjects/ humans / plots of land

#Covariates are observed at discrete levels. 

# Parameters associated with the particular levels of a covariate
# are called ''effects'' of the levels


# If the set of possible levels of the covariate is fixed and reproducible
# we model the covariate using fixed effects

#If the levels we observe represent a random sample from the 
# set of all possible levels we incorporate them into a random effects


#Fixed and random is a property of of the levels of the categorical 
# covariate than a property of the effects

#Obs: fixed effects are parameters and random effects arent
    # Random effects are random variables (unobserved)


# book: Statistical Methods in Research and Production
#O.L Davies 1947
#Examples of use of random effects to characterize batch to batch
#variability in chemical processes.


#------ Dyestuff data

#Davies Goldsmith 1972 

# an investigation to find out how much the variation from batch to batch
#in the quality of an intermediate product (H-acid) contributes to the 
#variation in the yield of the dyestuff (Naphthalene Black 12B) made from it. 


#Six samples of intermediate , representing different batches of 
# works manufacture, were obtained, and five preparations of the 
# dyestuff were made in the lab from each sample

# The equivalent yield of each preparation 
#as grams of standard colour was determined by dye-trial.

str(Dyestuff)

levels(Dyestuff$Batch)

#30 obs of yield

head(Dyestuff)
tail(Dyestuff)
summary(Dyestuff)

#Data are balanced

set.seed(1234543)
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
              ylab = "Batch", jitter.y = TRUE, pch = 21,
              xlab = "Yield of dyestuff (grams of standard color)",
              type = c("p", "a")))

#Variability within yields, even for preparations from same batch
#Batch to batch variability

#F vs C
#Why is there so much variability bw batches?

#Objective: predict the yield from future batches, taking into account the
#batch-to-batch variability and the within-batch variability.

#Should batches be fixed effect or random effect?



#----- Dyestuff 2

str(Dyestuff2)
summary(Dyestuff2)


print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
              ylab = "Batch", jitter.y = TRUE, pch = 21,
              xlab = "Simulated response (dimensionless)",
              type = c("p", "a")))


#Within variability is high
#between variability is low

#----- Dyestuff 1 model

fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)#REML criterion``
print(fm1)
summary(fm1)



#Variabilities: batch to batch and per observation (within )
              #residual: not explained, estimate of variance of residuals

#intercept: overall mean
# how to interpret: given the intercept you add a random effect with var 1764 
#obs: Std Dev is NOT standard error

#----- Dyestuff 2 model

((fm2 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2)))
summary(fm2)

#Variance==0 
#indicates that the level of bw group variability is not sufficient to warrant
#incorporating random effects in the model

#Estimates of variance can be 0 
#fm2 is then

summary(fm2a <- lm(Yield ~ 1, Dyestuff2))

#obs: Std dev of residual is the same


pr<-profile(fm1)
confint(pr, level=0.95)

ranef(fm1) #BLUPS: ''estimates of random effects''  

print(dotplot(ranef(fm1, condVar=TRUE), strip = FALSE)) #prediction intervals



#-------------  Sleepstudy

#Reaction times to a series of tests after sleep deprivation

str(sleepstudy)
summary(sleepstudy)

(split_plot <- ggplot(aes(Days, Reaction), data = sleepstudy) + 
    geom_point() + 
    facet_wrap(~ Subject) + geom_smooth(method = "lm")+
    xlab("days") + 
    ylab("Reaction score"))




print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(6,3), type = c("g", "p", "r"),
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)"))


xtabs(~ Subject + Days, sleepstudy)

#Random Efect and Slope per subject
M1 <- lmer ( Reaction ~1+ Days + ( 1+Days | Subject ) , data= sleepstudy )

summary(M1)
pr<-profile(M1)
confint(pr, level=0.90)

ranef(M1) #BLUPS: ''estimates of random effects'' 


M2 <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0+Days|Subject), sleepstudy, REML = 0)
summary(M2)
pr<-profile(M2)
confint(pr, level=0.90)

#M2 check correlation fixed

anova(M1,M2) #small chisq extra parameter isnt significant


df <- coef(lmList(Reaction ~ Days | Subject, sleepstudy))
fclow <- subset(df, `(Intercept)` < 251)
fchigh <- subset(df, `(Intercept)` > 251)
cc1 <- as.data.frame(coef(M1)$Subject)
names(cc1) <- c("A", "B")
df <- cbind(df, cc1)
ff <- fixef(M1)


print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(6,3), type = c("g", "p", "r"),
             coef.list = df[,3:4],
             panel = function(..., coef.list) {
               panel.xyplot(...)
               panel.abline(as.numeric(coef.list[packet.number(),]),
                            col.line = trellis.par.get("superpose.line")$col[2],
                            lty = trellis.par.get("superpose.line")$lty[2]
               )
               panel.abline(fixef(M2),
                            col.line = trellis.par.get("superpose.line")$col[4],
                            lty = trellis.par.get("superpose.line")$lty[4]
               )
             },
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)",
             key = list(space = "top", columns = 3,
                        text = list(c("Within-subject", "Mixed model", "Population")),
                        lines = list(col = trellis.par.get("superpose.line")$col[c(1:2,4)],
                                     lty = trellis.par.get("superpose.line")$lty[c(1:2,4)]))))


# Bayesian LMM
#Stan HMC

brms1 <- brm( Reaction ~1+ Days + ( 1+Days | Subject ) ,
             family=gaussian,
              data= sleepstudy,iter=10000, chains=2, cores=1  )

brms1$prior

print(brms1)


summary(M1)
summary(brms1)

plot(brms1)
launch_shinystan(brms1) 





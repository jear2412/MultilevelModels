#Radon Example

library(mvtnorm)
library(purrr)
library(tidyr)
library(stringr)
library(janitor)
library(tidybayes)
library(forcats)
library(rstan) 
options(mc.cores = parallel::detectCores())
options(buildtools.check = function(action) TRUE )
#rstan_options(auto_write = TRUE)
library("brms")
library("ggplot2")
library("dplyr")
library("loo")
library(nlme)
library(lme4)
library(GGally)
library(corrplot)
library(tidyverse)
library(bayesplot)
library("ggridges")
library(cmdstanr)
library(HLMdiag)
library(loo)

options(mc.cores = parallel::detectCores())

setwd("~/Documents/GitHub/R2D2M2/scripts")

#------ CDMSTAN

PATH_TO_CMDSTAN= "~/cmdstan"
set_cmdstan_path(PATH_TO_CMDSTAN)
cmdstan_path()
cmdstan_version()

#---- functions

Convert_Numeric = function(X) {
  L = levels(X)
  Y = as.numeric(factor(X, labels = seq(1:length(L))))
  return(Y)
}

#----- 

data("radon", package = "HLMdiag")
head(radon)

radon$basement = as.factor(radon$basement)
radon$county = as.factor(radon$county)

# log radon histogram
ggplot( radon,    mapping= aes( log.radon)) +
  geom_histogram(color="black", fill = "white") + ggtitle("logRadon histogram") + theme_minimal()

#per basement
radon %>%
  ggplot(aes(x = basement, y = log.radon)) +
  geom_point( size=0.5)

radon %>% 
  ggplot(aes(x = basement, y = log.radon)) +
  geom_boxplot()

radon %>% 
  ggplot(aes(x = log.radon, y = basement)) +
  geom_density_ridges( )


#per county
radon %>%
  ggplot(aes(x = county, y = log.radon)) +
  geom_point( size=0.5)

radon %>% 
  ggplot(aes(x = county, y = log.radon)) +
  geom_boxplot()

radon %>% 
  ggplot(aes(x = log.radon, y = county)) +
  geom_density_ridges( )

#--------- Model 0
K=1 # County is the grouping factor
D=3
N= dim(radon)[1]

df= radon %>% 
  select(log.radon , basement, uranium , county ) %>% 
  mutate( log.radon = as.numeric(log.radon))  %>% 
  mutate( basement = as.numeric(basement))  %>% 
  mutate( county = as.numeric(county))

#floor uranium
df$fur = df$basement* df$uranium

df$county =as.factor(df$county)

Lcounty= nlevels(as.factor(radon$county)) 

#Scale X and Z

X = df %>% select(basement, uranium, fur) 

X= scale(X, center= TRUE, scale=FALSE)

X = as_tibble(X) %>%
  add_column( x0= rep(1, N)   ,.before = 1  ) 

head(X)

Jcounty= df %>% 
  select(county) %>% 
  mutate( county = Convert_Numeric( as.factor(radon$county)   )) %>% 
  mutate_if(is.factor, as.numeric)

D_1 =1
R2D2_alpha= rep(1, D+1+D_1 )

#check Z

datM1 <- list( N=N, Y= df$log.radon  , D= D+1,X=  X , K=K,
               L_1= Lcounty, D_1= D_1, J_1= Jcounty$county  , 
               Z_1= t(X[c(1,2)]),
               R2D2_alpha = R2D2_alpha,
               R2D2_mean_R2=0.5, R2D2_prec_R2= 1, prior_only= 0)

#cmdstan 
file <- file.path("~/Documents/GitHub/R2D2M2/stan", "radonv0.stan") 
mods <- cmdstan_model(file)

fit <- mods$sample(
  data = datM1,
  chains = 1,
  refresh = 500,
  iter_sampling=3000,
  seed=1
)



fit_rstan <- rstan::read_stan_csv(fit$output_files())

r2 <-fit$draws("R2D2_R2")
tau2 <-fit$draws("R2D2_tau2")

sigma <-  fit$draws("sigma")
bs <-  fit$draws( "b")

r1s <-  fit$draws( "r_1")

phis <- fit$draws("R2D2_phi")

nbs <-  dimnames(bs)$variable
nrs <-  dimnames(r1s)$variable # [ group, d, level ]

#---- plots

p1 <- mcmc_hist(r2) 

p2 <- mcmc_intervals( bs,  outer_size = 1, inner_size=4, point_size = 6)+hline_0()

p3 <- mcmc_hist( sigma)

p4 <-  mcmc_intervals(phis  )

p1
p2
p3
p4

r1sm=as.data.frame(r1s)
pr1 <- mcmc_intervals( r1s)

mcmc_areas(r1sm[,1:20] )
mcmc_areas(bs)

M1 <-  lmer(  data=df,  log.radon ~ basement + uranium + fur + (1 + basement|| county)     )
M2 <-  brm( data=df,  log.radon ~ basement + uranium + fur + (1 + basement|| county), backend= "cmdstan")

summary(M1)
summary(M2)

s <-  summary(fit_rstan,c("Intercept" ,"b", "R2D2_R2","R2D2_tau2" , "sigma", "R2D2_phi" ), probs=c(0.025, 0.975) )
s$summary

#---- Posterior predictive

y= radon$log.radon
ytilde <- fit$draws("y_tilde")
ytilde= as.matrix(as.data.frame(ytilde))

I=  sample(1:  nrow(ytilde) , 100,  replace=FALSE) 

ppc1 <- ppc_dens_overlay(y, ytilde[I,])
ppc2 <-ppc_intervals(y, ytilde[I,])
ppc3 <-ppc_intervals(y, ytilde[I,],x = y) + abline_01()

ppc1
ppc2


loglik1 <- fit$draws("log_lik")
r_eff <- relative_eff(exp(loglik1), cores = 8)

loo_1 <- loo(loglik1,  r_eff = r_eff, cores = 8)
print(loo_1)
plot(loo_1)

#--------- Model 1

K=2
D= 2
N= dim(radon)[1]

data("radon", package = "HLMdiag")
head(radon)

radon$basement = as.factor(radon$basement)
radon$county = as.factor(radon$county)

Lbas= nlevels(as.factor(radon$basement)) 
Lcounty= nlevels(as.factor(radon$county)) 

df= radon %>% 
  select(log.radon , basement, uranium , county ) %>% 
  mutate( log.radon = as.numeric(log.radon))  %>% 
  mutate( basement = as.numeric(basement))  %>% 
  mutate( county = as.numeric(county))

df$basement= df$basement-1


#Scale X and Z

X = df %>% 
  select(uranium) 

X= scale(X, center= TRUE, scale=FALSE)


X = as_tibble(X) %>%
    add_column( basement= df$basement ) %>%
    add_column( x0= rep(1, N)   ,.before = 1  ) 

head(X)


Jbas= df %>% 
  select(basement) %>% 
  mutate( basement = Convert_Numeric( as.factor( radon$basement  )   )) %>% 
  mutate_if(is.factor, as.numeric)

Jcounty= df %>% 
  select(county) %>% 
  mutate( county = Convert_Numeric( as.factor(radon$county)   )) %>% 
  mutate_if(is.factor, as.numeric)

R2D2_alpha= rep(1, (D)+K+(D)*K)

datM1 <- list( N=N, Y= df$log.radon  , D= D+1, Dnum=1,  X=  X , K=K,
                L_1= Lbas, D_1= D+1, J_1= Jbas$basement  , 
                L_2= Lcounty, D_2= D+1, J_2= Jcounty$county,
                Z_1= t(X),Z_2= t(X),
                R2D2_alpha = R2D2_alpha,
                R2D2_mean_R2=0.5, R2D2_prec_R2= 1, prior_only= 0)

#cmdstan 
file <- file.path("~/Documents/GitHub/R2D2M2/stan", "radonv1.stan") 
mods <- cmdstan_model(file)

mods$print()
mods$exe_file()

fit <- mods$sample(
  data = datM1,
  chains = 1,
  refresh = 500,
  iter_sampling=4000,
  adapt_delta = 0.95,
  seed=1
)

fit_rstan <- rstan::read_stan_csv(fit$output_files())

r2 <-fit$draws("R2D2_R2")
tau2 <-fit$draws("R2D2_tau2")

sigma <-  fit$draws("sigma")
bs <-  fit$draws( "b")

r1s <-  fit$draws( "r_1")
r2s <-  fit$draws( "r_2")
phis <- fit$draws("R2D2_phi")

nbs <-  dimnames(bs)$variable
#nrs <-  dimnames(rs)$variable # [ group, d, level ]

#---- plots

p1 <- mcmc_hist(r2) 

p2 <- mcmc_intervals( bs,  outer_size = 1, inner_size=4, point_size = 6)+hline_0()

p3 <- mcmc_hist( sigma)

p4 <-  mcmc_intervals(phis  )

gridname <- "Radon"
grid1 <- bayesplot_grid(p1, p2, p3, p4,
                       subtitles = rep( gridname , 4))
grid1

pr1 <- mcmc_intervals( r1s)
pr2 <- mcmc_intervals( r2s)



#---- Posterior predictive

y= radon$log.radon
ytilde <- fit$draws("y_tilde")
ytilde= as.matrix(as.data.frame(ytilde))

I=  sample(1:  nrow(ytilde) , 100,  replace=FALSE) 

ppc1 <- ppc_dens_overlay(y, ytilde[I,])
ppc2 <-ppc_intervals(y, ytilde[I,])
ppc3 <-ppc_intervals(y, ytilde[I,],x = y) + abline_01()

#per subject
Js= as.factor(as.matrix(Jbas ))
ppc4 <-ppc_scatter_avg_grouped(y,ytilde[I,],Js,alpha = 0.7)+ abline_01()

gridname <-paste("Predictive Obs") 
grid2 <- bayesplot_grid(ppc1, ppc2, ppc3, ppc4 ,
                       subtitles = rep( gridname , 4))

grid2
#------ loocv

loglik1 <- fit$draws("log_lik")
r_eff <- relative_eff(exp(loglik1), cores = 8)

# preferably use more than 2 cores (as many cores as possible)
# will use value of 'mc.cores' option if cores is not specified
loo_1 <- loo(loglik1,  r_eff = r_eff, cores = 8)
print(loo_1)
plot(loo_1)

#---------------------- Model 2

data("radon", package = "HLMdiag")
head(radon)

radon$basement = as.factor(radon$basement)
radon$county = as.factor(radon$county)

Lbas= nlevels(as.factor(radon$basement)) 
Lcounty= nlevels(as.factor(radon$county)) 

df= radon %>% 
  select(log.radon , basement, uranium , county ) %>% 
  mutate( log.radon = as.numeric(log.radon))  %>% 
  mutate( basement = as.numeric(basement))  %>% 
  mutate( county = as.numeric(county))

df$basement= df$basement-1

#Scale X and Z
Xc = df %>% 
  select(basement) 

head(X)

Jcounty= df %>% 
  select(county) %>% 
  mutate( county = Convert_Numeric( as.factor(radon$county)   )) %>% 
  mutate_if(is.factor, as.numeric)

R2D2_alpha= rep(0.5, (D)+K+(D)*K)

means_X=c(1)
sds_X= c(1, sd(df$uranium ) , sd(df$uranium*df$basement) )


K= 2
D= 1
N= dim(radon)[1]


R2D2_alpha=rep(1, 3)
datM2 <- list( N=N, Y= df$log.radon  , D= D+1, Dnum=1,   K=K,
               Xc= Xc, sds_X = sds_X,
               L_1= Lcounty, D_1= 1, J_1= Jcounty$county  , 
               L_2= Lcounty, D_2= 1, J_2= Jcounty$county,
               Z_1= df$uranium  ,Z_2= df$uranium*df$basement ,
               R2D2_alpha = R2D2_alpha,
               R2D2_mean_R2=0.5, R2D2_prec_R2= 1, prior_only= 0)

#cmdstan 
file <- file.path("~/Documents/GitHub/R2D2M2/stan", "radonv2.stan") 
mods <- cmdstan_model(file)

mods$print()
mods$exe_file()

fit2 <- mods$sample(
  data = datM2,
  chains = 1,
  refresh = 500,
  iter_sampling=4000,
  adapt_delta = 0.95,
  seed=1
)

fit_rstan2 <- rstan::read_stan_csv(fit2$output_files())

r22 <-fit2$draws("R2D2_R2")
tau22 <-fit2$draws("R2D2_tau2")

sigma2 <-  fit2$draws("sigma")
bs2 <-  fit2$draws( "b")

r1s2 <-  fit2$draws( "r_1")
r2s2 <-  fit2$draws( "r_2")
phis2 <- fit2$draws("R2D2_phi")

nbs2 <-  dimnames(bs)$variable
#nrs <-  dimnames(rs)$variable # [ group, d, level ]

#---- plots

p1_2 <- mcmc_hist(r22) #expected, should we make a more restrictive R2 prior?

p2_2 <- mcmc_intervals( bs2,  outer_size = 1, inner_size=4, point_size = 6)+hline_0()

p3_2 <- mcmc_hist( sigma2)

p4_2 <-  mcmc_intervals(phis2  )

gridname <- "Radon 2"
grid_2 <- bayesplot_grid(p1_2, p2_2, p3_2, p4_2,
                       subtitles = rep( gridname , 4))
grid_2

pr1_2 <- mcmc_intervals( r1s2)
pr2_2 <- mcmc_intervals( r2s2)



#---- Posterior predictive

y= radon$log.radon
ytilde2 <- fit2$draws("y_tilde")
ytilde2= as.matrix(as.data.frame(ytilde2))

I=  sample(1:  nrow(ytilde2) , 100,  replace=FALSE) 

ppc1_2 <- ppc_dens_overlay(y, ytilde2[I,])
ppc2_2 <-ppc_intervals(y, ytilde2[I,])
ppc3_2 <-ppc_intervals(y, ytilde2[I,],x = y) + abline_01()

#per subject

gridname <-paste("Predictive Obs") 
grid_ppc2 <- bayesplot_grid(ppc1_2, ppc2_2, ppc3_2,
                       subtitles = rep( gridname , 3))

#------ loocv

loglik1_2 <- fit2$draws("log_lik")
r_eff_2 <- relative_eff(exp(loglik1_2), cores = 8)

# preferably use more than 2 cores (as many cores as possible)
# will use value of 'mc.cores' option if cores is not specified
loo_1_2 <- loo(loglik1_2,  r_eff = r_eff_2, cores = 8)
print(loo_1_2)

plot(loo_1_2)

#------ Summary

# Model 1
grid1
grid2

# basement, uranium 
# 2 grouping factors: basement(0,1), county (85 levels)

# more reasonable basement as fixed 
pr1
pr2 #uranium per county

# Model 2

grid_2
grid_ppc2


pr1_2 # uranium per county
pr2_2 # urannium per county * basement


# loo 
print(loo_1)
print(loo_1_2)

plot(loo_1)
plot(loo_1_2)

loo_compare(loo_1, loo_1_2)

# 





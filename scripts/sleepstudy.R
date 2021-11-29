

# Sleepstudy example

library(lme4)
library(tidyverse)
library(ggplot2)
library(brms) # Bayesian Mixed Models
library(cmdstanr)

# credit for plots goes to 
# # https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
# I did minor edits on them 

#---- cmdstan 
#ignore if you dont have cmdstan installed

options(mc.cores = parallel::detectCores())
PATH_TO_CMDSTAN= "~/cmdstan"

set_cmdstan_path(PATH_TO_CMDSTAN)
cmdstan_path()
cmdstan_version()


#-------------  Sleepstudy

#Reaction times (ms) to a series of tests after sleep deprivation.
# subjects sleep max 3 h for a number of days

str(sleepstudy)
summary(sleepstudy)

# Convert to tibble for better printing. Convert factors to strings
sleepstudy <- sleepstudy %>% 
  as_tibble() %>% 
  mutate(Subject = as.character(Subject))

df_sleep<- sleepstudy
xlab <- "Days of sleep deprivation"
ylab <- "Average reaction time (ms)"


#------- Complete Pooling

#--- plot of data
p1<- ggplot(aes(Days, Reaction), data = sleepstudy) + 
  geom_point(size=1) + 
  labs(x = xlab, y = ylab)+ ggtitle("Sleep study") + 
  scale_x_continuous(breaks = 0:4 * 2)+theme_bw()

#--- plot of data and complete pooling regression

p2<-ggplot(aes(Days, Reaction), data = sleepstudy) + 
  geom_point(size=1) + 
  labs(x = xlab, y = ylab)+ ggtitle("Sleep study") + 
  geom_smooth(method = "lm", color="blue") +
  scale_x_continuous(breaks = 0:4 * 2)+theme_bw()


#------ no pooling
# unaware other participants exist
# separate line for each cluster of data

p3<-ggplot(aes(Days, Reaction), data = sleepstudy) + 
  geom_point(size=1) + 
  facet_wrap("Subject")+ 
  labs(x = xlab, y = ylab)+ 
  scale_x_continuous(breaks = 0:4 * 2)+ggtitle("Sleep study: no pooling")+theme_bw()


p4<-ggplot(df_sleep) + 
  aes(x = Days, y = Reaction) + 
  stat_smooth(method = "lm", se = TRUE ) +
  # Put the points on top of lines
  geom_point() +
  facet_wrap("Subject") +
  labs(x = xlab, y = ylab)+ 
  scale_x_continuous(breaks = 0:4 * 2)+ggtitle("Sleep study: no pooling")+theme_bw()

df_no_pooling <- lmList(Reaction ~ Days | Subject, df_sleep) %>% 
  coef() %>% 
  # Subject IDs are stored as row-names. Make them an explicit column
  rownames_to_column("Subject") %>% 
  rename(Intercept = `(Intercept)`, Slope_Days = Days) %>% 
  add_column(Model = "No pooling") 
  

#---- Complete vs No pooling

m_pooled <- lm(Reaction ~ Days, df_sleep) 

# Repeat the intercept and slope terms for each participant
df_pooled <- tibble(
  Model = "Complete pooling",
  Subject = unique(df_sleep$Subject),
  Intercept = coef(m_pooled)[1], 
  Slope_Days = coef(m_pooled)[2]
)

df_models <- bind_rows(df_pooled, df_no_pooling) %>% 
  left_join(df_sleep, by = "Subject")


model_comparison <- ggplot(df_models) + 
  aes(x = Days, y = Reaction) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_Days, color = Model),
    size = .75
  ) + 
  geom_point() +
  facet_wrap("Subject") +
  labs(x = xlab, y = ylab) + 
  scale_x_continuous(breaks = 0:4 * 2) + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + ggtitle("Sleep study: No pooling")+theme_bw()+
  theme(legend.position = c(0.75, 0.1), legend.title = element_blank())


p5<- model_comparison

#---------------- Partial Pooling

m <- lmer(Reaction ~ 1 + Days + (1 + Days | Subject), df_sleep)

df_partial_pooling <- coef(m)[["Subject"]] %>% 
  rownames_to_column("Subject") %>% 
  as_tibble() %>% 
  rename(Intercept = `(Intercept)`, Slope_Days = Days) %>% 
  add_column(Model = "Partial pooling")

df_models <- bind_rows(df_pooled, df_no_pooling, df_partial_pooling) %>% 
  left_join(df_sleep, by = "Subject")


model_comparison <- ggplot(df_models) + 
  aes(x = Days, y = Reaction) + 
  geom_abline(
    aes(intercept = Intercept, slope = Slope_Days, color = Model),
    size = 0.75
  ) + 
  geom_point() +
  facet_wrap("Subject") +
  labs(x = xlab, y = ylab) + 
  scale_x_continuous(breaks = 0:4 * 2) + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + ggtitle("Sleep study: Pooling")+theme_bw()+
  theme(legend.position = c(0.75, 0.1), legend.title = element_blank())
  

p6<- model_comparison


file <- file.path("~/Documents/GitHub/R2D2M2/talk/plots/")
paste()

ggsave(p1,filename= paste(file,"1.pdf", sep='') ,device="pdf")
ggsave(p2,filename= paste(file,"2.pdf", sep='') ,device="pdf")
ggsave(p3,filename= paste(file,"3.pdf", sep='') ,device="pdf")
ggsave(p4,filename= paste(file,"4.pdf", sep='') ,device="pdf")
ggsave(p5,filename= paste(file,"5.pdf", sep='') ,device="pdf")
ggsave(p6,filename= paste(file,"6.pdf", sep='') ,device="pdf")


#---- Bivariate shrinkage

# Also visualize the point for the fixed effects
df_fixef <- tibble(
  Model = "Partial pooling (average)",
  Intercept = fixef(m)[1],
  Slope_Days = fixef(m)[2]
)

# Complete pooling / fixed effects are center of gravity in the plot
df_gravity <- df_pooled %>% 
  distinct(Model, Intercept, Slope_Days) %>% 
  bind_rows(df_fixef)

df_pulled <- bind_rows(df_no_pooling, df_partial_pooling)

p7<-ggplot(df_pulled) + 
  aes(x = Intercept, y = Slope_Days, color = Model, shape = Model) + 
  geom_point(size = 2.5) + 
  geom_point(
    data = df_gravity, 
    size = 5,
    # Prevent size-5 point from showing in legend keys
    show.legend = FALSE
  ) + 
  # Draw an arrow connecting the observations between models
  geom_path(
    aes(group = Subject, color = NULL), 
    arrow = arrow(length = unit(.02, "npc")),
    show.legend = FALSE
  ) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = Subject, color = NULL), 
    data = df_no_pooling,
    show.legend = FALSE
  ) + 
  theme(
    legend.position = "bottom", 
    legend.justification = "right"
  ) + 
  ggtitle("Pooling of regression parameters") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_brewer(palette = "Dark2") 

ggsave(p7,filename= paste(file,"7.pdf", sep='') ,device="pdf")




#---------- Likelihood modeling

#Random Efect and Slope per subject
M1 <- lmer ( Reaction ~ 1+ Days + ( 1+Days | Subject ) , data= sleepstudy ) 

summary(M1)
pr<-profile(M1)
confint(pr, level=0.90)
ranef(M1) #BLUPS: ''estimates of random effects'' 

M2 <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0+Days|Subject), sleepstudy, REML = 0)
summary(M2)
pr<-profile(M2)
confint(pr, level=0.90)
ranef(M1) #BLUPS: ''estimates of random effects'' 

#M2 check correlation fixed

anova(M1,M2) #small chisq extra parameter isnt significant

#---------- Bayesian LMM

#Stan HMC

brms1 <- brm( Reaction ~1+ Days + ( 1+Days | Subject ) ,
              family=gaussian,
              data= sleepstudy,iter=4000, chains=2, cores=8, backend = "cmdstan"  )

brms1$prior

print(brms1)


summary(M1)
summary(brms1)

plot(brms1)
launch_shinystan(brms1) 





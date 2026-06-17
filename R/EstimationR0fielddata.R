########################################################################################################
#                                                                                                      #
#                     Methods to estimate R0 from field data                                           #
#                     Course Advanced Veterinary Epidemiology                                          #
#                                                                                                      #
#                     Version March 2023                                                               #
#                     Copyright Egil Fischer e.a.j.fischer@uu.nl                                       #
#                                                                                                      #
########################################################################################################

#Load required libraries if packages are not installed use install.packages("name of package") 
library(graphics)
library(bbmle)
library(ggplot2)

#set working directory to the location of the CSV-files
#either change the path in the next line or use menu: Session > Set Working Directory...
setwd("C:/Surfdrive/Onderwijs/Cursorisch/MSc Vet Epi/2019/EpidemiologyOfAnimalInfectiousDiseases/Estimation from field data")

########################################################################################################
#                                                                                                      #
#                     Maximum likelihood estimation of exponential phase                               #
#                                                                                                      #
#                                                                                                      #
########################################################################################################

#########################################################################################################
# This is a resampled dataset of part of the SARS epidemic in Hong Kong. Each record gives the number   #
# of new cases per day since start of the epidemic. The dataset contains only cases of the first 55 days#
# of the epidemic.                                                                                      #
#########################################################################################################

#import data
rm(sars)
sars <- read.csv("SARS.csv", header = TRUE, sep = ",")[,1:2]

#descriptive data
ggplot(data = sars,aes(x = Time, y = Incidence))+geom_point()

#Estimation of R0 using the exponential phase of the epidemic
end.exp.phase <- 30 #decide when the exponential phase ends

#log-likelihood function for the selected data
#maximize the log-likelihood
#We use the suppressWarnings to avoid warnings when the algorithm tries to go outside the possible range of values
#the log-likelihood function assumes 
#1. Poisson distributed counts around an expected value
#2. An exponential growth of the number of cases
dLL <-  function(y0, r,  log = TRUE) {
  -sum(suppressWarnings(log(dpois(x = sars$Incidence[sars$Time<= end.exp.phase],
                                  lambda = y0 * exp(r *  sars$Time[sars$Time<= end.exp.phase]) , log = FALSE))))
  }
fit <- mle2(minuslogl = dLL, start = list(y0 = 1, r = 0.1))
fit

#get the estimated parameters
y0.est <- coef(fit)[1]
r.est <-coef(fit)[2]

#also possible to use generalized linear models 
fit.glm <- glm(formula = Incidence ~ Time, data = sars[sars$Time<= end.exp.phase,],family = "poisson")

#the GLM used the formulate ln(y) ~ln(y0 e^rt), which is ln(y) ~ ln(y0) + r t. Interpretation of the coefficients is thus:
#estimates differ slightly due to numerical differences in the methods.
y0.est <- exp(coef(fit.glm)[1])
r.est <- coef(fit.glm)[2]

#Plot data and fitted model#
#data
sars.plot <- ggplot()+
  geom_point(data = sars,aes(x = Time, y = log(Incidence)))+xlab ("Time")+ ylab("log Cases")
#show plot
sars.plot

#add fitted model in exponential phase
sars.plot+ 
  geom_line(data= data.frame(Time = seq(1,end.exp.phase,1), 
                 Incidence = sapply(X = seq(1,end.exp.phase,1),function(t){log(y0.est*exp(r.est * t))})),
            aes(Time,Incidence), 
  linetype = 1,linewidth = 2, colour = "red")
#add extrapolation after the exponential phase
sars.plot+ 
   geom_line(data= data.frame(Time = seq(1,end.exp.phase,1), 
                              Incidence = sapply(X = seq(1,end.exp.phase,1),function(t){log(y0.est*exp(r.est * t))})),
             aes(Time,Incidence), 
             linetype = 1,linewidth = 2, colour = "red")+
  geom_line(data= data.frame(Time = seq(end.exp.phase,55,1), 
                             Incidence = sapply(X = seq(end.exp.phase,55,1),function(t){log(y0.est*exp(r.est * t))})),
            aes(Time,Incidence), 
            linetype = 1,linewidth = 2, colour = "pink")

#estimate R0
Tg <- 16.
#R0 
1 + r.est*Tg



########################################################################################################
#                                                                                                      #
#                     Final Size                                                                       #
#                                                                                                      #
########################################################################################################

#########################################################################################################
# This is a dataset of five outbreaks of equine influenza A virus (H3N8) in Japan, 1971.                #
# For each race course, the total number of horses is recorded, as well as the number of horses         # 
# infected during the outbreak (final size)                                                             #         #
#########################################################################################################
race <- read.csv("RacingHorses.csv", header = TRUE, sep = ",")
#plot-data
ggplot()+
  geom_point(aes(y = (race$Infected/race$Population), x = factor(race$Race.course)))+
  xlab( "Race course")+ ylab("Final Size")+ylim(c(0.6,1))

#Calculate Final Size and the standard deviation of the Final Size
race$FS <- (race$Infected/race$Population) #calculate the final size
race$sdFS <- sqrt(race$FS*(1-race$FS)/race$Population)  #calculate the standard deviation

#plot the final size
ggplot(race)+
  geom_point(aes(x = Race.course, y = FS))+ #mean
  geom_errorbar(aes(x = Race.course,ymin = FS-1.96 * sdFS,ymax = FS+1.96 * sdFS))+ #confidence interval
  xlab("Race course")+ylab("Final Size")+ #axis labels
  ylim(c(0,1.))


#Calculate R0
race$R0 <- -log(1-race$FS)/race$FS
race$R0.ul<- -log(1-race$FS - 1.96 * race$sdFS)/(race$FS - 1.96 * race$sdFS) #upper limit of R0 corresponds to the smallest fraction susceptibles
race$R0.ll<- -log(1-race$FS + 1.96 * race$sdFS)/(race$FS + 1.96 * race$sdFS) #lower limit of R0 corresponds to the largest fraction susceptibles

#plot the basic reproduction ratio
ggplot(race)+
  geom_point(aes(x = Race.course, y = R0))+ #mean
  geom_errorbar(aes(x = Race.course,ymin = R0.ll,ymax = R0.ul))+ #confidence interval
  ylim(c(0,5))

########################################################################################################
#                                                                                                      #
#                     Endemic equilibrium                                                              #
#                                                                                                      #
########################################################################################################

#########################################################################################################
# This is a dataset of Salmonella Dublin seroprevalence at two Dutch dairy herds                        #
# (Van Schaik et al 2007: Vet Res 38, pp 861-9). The dataset contains the sampling day                  #
# (since the first sample), the herd size at that time, and the number of seropositive animals.         #
#########################################################################################################
salmo <- read.csv("Salmonella.csv", header = TRUE, sep = ",")
attach(salmo)
salmo.plot <- ggplot() + 
  geom_path(data = salmo, aes(x = Day, y =Seropositive/Herd.size, colour = factor(Herd)))+
  geom_point(data = salmo,aes(x = Day, y =Seropositive/Herd.size, colour = factor(Herd)))+
  scale_colour_manual(name = "Herd", values = c("red","blue"))

#determine the beginning of the endemic equilibrium 
start.end.eq.herd1 <- 13
start.end.eq.herd2 <- 50

#estimate the endemic equilibrium
end.eq.herd1 <- mean(Seropositive[Herd == 1 & Day >= start.end.eq.herd1 ]/Herd.size[Herd == 1 & Day >= start.end.eq.herd1])
end.eq.herd2 <- mean(Seropositive[Herd == 2 & Day >= start.end.eq.herd2 ]/Herd.size[Herd == 2 & Day >= start.end.eq.herd2])

#plot the endemic equilibria to the graph
salmo.plot + 
  geom_line(aes(y = c(end.eq.herd1,end.eq.herd1), x = c(start.end.eq.herd1, max(Day[Herd ==1]) + 1)), linetype = 2, colour = "red")+
  geom_line(aes(y = c(end.eq.herd2,end.eq.herd2), x = c(start.end.eq.herd2, max(Day[Herd ==2]) + 1)), linetype = 2, colour = "blue")

#literature values for other parameters
lat<- 2
inf<- 10
imm<- 90

#Estimate R0
R0.end.herd1 <- 1/(1-(1+(lat+inf)/imm)*end.eq.herd1)
R0.end.herd1
R0.end.herd2 <- 1/(1-(1+(lat+inf)/imm)*end.eq.herd2)
R0.end.herd2 
detach(salmo)



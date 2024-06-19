###############################################################################################
#                                                                                             #
#                              load packages                                                  #
#                                                                                             #
###############################################################################################

library(deSolve)
library(rootSolve)
###############################################################################################
#                                                                                             #
#                         SIR -     Deterministic without demography                          #
#                                                                                             #
###############################################################################################

#Ordinary differential Equations
sir <- function(Time, State, Pars){
  with(as.list(c(State,Pars)),{
    dx <- -1 * beta * X * Y / sum(State)
    dy <- beta * X * Y / sum(State)  - gamma * Y
    dz <- gamma * Y
    return(list(c(dx,dy,dz)))
  })
}

#Parameters
pars <- c(beta = 0.6, #beta = transmission constant
          gamma = 0.2 #gamma = recovery constant
)
#initial values
N0 = 10000 #population size
Y0 = .1 #initial number infected 
Z0 = 0# initial number recovered, 
X0 = N0 - Y0 - Z0
init <- c(X = X0,Y = Y0,Z = Z0)


#Data storage time
dt = 0.5#timestep for storing data
sim_duration  = 150 #length of the simulation
times <- seq(0, sim_duration, by = dt)


#Solve the ordinary differential equations
ode.out <- ode(init, times, sir, pars) 

#plot X, Y, and Z against time
plot(x = ode.out[,1], y = ode.out[,2], type = "l", xlab = "time", ylab = "animals", col = "purple",lwd =2)
lines(x = ode.out[,1], y = ode.out[,3], type = "l", xlab = "time", ylab = "animals", col = "red",lwd =2)
lines(x = ode.out[,1], y = ode.out[,4], type = "l", xlab = "time", ylab = "animals", col = "palegreen1",lwd =2)
legend(x = 60, y = 15000, c("X","Y","Z"),lwd =2, col = c("purple", "red", "palegreen1"))


#plot  X against Y
plot(x = ode.out[,2], y = ode.out[,3], type = "l", xlab ="X", ylab ="Y")
s <- seq(length(ode.out[,2])-1)
arrows(x0 = ode.out[s,2], y0 = ode.out[s,3], x1 = ode.out[s + 1,2], y1 = ode.out[s + 1 ,3],length = .075)

#plot final size againt R0.
# For simplicity we set gamma = 1 so that beta = R0
finalsize <- NULL
for(R0 in seq(0,2,0.05)){
  #add the outcomes to the results
finalsize <- rbind(finalsize, data.frame(R0= R0, #R0 value
                                         finalsize = runsteady(y = init, func = sir, parms = c(beta = R0, 
                                          gamma = 1))$y["Z"])) #final size at equilibrium (Y ==0)
}
plot(finalsize$R0,finalsize$finalsize, type ="l")

###############################################################################################
#                                                                                             #
#                         SEIR -     Deterministic without demography                          #
#                                                                                             #
###############################################################################################

#Ordinary differential Equations
seir <- function(Time, State, Pars){
  with(as.list(c(State,Pars)),{
    dx <- -1 * beta * X * Y / sum(State)
    dw <- beta * X * Y / sum(State) - sigma * W
    dy <- sigma * W  - gamma * Y
    dz <- gamma * Y
    return(list(c(dx,dw,dy,dz)))
  })
}

#Parameters
pars <- c(beta = 0.6, #beta = transmission constant
          gamma = 0.2, #gamma = recovery constant
          sigma = 0.5 #sigma = transient latent to infectious constant
)
#initial values
N0 = 30000 #population size
W0 = 0 #initial number latent
Y0 = 1 #initial number infected 
Z0 = 0# initial number recovered, 
X0 = N0 -  W0- Y0 - Z0 
init <- c(X = X0,W = W0, Y = Y0,Z = Z0)


#Data storage time
dt = 0.5#timestep for storing data
sim_duration  = 150 #length of the simulation
times <- seq(0, sim_duration, by = dt)


#Solve the ordinary differential equations
ode.out <- ode(init, times, seir, pars) 

#plot X, W, Y, and Z against time
plot(x = ode.out[,1], y = ode.out[,2], type = "l", xlab = "time", ylab = "animals", col = "purple",lwd =2, ylim = c(0,N0))
lines(x = ode.out[,1], y = ode.out[,3], type = "l", xlab = "time", ylab = "animals", col = "blue",lwd =2)
lines(x = ode.out[,1], y = ode.out[,4], type = "l", xlab = "time", ylab = "animals", col = "red",lwd =2)
lines(x = ode.out[,1], y = ode.out[,5], type = "l", xlab = "time", ylab = "animals", col = "palegreen1",lwd =2)
legend(x = 100, y = 22000, c("X","W","Y","Z"),lwd =2, col = c("purple", "blue","red", "palegreen1"))


#plot  X against Y
plot(x = ode.out[,2], y = ode.out[,3], type = "l", xlab ="X", ylab ="Y")
s <- seq(length(ode.out[,2])-1)
arrows(x0 = ode.out[s,2], y0 = ode.out[s,3], x1 = ode.out[s + 1,2], y1 = ode.out[s + 1 ,3],length = .075)



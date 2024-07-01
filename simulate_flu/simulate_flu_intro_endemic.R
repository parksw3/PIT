library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../R/model_gog.R")

R0 <- 1.8
s <- 0.07
n <- 182
nscale <- 10
gamma <- 1/5
mu <- 1/30/365

x <- 166
y <- 182

sigmamat <- pmax(1-outer(1:n, 1:n, function(x, y) abs(x-y)*s/nscale),0)

yini <-  c(rep(1, n), rep(0, n))
names(yini)[1:n] <- paste0("S", 1:n)
names(yini)[(n+1):(2*n)] <- paste0("I", 1:n)		

# mu - beta * S * I - mu * S = 0
# mu - (mu + gamma) * I - mu/R0 = 0
# mu * (1-1/R0)/(mu+gamma) = I

# mu - beta * S * sigma * I - mu * S = 0
# mu -  mu * (R0-1) * S * sigma - mu * S = 0
# mu =(mu + mu * (R0-1) * sigma) * S
# S = mu/(mu + mu * (R0-1) * sigma)

yini_x <- yini_y <- yini
yini_x[1:n] <- mu/(mu+mu*(R0-1)*sigmamat[,x])
yini_x["I166"] <- (mu - mu/R0)/(gamma)

simulate_flu_intro_endemic_x <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																						 gamma=1/5, mu=1/30/365,
																						 s=0.07, n=n, nscale=10,
																						 yini=yini_x, tvec=(-200*365-180):(30*365+185))

yini_xy <- tail(head(unlist(simulate_flu_intro_endemic_x[simulate_flu_intro_endemic_x$time==3*365-180,]),-2),-1)
yini_xy["I182"] <- 1e-6

simulate_flu_intro_endemic_xy  <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																							 gamma=1/5, mu=1/30/365,
																							 s=0.07, n=n, nscale=10,
																							 yini=yini_xy, tvec=(2*365-180):(30*365+185))

yini_y[1:n] <- mu/(mu+mu*(R0-1)*sigmamat[,y])
yini_y["I182"] <- (mu - mu/R0)/(gamma)

simulate_flu_intro_endemic_y <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																						 gamma=1/5, mu=1/30/365,
																						 s=0.07, n=n, nscale=10,
																						 yini=yini_y, tvec=(-200*365-180):(30*365+185))

yini_yx <- tail(head(unlist(simulate_flu_intro_endemic_y[simulate_flu_intro_endemic_y$time==3*365-180,]),-2),-1)
yini_yx["I166"] <- 1e-6

simulate_flu_intro_endemic_yx  <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																							 gamma=1/5, mu=1/30/365,
																							 s=0.07, n=n, nscale=10,
																							 yini=yini_yx, tvec=(2*365-180):(30*365+185))

simulate_flu_intro_endemic_x <- simulate_flu_intro_endemic_x %>%
	filter(time %in% (1*365-180):(30*365+185))

simulate_flu_intro_endemic_y <- simulate_flu_intro_endemic_y %>%
	filter(time %in% (1*365-180):(30*365+185))

save("simulate_flu_intro_endemic_x", "simulate_flu_intro_endemic_y",
		 "simulate_flu_intro_endemic_xy", "simulate_flu_intro_endemic_yx",
		 file="simulate_flu_intro_endemic.rda")

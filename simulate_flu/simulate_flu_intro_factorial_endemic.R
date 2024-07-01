library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../R/model_gog.R")

n <- 2
mu <- 1/50/365
gamma <- 1/5
nscale <- 1

Reff <- 1.05

R0vec <- c(1.5, 3, 6, 12)
R0ratiovec <- exp(seq(log(0.5), log(2), length.out=101))
svec <- seq(0, 1, length.out=101)

paramdata <- expand.grid(R0vec, R0ratiovec, svec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
	print(i)
	pp <- paramdata[i,]
	R0 <- c(pp[[1]], pp[[1]]*pp[[2]])
	s <- pp[[3]]
	
	S_11 <- 1/R0[1]
	S_21 <- 1/R0[1] + (1-1/R0[1]) * s
	
	S_22 <- 1/R0[2]
	S_12 <- 1/R0[2] + (1-1/R0[2]) * s
	
	R_1 <- R0[1]
	R_2 <- R0[2]
	
	niche <- sqrt(S_11*S_22/(S_12*S_21))
	
	fitness <- R_2/R_1*sqrt(S_22*S_21/(S_11*S_12))
	
	sigmamat <- pmax(1-outer(1:n, 1:n, function(x, y) abs(x-y)*s/nscale),0)
	
	yini <-  c(rep(1, n), rep(0, n))
	names(yini)[1:n] <- paste0("S", 1:n)
	names(yini)[(n+1):(2*n)] <- paste0("I", 1:n)		
	
	yini_xy <- yini_yx <- yini
	yini_xy[1] <- 1/R0[1]
	yini_xy[2] <- 1/R0[2]*Reff
	yini_xy["I1"] <- (mu - mu/R0[1])/(gamma+mu)
	yini_xy["I2"] <- 1e-7
	
	ss_xy <- simulate_gog(R0=R0, Reff=1, phi=0,
												gamma=gamma, mu=mu,
												s=s, n=n, nscale=nscale,
												yini=yini_xy, tvec=0:(365*5))
	
	reslist[[i]] <- data.frame(
		R0=pp[[1]],
		R0ratio=pp[[2]],
		Reffinv=yini_xy[[2]]*R0[2],
		cross=1-s,
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		xy_trough_x=min(ss_xy$I1),
		xy_trough_y=min(ss_xy$I2)
	)
}

simulate_flu_intro_factorial_endemic <- reslist %>%
	bind_rows

save("simulate_flu_intro_factorial_endemic", file="simulate_flu_intro_factorial_endemic.rda")

library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../R/model_white.R")

tmax <- 100

immune_dur_vec <- exp(seq(log(1), log(20), length.out=21))
alphavec <- seq(1, 1.8, length.out=21)
eta <- 0.4126
theta <- 0
sigma_ho <- 0.3569
sigma_he <- 0.8426
alpha_base <- 0.9159

paramdata <- expand.grid(immune_dur_vec, alphavec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
	print(i)
	pp <- paramdata[i,]
	immune_dur <- pp[[1]]
	alpha <- pp[[2]]
	delta <- 1/immune_dur/52
	
	out1 <- simulate_white_resident1(alpha=1/alpha,
																	 eta=eta,
																	 delta=delta,
																	 b_1=99.51/52*alpha_base*alpha,
																	 theta=0,
																	 tmin=0,
																	 tmax=100)
	
	yini12 <- unlist(tail(out1,1)[2:13])
	yini12["Ip2"] <- 1e-100
	
	out12 <- simulate_white(alpha=1/alpha,
													eta=eta,
													delta=delta,
													b_1=99.51/52*alpha_base*alpha,
													theta=0,
													tmin=0,
													tmax=3,
													yini=yini12)
	
	out2 <- simulate_white_resident2(alpha=1/alpha,
																	 eta=eta,
																	 delta=delta,
																	 b_1=99.51/52*alpha_base*alpha,
																	 theta=0,
																	 tmin=0,
																	 tmax=100)
	
	yini21 <- unlist(tail(out2,1)[2:13])
	yini21["Ip1"] <- 1e-100
	
	out21 <- simulate_white(alpha=1/alpha,
													eta=eta,
													delta=delta,
													b_1=99.51/52*alpha_base*alpha,
													theta=0,
													tmin=0,
													tmax=3,
													yini=yini21)
	
	out12_filter <- head(tail(out12, -53), -1)
	out21_filter <- head(tail(out21, -53), -1)
	
	out1_filter <- out1[98 < out1$time/52 & out1$time/52 < 100,]
	out2_filter <- out2[98 < out2$time/52 & out2$time/52 < 100,]
	
	S_21 <- out1_filter$S + sigma_he * out1_filter$R1
	
	S_12 <- out2_filter$S + sigma_he * out2_filter$R2
	
	R_1 <- out21_filter$R_1 * 
		(out21_filter$Ip1 + eta * (out21_filter$Is_he1 + out21_filter$Is_ho1 + out21_filter$It1))/
		(out21_filter$Ip1 + out21_filter$Is_he1 + out21_filter$Is_ho1 + out21_filter$It1)
	
	R_2 <- out12_filter$R_2 * 
		(out12_filter$Ip2 + eta * (out12_filter$Is_he2 + out12_filter$Is_ho2 + out12_filter$It2))/
		(out12_filter$Ip2 + out12_filter$Is_he2 + out12_filter$Is_ho2 + out12_filter$It2)
	
	niche <- mean(sqrt(1/(R_1 * R_2 * S_12*S_21)))
	
	fitness <- mean(sqrt(R_2/R_1*S_21/(S_12)))
	
	yini21 <- unlist(tail(out2,1)[2:13])
	yini21["Ip1"] <- 1e-7
	
	## always want to simulate a more fit variant invading a less fit variant
	sim21 <- simulate_white(alpha=1/alpha,
													eta=eta,
													delta=delta,
													b_1=99.51/52*alpha_base*alpha,
													theta=0,
													tmin=0,
													tmax=10,
													yini=yini21)
	
	S_AB1 <- sim21$S + sigma_ho * sim21$R1 + sigma_he * sim21$R2 + sigma_ho * sigma_he * sim21$R12
	S_AB2 <- sim21$S + sigma_ho * sim21$R2 + sigma_he * sim21$R1 + sigma_ho * sigma_he * sim21$R12
	
	Reff_AB1 <- sim21$R_1 * (sim21$Ip1 + eta * (sim21$Is_he1 + sim21$Is_ho1 + sim21$It1))/(sim21$Ip1 + (sim21$Is_he1 + sim21$Is_ho1 + sim21$It1))
	Reff_AB2 <- sim21$R_2 * (sim21$Ip2 + eta * (sim21$Is_he2 + sim21$Is_ho2 + sim21$It2))/(sim21$Ip2 + (sim21$Is_he2 + sim21$Is_ho2 + sim21$It2))
	Reff_AB1[1:5] <- NA
	
	RR <- Reff_AB2 * S_AB2
	
	ret <- which(head(RR,-1) < 0.99 & tail(RR,-1) > 0.99)[1]
	if (is.na(ret)) {
		ret <- nrow(sim21)
		meetmax <- TRUE
	} else {
		meetmax <- FALSE
	}
	
	returntime <- ret - which(RR < 0.99)[1]
	overcomp <- 1-min(RR[1:ret])
	
	reslist[[i]] <- data.frame(
		immune_dur=immune_dur,
		alpha=alpha,
		A_trough=min(sim21$prevalence1),
		B_trough=min(sim21$prevalence2),
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		returntime=returntime,
		overcomp=overcomp,
		meetmax=meetmax
	)
}

simulate_rsv_factorial <- reslist %>%
	bind_rows

save("simulate_rsv_factorial", file="simulate_rsv_factorial.rda")

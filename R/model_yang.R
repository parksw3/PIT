model_yang <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		beta1 <- R01 * (1 + theta * cos(2 * pi * t/52)) * (mu + 1/D1)
		beta2 <- R02 * (1 + theta * cos(2 * pi * t/52)) * (mu + 1/D2)
		
		dS1 <- - beta1 * I1 * S1 - c12 * beta2 * I2 * S2 + mu * (1-S1) + (1-S1-I1)/L1
		dI1 <- beta1 * I1 * S1 - I1/D1 - mu * I1
		
		dS2 <- - c21 * beta1 * I1 * S1 - beta2 * I2 * S2 + mu * (1-S1) + (1-S2-I2)/L2
		dI2 <- beta2 * I2 * S2 - I2/D2 - mu * I2
		
		list(c(dS1, dI1, dS2, dI2),
				 R_1=beta1/(mu + 1/D1),
				 R_2=beta2/(mu + 1/D2))
	})
}

simulate_yang <- function(R01=1.44,
													R02=1.60,
													theta=0.1,
													D1=2.64/7,
													D2=3.03/7,
													L1=3.12*52,
													L2=2.28*52,
													mu=1/80/52,
													c12=0.17,
													c21=0.4,
													yini,
													tmin=0,
													tmax=100) {
	parms <- c(R01=R01,
						 R02=R02,
						 theta=theta,
						 D1=D1,
						 D2=D2,
						 L1=L1,
						 L2=L2,
						 mu=mu,
						 c12=c12,
						 c21=c21)
	
	if (missing(yini)) {
		yini <- c(S1=1-1e-6, I1=1e-6, S2=1-1e-6, I2=1e-6)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(rk(yini, times, model_yang, parms))
	
	out	
}

simulate_yang_resident1 <- function(R01=1.44,
																		R02=1.60,
																		theta=0.1,
																		D1=2.64/7,
																		D2=3.03/7,
																		L1=3.12*52,
																		L2=2.28*52,
																		mu=1/80/52,
																		c12=0.17,
																		c21=0.4,
																		tmin=0,
																		tmax=100) {
	yini <- c(S1=1-1e-6, I1=1e-6, S2=1, I2=0) 
	
	out <- simulate_yang(R01=R01, R02=R02, theta=theta, D1=D1, D2=D2, L1=L1, L2=L2, mu=mu,
											 c12=c12, c21=c21, yini=yini, tmin=tmin, tmax=tmax)
	
	out
}

simulate_yang_resident2 <- function(R01=1.44,
																		R02=1.60,
																		theta=0.1,
																		D1=2.64/7,
																		D2=3.03/7,
																		L1=3.12*52,
																		L2=2.28*52,
																		mu=1/80/52,
																		c12=0.17,
																		c21=0.4,
																		tmin=0,
																		tmax=100) {
	yini <- c(S1=1, I1=0, S2=1-1e-6, I2=1e-6) 
	
	out <- simulate_yang(R01=R01, R02=R02, theta=theta, D1=D1, D2=D2, L1=L1, L2=L2, mu=mu,
											 c12=c12, c21=c21, yini=yini, tmin=tmin, tmax=tmax)
	
	out
}

simulate_yang_niche <- function(R01=1.44,
																 R02=1.60,
																 theta=0.1,
																 D1=2.64/7,
																 D2=3.03/7,
																 L1=3.12*52,
																 L2=2.28*52,
																 mu=1/80/52,
																 c12=0.17,
																 c21=0.4,
																 tmin=0,
																 tmax=200,
																 invmin=180,
																 invmax=200,
																 system="Flu H1N1, H3N2") { 
	out1 <- simulate_yang_resident1(
		R01=R01, R02=R02, theta=theta, D1=D1, D2=D2, L1=L1, L2=L2, mu=mu,
		c12=c12, c21=c21, tmin=tmin, tmax=tmax)
	
	out2 <- simulate_yang_resident2(
		R01=R01, R02=R02, theta=theta, D1=D1, D2=D2, L1=L1, L2=L2, mu=mu,
		c12=c12, c21=c21, tmin=tmin, tmax=tmax)
	
	out1_filter <- out1[invmin < out1$time/52 & out1$time/52 < invmax,]
	out2_filter <- out2[invmin < out2$time/52 & out2$time/52 < invmax,]
	
	S_11 <- out1_filter$S1
	S_21 <- out1_filter$S2
	
	S_22 <- out2_filter$S2
	S_12 <- out2_filter$S1
	
	R_1 <- out1_filter$R_1
	R_2 <- out2_filter$R_2
	
	niche <- mean(sqrt(S_11*S_22/(S_12*S_21)))
	
	fitness <- mean(R_2/R_1*sqrt(S_22*S_21/(S_11*S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		system=system,
		meanR0=mean(sqrt(R_1*R_2)),
		nulldiff=pmax(mean(sqrt(R_1/R_2)), mean(sqrt(R_2/R_1)))
	)
}

# parameters from https://royalsocietypublishing.org/doi/pdf/10.1098/rspb.2005.3454
model_rohani <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		beta1 <- b01 * (1 + theta * cos(2*pi*t/52))
		beta2 <- b02 * (1 + theta * cos(2*pi*t/52))
		
		lambda1 <- beta1 * (I1T + I1R)
		lambda2 <- beta2 * (I2T + I2R)
		
		dS12 <- mu - (lambda1 + lambda2 + mu) * S12
		dS1R <- delta2 * C2T - (lambda1 + mu) * S1R
		dS2R <- delta1 * C1T - (lambda2 + mu) * S2R
		dE1T <- lambda1 * S12 - (mu + sigma1) * E1T
		dE1R <- lambda1 * S1R - (mu + sigma1) * E1R
		dE2T <- lambda2 * S12 - (mu + sigma2) * E2T
		dE2R <- lambda2 * S2R - (mu + sigma2) * E2R
		dI1T <- sigma1 * E1T - (mu + gamma1) * I1T
		dI1R <- sigma1 * E1R - (mu + gamma1) * I1R
		dI2T <- sigma2 * E2T - (mu + gamma2) * I2T
		dI2R <- sigma2 * E2R - (mu + gamma2) * I2R
		dC1T <- gamma1 * I1T - (mu + delta1) * C1T
		dC2T <- gamma2 * I2T - (mu + delta2) * C2T
		
		list(c(dS12, 
					 dS1R, dS2R, dE1T, dE1R, dE2T, dE2R,
					 dI1T, dI1R, dI2T, dI2R, dC1T, dC2T),
				 R_1=beta1/(gamma1+mu) * sigma1/(sigma1+mu),
				 R_2=beta2/(gamma2+mu) * sigma2/(sigma2+mu))
	})
}

simulate_rohani <- function(b01=1250/52,
														sigma1=7/8,
														gamma1=7/5,
														b02=625/52,
														sigma2=7/8,
														gamma2=7/10,
														delta1=7/7, 
														delta2=7/14,
														theta=0.1,
														mu=1/50/52,
														yini,
														tmin=0,
														tmax=100) {
	parms <- c(b01=b01,
						 sigma1=sigma1,
						 gamma1=gamma1,
						 b02=b02,
						 sigma2=sigma2,
						 gamma2=gamma2,
						 delta1=delta1,
						 delta2=delta2,
						 theta=theta,
						 mu=mu)
	
	if (missing(yini)) {
		yini <- c(S12=1-2e-6, 
							S1R=0, S2R=0, E1T=0, E1R=0, E2T=0, E2R=0,
							I1T=1e-6, I1R=0, I2T=1e-6, I2R=0, C1T=0, C2T=0)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(rk(yini, times, model_rohani, parms))
	
	out
}

simulate_rohani_resident1 <- function(b01=1250/52,
																			sigma1=7/8,
																			gamma1=7/5,
																			b02=625/52,
																			sigma2=7/8,
																			gamma2=7/10,
																			delta1=7/7, 
																			delta2=7/14,
																			theta=0.1,
																			mu=1/50/52,
																			yini,
																			tmin=0,
																			tmax=100) {
	yini <- c(S12=1-1e-6, 
						S1R=0, S2R=0, E1T=0, E1R=0, E2T=0, E2R=0,
						I1T=1e-6, I1R=0, I2T=0, I2R=0, C1T=0, C2T=0)
	
	out <- simulate_rohani(b01=b01,
												 sigma1=sigma1,
												 gamma1=gamma1,
												 b02=b02,
												 sigma2=sigma2,
												 gamma2=gamma2,
												 delta1=delta1,
												 delta2=delta2,
												 theta=theta,
												 mu=mu,
												 yini=yini,
												 tmin=tmin,
												 tmax=tmax)
	
	out
}

simulate_rohani_resident2 <- function(b01=1250/52,
																			sigma1=7/8,
																			gamma1=7/5,
																			b02=625/52,
																			sigma2=7/8,
																			gamma2=7/10,
																			delta1=7/7, 
																			delta2=7/14,
																			theta=0.1,
																			mu=1/50/52,
																			yini,
																			tmin=0,
																			tmax=100) {
	yini <- c(S12=1-1e-6, 
						S1R=0, S2R=0, E1T=0, E1R=0, E2T=0, E2R=0,
						I1T=0, I1R=0, I2T=1e-6, I2R=0, C1T=0, C2T=0)
	
	out <- simulate_rohani(b01=b01,
												 sigma1=sigma1,
												 gamma1=gamma1,
												 b02=b02,
												 sigma2=sigma2,
												 gamma2=gamma2,
												 delta1=delta1,
												 delta2=delta2,
												 theta=theta,
												 mu=mu,
												 yini=yini,
												 tmin=tmin,
												 tmax=tmax)
	
	out
}

simulate_rohani_niche <- function(b01=1250/52,
																	sigma1=7/8,
																	gamma1=7/5,
																	b02=625/52,
																	sigma2=7/8,
																	gamma2=7/10,
																	delta1=7/7, 
																	delta2=7/14,
																	theta=0.2,
																	mu=1/50/52,
																	tmin=0,
																	tmax=100,
																	invmin=90,
																	invmax=100,
																	system="Measles, whooping cough") { 
	out1 <- simulate_rohani_resident1(b01=b01,
																		sigma1=sigma1,
																		gamma1=gamma1,
																		b02=b02,
																		sigma2=sigma2,
																		gamma2=gamma2,
																		delta1=delta1,
																		delta2=delta2,
																		theta=theta,
																		mu=mu,
																		yini=yini,
																		tmin=tmin,
																		tmax=tmax)
	
	out2<- simulate_rohani_resident2(b01=b01,
																		sigma1=sigma1,
																		gamma1=gamma1,
																		b02=b02,
																		sigma2=sigma2,
																		gamma2=gamma2,
																		delta1=delta1,
																		delta2=delta2,
																		theta=theta,
																		mu=mu,
																		yini=yini,
																		tmin=tmin,
																		tmax=tmax)
	
	out1_filter <- out1[invmin < out1$time/52 & out1$time/52 < invmax,]
	out2_filter <- out2[invmin < out2$time/52 & out2$time/52 < invmax,]
	
	S_11 <- out1_filter$S12
	S_21 <- out1_filter$S12 + out1_filter$S2R
	
	S_22 <- out2_filter$S12
	S_12 <- out2_filter$S12 + out2_filter$S1R
	
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

model_white <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		beta_1 <- b_1 * (1 + theta * cos(2 * pi * (t-phi)/52))
		beta_2 <- alpha * beta_1
		
		lambda_1 <- beta_1 * (Ip1 + eta * (Is_he1 + Is_ho1 + It1))
		lambda_2 <- beta_2 * (Ip2 + eta * (Is_he2 + Is_ho2 + It2))
		
		dS <- mu - (lambda_1 + lambda_2 + mu) * S + delta * (R1 + R2)
		dIp1 <- lambda_1 * S - (gamma + mu) * Ip1
		dIp2 <- lambda_2 * S - (gamma + mu) * Ip2
		dR1 <- gamma * Ip1 - (sigma_ho * lambda_1 + sigma_he * lambda_2 + delta  + mu) * R1 + gamma * Is_ho1 + delta * R12
		dR2 <- gamma * Ip2 - (sigma_ho * lambda_2 + sigma_he * lambda_1 + delta  + mu) * R2 + gamma * Is_ho2 + delta * R12
		dIs_he1 <- sigma_he * lambda_1 * R2 - (gamma + mu) * Is_he1
		dIs_he2 <- sigma_he * lambda_2 * R1 - (gamma + mu) * Is_he2
		dIs_ho1 <- sigma_ho * lambda_1 * R1 - (gamma + mu) * Is_ho1
		dIs_ho2 <- sigma_ho * lambda_2 * R2 - (gamma + mu) * Is_ho2
		dR12 <- gamma * (Is_he1 + Is_he2 + It1 + It2) - (sigma_he * sigma_ho * (lambda_1 + lambda_2) + 2 * delta + mu) * R12
		dIt1 <- sigma_he * sigma_ho * lambda_1 * R12 - (gamma + mu) * It1
		dIt2 <- sigma_he * sigma_ho * lambda_2 * R12 - (gamma + mu) * It2
		
		prevalence1 <- Ip1 + Is_he1 + Is_ho1 + It1
		prevalence2 <- Ip2 + Is_he2 + Is_ho2 + It2
		
		list(c(dS, dIp1, dIp2, dR1, dR2, dIs_he1, dIs_he2, dIs_ho1, dIs_ho2, dR12, dIt1, dIt2),
				 R_1=beta_1/(gamma+mu), R_2=beta_2/(gamma+mu),
				 prevalence1=prevalence1, prevalence2=prevalence2)
	})
}

simulate_white <- function(alpha=0.9159,
													 eta=0.4126,
													 sigma_ho=0.3569,
													 sigma_he=0.8426,
													 gamma=40.56/52,
													 delta=0.51/52,
													 b_1=99.51/52,
													 mu=0.012/52,
													 phi=0.97*52,
													 theta=0.347,
													 yini,
													 tmin=0,
													 tmax=100) {
	parms <- c(alpha=alpha,
						 eta=eta,
						 sigma_ho=sigma_ho,
						 sigma_he=sigma_he,
						 gamma=gamma,
						 delta=delta,
						 b_1=b_1,
						 mu=mu,
						 phi=phi,
						 theta=theta)
	
	if (missing(yini)) {
		yini <- c(S=1-2e-6, Ip1=1e-6, Ip2=1e-6,
							R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(rk(yini, times, model_white, parms))
	
	out
}

simulate_white_resident1 <- function(alpha=0.9159,
																		 eta=0.4126,
																		 sigma_ho=0.3569,
																		 sigma_he=0.8426,
																		 gamma=40.56/52,
																		 delta=0.51/52,
																		 b_1=99.51/52,
																		 mu=0.012/52,
																		 phi=0.97*52,
																		 theta=0.347,
																		 tmin=0,
																		 tmax=50) {
	yini <- c(S=1-1e-6, Ip1=1e-6, Ip2=0,
						R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)
	
	out <- simulate_white(alpha=alpha,
												eta=eta,
												sigma_ho=sigma_ho,
												sigma_he=sigma_he,
												gamma=gamma,
												delta=delta,
												b_1=b_1,
												mu=mu,
												phi=phi,
												theta=theta,
												yini=yini,
												tmin=tmin,
												tmax=tmax)
	
	out
}

simulate_white_resident2 <- function(alpha=0.9159,
																		 eta=0.4126,
																		 sigma_ho=0.3569,
																		 sigma_he=0.8426,
																		 gamma=40.56/52,
																		 delta=0.51/52,
																		 b_1=99.51/52,
																		 mu=0.012/52,
																		 phi=0.97*52,
																		 theta=0.347,
																		 tmin=0,
																		 tmax=50) {
	R0 <- b_1/(gamma+mu) * alpha
	
	yini <- c(S=1-1e-6, Ip1=0, Ip2=1e-6,
						R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)
	
	out <- simulate_white(alpha=alpha,
												eta=eta,
												sigma_ho=sigma_ho,
												sigma_he=sigma_he,
												gamma=gamma,
												delta=delta,
												b_1=b_1,
												mu=mu,
												phi=phi,
												theta=theta,
												yini=yini,
												tmin=tmin,
												tmax=tmax)
	
	out
}

simulate_white_niche <- function(alpha=0.9159,
																		eta=0.4126,
																		sigma_ho=0.3569,
																		sigma_he=0.8426,
																		gamma=40.56/52,
																		delta=0.51/52,
																		b_1=99.51/52,
																		mu=0.012/52,
																		phi=0.97*52,
																		theta=0.347,
																 tmin=0,
																 tmax=100,
																 invmin=90,
																 invmax=100) { 
	out1 <- simulate_white_resident1(alpha=alpha,
																	 eta=eta,
																	 sigma_ho=sigma_ho,
																	 sigma_he=sigma_he,
																	 gamma=gamma,
																	 delta=delta,
																	 b_1=b_1,
																	 mu=mu,
																	 phi=phi,
																	 theta=theta,
																	 tmin=tmin,
																	 tmax=tmax)
	
	yini12 <- unlist(tail(out1,1)[2:13])
	yini12["Ip2"] <- 1e-100
	
	out12 <- simulate_white(alpha=alpha,
													eta=eta,
													sigma_ho=sigma_ho,
													sigma_he=sigma_he,
													gamma=gamma,
													delta=delta,
													b_1=b_1,
													mu=mu,
													phi=phi,
													theta=theta,
													tmin=0,
													tmax=11,
													yini=yini12)
	
	out2 <- simulate_white_resident2(alpha=alpha,
																	 eta=eta,
																	 sigma_ho=sigma_ho,
																	 sigma_he=sigma_he,
																	 gamma=gamma,
																	 delta=delta,
																	 b_1=b_1,
																	 mu=mu,
																	 phi=phi,
																	 theta=theta,
																	 tmin=tmin,
																	 tmax=tmax)
	
	yini21 <- unlist(tail(out2,1)[2:13])
	yini21["Ip1"] <- 1e-100
	
	out21 <- simulate_white(alpha=alpha,
													eta=eta,
													sigma_ho=sigma_ho,
													sigma_he=sigma_he,
													gamma=gamma,
													delta=delta,
													b_1=b_1,
													mu=mu,
													phi=phi,
													theta=theta,
													tmin=0,
													tmax=11,
													yini=yini21)
	
	out12_filter <- head(tail(out12, -53), -1)
	out21_filter <- head(tail(out21, -53), -1)
	
	out1_filter <- out1[invmin < out1$time/52 & out1$time/52 < invmax,]
	out2_filter <- out2[invmin < out2$time/52 & out2$time/52 < invmax,]
	
	# S_21 means susceptibility of 2 when 1 is at resident
	# so we want out_12 (2 invading when 1 is resident)
	# bad naming... but OK
	S_21 <- out12_filter$S + sigma_he * out12_filter$R1
	
	S_12 <- out21_filter$S + sigma_he * out21_filter$R2

	R_1 <- out21_filter$R_1 * (out21_filter$Ip1 + eta * out21_filter$Is_he1)/(out21_filter$Ip1 + out21_filter$Is_he1)
		
	R_2 <- out21_filter$R_2 * (out12_filter$Ip2 + eta * out12_filter$Is_he2)/(out12_filter$Ip2 + out12_filter$Is_he2)
	
	niche <- mean(sqrt(1/(R_1*R_2*S_12*S_21)))
	
	fitness <- mean(sqrt(R_2/R_1*S_21/(S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		system="RSV A, B",
		meanR0=mean(sqrt(R_1*R_2)),
		nulldiff=mean(sqrt(R_1/R_2))
	)
}

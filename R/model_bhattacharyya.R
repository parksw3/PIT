model_bhattacharyya <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		beta1 <- b01 * (1 + theta1 * sin(2*pi*(t-phi1)/52)) * q1
		beta2 <- b02 * (1 + theta2 * sin(2*pi*(t-phi2)/52)) * q2
		
		prevalence1 <- I1 + J1
		prevalence2 <- I2 + J2
		
		lambda1 <- beta1 * prevalence1
		lambda2 <- beta2 * prevalence2
		
		dS <- mu - mu * S - (lambda1 + lambda2) * S + rho1 * R1 + rho2 * R2
		dI1 <- lambda1 * S - (gamma1 + mu) * I1
		dI2 <- lambda2 * S - (gamma2 + mu) * I2
		dR1 <- gamma1 * I1 - lambda2 * epsilon21 * R1 - rho1 * R1 + rho2 * R - mu * R1
		dR2 <- gamma2 * I2 - lambda1 * epsilon12 * R2 - rho2 * R2 + rho1 * R - mu * R2
		dJ1 <- lambda1 * epsilon12 * R2 - (gamma1 + mu) * J1
		dJ2 <- lambda2 * epsilon21 * R1 - (gamma2 + mu) * J2
		dR <- gamma1 * J1 + gamma2 * J2 - rho1 * R - rho2 * R - mu * R
		
		list(c(dS, dI1, dI2, dR1, dR2, dJ1, dJ2, dR),
				 prevalence1=prevalence1,
				 prevalence2=prevalence2,
				 R_1=beta1/(gamma1+mu),
				 R_2=beta2/(gamma2+mu))
	})
}

simulate_bhattacharyya <- function(b01=3.4,
																	 b02=2.9,
																	 theta1=0.4,
																	 theta2=0.33,
																	 phi1=-17.5,
																	 phi2=-29.5,
																	 epsilon12=0.9,
																	 epsilon21=0.54,
																	 gamma1=7/10,
																	 gamma2=7/10,
																	 rho1=1/52,
																	 rho2=1/52,
																	 q1=0.5,
																	 q2=0.5,
																	 mu=1/70/52,
																	 yini,
																	 tmin=0,
																	 tmax=30) {
	parms <- c(b01=b01, b02=b02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2, epsilon12=epsilon12, epsilon21=epsilon21,
						 gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu)
	
	if (missing(yini)) {
		yini <- c(S=1-2e-6, I1=1e-6, I2=1e-6, R1=0, R2=0, J1=0, J2=0, R=0)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(rk(yini, times, model_bhattacharyya, parms))
	
	out
}

simulate_bhattacharyya_resident1 <- function(b01=3.4,
																						 b02=2.9,
																						 theta1=0.4,
																						 theta2=0.33,
																						 phi1=-17.5,
																						 phi2=-29.5,
																						 epsilon12=0.9,
																						 epsilon21=0.54,
																						 gamma1=7/10,
																						 gamma2=7/10,
																						 rho1=1/52,
																						 rho2=1/52,
																						 q1=0.5,
																						 q2=0.5,
																						 mu=1/70/52,
																						 tmin=0,
																						 tmax=30) {
	yini <- c(S=1-1e-6, I1=1e-6, I2=0, R1=0, R2=0, J1=0, J2=0, R=0)
	
	out <- simulate_bhattacharyya(
		b01=b01, b02=b02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2, epsilon12=epsilon12, epsilon21=epsilon21,
		gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu,
		yini=yini,
		tmin=tmin,
		tmax=tmax)
	
	out
}

simulate_bhattacharyya_resident2 <- function(b01=3.4,
																						 b02=2.9,
																						 theta1=0.4,
																						 theta2=0.33,
																						 phi1=-17.5,
																						 phi2=-29.5,
																						 epsilon12=0.9,
																						 epsilon21=0.54,
																						 gamma1=7/10,
																						 gamma2=7/10,
																						 rho1=1/52,
																						 rho2=1/52,
																						 q1=0.5,
																						 q2=0.5,
																						 mu=1/70/52,
																						 tmin=0,
																						 tmax=30) {
	yini <- c(S=1-1e-6, I1=0, I2=1e-6, R1=0, R2=0, J1=0, J2=0, R=0)
	
	out <- simulate_bhattacharyya(
		b01=b01, b02=b02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2, epsilon12=epsilon12, epsilon21=epsilon21,
		gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu,
		yini=yini,
		tmin=tmin,
		tmax=tmax)
	
	out
}

simulate_bhattacharyya_niche <- function(b01=3.4,
																				 b02=2.9,
																				 theta1=0.4,
																				 theta2=0.33,
																				 phi1=-17.5,
																				 phi2=-29.5,
																				 epsilon12=0.9,
																				 epsilon21=0.54,
																				 gamma1=7/10,
																				 gamma2=7/10,
																				 rho1=1/52,
																				 rho2=1/52,
																				 q1=0.5,
																				 q2=0.5,
																				 mu=1/70/52,
																				 tmin=0,
																				 tmax=30,
																				 invmin=20,
																				 invmax=30,
																				 system="RSV, HPIV-1") { 
	out1 <- simulate_bhattacharyya_resident1(
		b01=b01, b02=b02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2, epsilon12=epsilon12, epsilon21=epsilon21,
		gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu,
		tmin=tmin,
		tmax=tmax)
	
	out2 <- simulate_bhattacharyya_resident2(
		b01=b01, b02=b02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2, epsilon12=epsilon12, epsilon21=epsilon21,
		gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu,
		tmin=tmin,
		tmax=tmax)
	
	out1_filter <- out1[invmin < out1$time/52 & out1$time/52 < invmax,]
	out2_filter <- out2[invmin < out2$time/52 & out2$time/52 < invmax,]
	
	S_11 <- out1_filter$S
	S_21 <- out1_filter$S + out1_filter$R1 * epsilon21
	
	S_22 <- out2_filter$S
	S_12 <- out2_filter$S + out2_filter$R2 * epsilon12
	
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

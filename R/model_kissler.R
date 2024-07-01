beta.val <- function(t, amplitude=1, baseline=1.5, phi.val=-4, gamma.val=1){
	gamma.val*(baseline+amplitude * cos(2*pi*(t-phi.val)/52))
}

# =============================================================================
# State update functions 
# =============================================================================

S1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1S2
})}
S1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*S1S2
})}
E1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1S2
})}
E1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*E1S2
})}
S1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1E2
})}
S1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*S1E2
})}
E1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1E2 
})}
E1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1E2 
})}
I1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1S2
})}
I1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*I1S2
})}
S1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1I2
})}
S1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*S1I2
})}
R1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1S2
})}
R1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*R1S2 
})}
I1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1E2
})}
I1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*I1E2
})}
E1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1I2
})}
E1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*E1I2
})}
S1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1R2
})}
S1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*S1R2
})}
R1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1E2
})}
R1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*R1E2
})}
I1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1I2 
})}
I1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1I2
})}
E1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1R2 
})}
E1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*E1R2
})}
R1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1I2
})}
R1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*R1I2
})}
I1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1R2
})}
I1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*I1R2
})}
R1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1R2
})}
R1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*R1R2
})}

model_kissler <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		dS1S2 <- mu -S1S2c1(t,y,parms) - S1S2c2(t,y,parms) + 
			R1S2c1(t,y,parms) + S1R2c2(t,y,parms) - mu * S1S2
		dE1S2 <- -E1S2c1(t,y,parms) - E1S2c2(t,y,parms) + 
			S1S2c1(t,y,parms) + E1R2c2(t,y,parms) - mu * E1S2
		dS1E2 <- -S1E2c1(t,y,parms) - S1E2c2(t,y,parms) + 
			R1E2c1(t,y,parms) + S1S2c2(t,y,parms) - mu * S1E2
		dE1E2 <- -E1E2c1(t,y,parms) - E1E2c2(t,y,parms) + 
			S1E2c1(t,y,parms) + E1S2c2(t,y,parms) - mu * E1E2
		dI1S2 <- -I1S2c1(t,y,parms) - I1S2c2(t,y,parms) + 
			E1S2c1(t,y,parms) + I1R2c2(t,y,parms) - mu * I1S2
		dS1I2 <- -S1I2c1(t,y,parms) - S1I2c2(t,y,parms) + 
			R1I2c1(t,y,parms) + S1E2c2(t,y,parms) - mu * S1I2
		dR1S2 <- -R1S2c1(t,y,parms) - R1S2c2(t,y,parms) + 
			I1S2c1(t,y,parms) + R1R2c2(t,y,parms) - mu * R1S2
		dI1E2 <- -I1E2c1(t,y,parms) - I1E2c2(t,y,parms) + 
			E1E2c1(t,y,parms) + I1S2c2(t,y,parms) - mu * I1E2
		dE1I2 <- -E1I2c1(t,y,parms) - E1I2c2(t,y,parms) + 
			S1I2c1(t,y,parms) + E1E2c2(t,y,parms) - mu * E1I2
		dS1R2 <- -S1R2c1(t,y,parms) - S1R2c2(t,y,parms) + 
			R1R2c1(t,y,parms) + S1I2c2(t,y,parms) - mu * S1R2
		dR1E2 <- -R1E2c1(t,y,parms) - R1E2c2(t,y,parms) + 
			I1E2c1(t,y,parms) + R1S2c2(t,y,parms) - mu * R1E2
		dI1I2 <- -I1I2c1(t,y,parms) - I1I2c2(t,y,parms) + 
			E1I2c1(t,y,parms) + I1E2c2(t,y,parms) - mu * I1I2
		dE1R2 <- -E1R2c1(t,y,parms) - E1R2c2(t,y,parms) + 
			S1R2c1(t,y,parms) + E1I2c2(t,y,parms) - mu * E1R2
		dR1I2 <- -R1I2c1(t,y,parms) - R1I2c2(t,y,parms) + 
			I1I2c1(t,y,parms) + R1E2c2(t,y,parms) - mu * R1I2
		dI1R2 <- -I1R2c1(t,y,parms) - I1R2c2(t,y,parms) + 
			E1R2c1(t,y,parms) + I1I2c2(t,y,parms) - mu * I1R2
		dR1R2 <- -R1R2c1(t,y,parms) - R1R2c2(t,y,parms) + 
			I1R2c1(t,y,parms) + R1I2c2(t,y,parms) - mu * R1R2
		
		prevalence1 <- I1S2 + I1E2 + I1I2 + I1R2
		prevalence2 <- S1I2 + E1I2 + I1I2 + R1I2
		
		return(list(c(dS1S2,dE1S2,dS1E2,dE1E2,dI1S2,dS1I2,dR1S2,dI1E2,dE1I2,dS1R2,dR1E2,dI1I2,dE1R2,dR1I2,dI1R2,dR1R2),
								prevalence1=prevalence1,
								prevalence2=prevalence2))
	})
}

simulate_kissler <- function(sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/45)
														 sigma2.val = 1/40, # Waning immunity rate, strain 2, weeks (def. 1/45)
														 nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
														 gamma.val = 1/(5/7),   # Rate of recovery, weeks (def. 1/1)
														 chi12.val = 0.7,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
														 chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
														 amplitude = 0.3,     # R0 seasonality amplitude (dev. 0.31)
														 baseline = 1.7,    # R0 seasonality baseline (def. 1.7)
														 phi.val =  -4,
														 mu=1/80/52,
														 yini,
														 tmin=0,
														 tmax=30) {
	
	parms <- c(sigma1.val = sigma1.val,
						 sigma2.val = sigma2.val,
						 nu.val = nu.val,
						 gamma.val = gamma.val,
						 chi12.val = chi12.val,
						 chi21.val = chi21.val,
						 amplitude = amplitude,
						 baseline = baseline,
						 phi.val = phi.val,
						 mu=mu)
	
	if (missing(yini)) {
		yini <- c(S1S2=1-2e-7,E1S2=1e-7,S1E2=1e-7,E1E2=0,I1S2=0,S1I2=0,R1S2=0,I1E2=0,E1I2=0,S1R2=0,R1E2=0,I1I2=0,E1R2=0,R1I2=0,I1R2=0,R1R2=0)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(rk(yini, times, model_kissler, parms))
	
	out
}

simulate_kissler_resident2 <- function(sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/45)
																			 sigma2.val = 1/40, # Waning immunity rate, strain 2, weeks (def. 1/45)
																			 nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
																			 gamma.val = 1/(5/7),   # Rate of recovery, weeks (def. 1/1)
																			 chi12.val = 0.7,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
																			 chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
																			 amplitude = 0.3,     # R0 seasonality amplitude (dev. 0.31)
																			 baseline = 1.7,    # R0 seasonality baseline (def. 1.7)
																			 phi.val =  -4,
																			 mu=1/80/52,
																			 tmin=0,
																			 tmax=100) {
	yini <- c(S1S2=1-1e-6,E1S2=0,S1E2=1e-6,E1E2=0,
						I1S2=0,S1I2=0,R1S2=0,I1E2=0,
						E1I2=0,S1R2=0,R1E2=0,I1I2=0,
						E1R2=0,R1I2=0,I1R2=0,R1R2=0)
	
	out <- simulate_kissler(sigma1.val = sigma1.val,
													sigma2.val = sigma2.val,
													nu.val = nu.val,
													gamma.val = gamma.val,
													chi12.val = chi12.val,
													chi21.val = chi21.val,
													amplitude = amplitude,
													baseline = baseline,
													phi.val = phi.val,
													mu=mu,
													yini=yini,
													tmin=tmin,
													tmax=tmax)
	
	out
}

simulate_kissler_resident1 <- function(sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/45)
																		 sigma2.val = 1/40, # Waning immunity rate, strain 2, weeks (def. 1/45)
																		 nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
																		 gamma.val = 1/(5/7),   # Rate of recovery, weeks (def. 1/1)
																		 chi12.val = 0.7,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
																		 chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
																		 amplitude = 0.3,     # R0 seasonality amplitude (dev. 0.31)
																		 baseline = 1.7,    # R0 seasonality baseline (def. 1.7)
																		 phi.val =  -4,
																		 mu=1/80/52,
																		 tmin=0,
																		 tmax=100) {
	yini <- c(S1S2=1-1e-6,E1S2=1e-6,S1E2=0,E1E2=0,
						I1S2=0,S1I2=0,R1S2=0,I1E2=0,
						E1I2=0,S1R2=0,R1E2=0,I1I2=0,
						E1R2=0,R1I2=0,I1R2=0,R1R2=0)
	
	out <- simulate_kissler(sigma1.val = sigma1.val,
		sigma2.val = sigma2.val,
		nu.val = nu.val,
		gamma.val = gamma.val,
		chi12.val = chi12.val,
		chi21.val = chi21.val,
		amplitude = amplitude,
		baseline = baseline,
		phi.val = phi.val,
		mu=mu,
		yini=yini,
		tmin=tmin,
		tmax=tmax)
	out
}

simulate_kissler_niche <- function(sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/45)
																			sigma2.val = 1/40, # Waning immunity rate, strain 2, weeks (def. 1/45)
																			nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
																			gamma.val = 1/(5/7),   # Rate of recovery, weeks (def. 1/1)
																			chi12.val = 0.7,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
																			chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
																			amplitude = 0.3,     # R0 seasonality amplitude (dev. 0.31)
																			baseline = 1.7,    # R0 seasonality baseline (def. 1.7)
																			phi.val =  -4,
																			mu=1/80/52,
																			tmin=0,
																			tmax=100,
																	 invmin=90,
																	 invmax=100) {
	out1 <- simulate_kissler_resident1(sigma1.val = sigma1.val,
														 sigma2.val = sigma2.val,
														 nu.val = nu.val,
														 gamma.val = gamma.val,
														 chi12.val = chi12.val,
														 chi21.val = chi21.val,
														 amplitude = amplitude,
														 baseline = baseline,
														 phi.val = phi.val,
														 mu=mu,
														 tmin=tmin,
														 tmax=tmax)
	
	out2 <- simulate_kissler_resident2(sigma1.val = sigma1.val,
																		 sigma2.val = sigma2.val,
																		 nu.val = nu.val,
																		 gamma.val = gamma.val,
																		 chi12.val = chi12.val,
																		 chi21.val = chi21.val,
																		 amplitude = amplitude,
																		 baseline = baseline,
																		 phi.val = phi.val,
																		 mu=mu,
																		 tmin=tmin,
																		 tmax=tmax)
	
	out1_filter <- out1[invmin < out1$time/52 & out1$time/52 < invmax,]
	out2_filter <- out2[invmin < out2$time/52 & out2$time/52 < invmax,]
	
	S_11 <- out1_filter$S1S2 
	S_21 <- out1_filter$S1S2 + (1-chi12.val) * (out1_filter$E1S2 + out1_filter$I1S2 + out1_filter$R1S2)
	
	S_22 <- out2_filter$S1S2 
	S_12 <- out2_filter$S1S2 + (1-chi21.val) * (out2_filter$S1E2 + out2_filter$S1I2 + out2_filter$S1R2)
	
	niche <- mean(sqrt(S_11*S_22/(S_12*S_21)))
	
	fitness <- mean(sqrt(S_22*S_21/(S_11*S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		system="HCoV HKU1, OC43",
		meanR0=baseline,
		nulldiff=1
	)
}

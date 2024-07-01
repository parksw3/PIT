model_cobey <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		prevalence <- I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I10
		lambda <- beta * prevalence
		
		gamma1 <- 1/(kappa + (nu - kappa) * exp(-0))
		gamma2 <- 1/(kappa + (nu - kappa) * exp(-epsilon))
		gamma3 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*2)))
		gamma4 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*3)))
		gamma5 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*4)))
		gamma6 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*5)))
		gamma7 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*6)))
		gamma8 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*7)))
		gamma9 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*8)))
		gamma10 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*9)))
		
		dS1 <- mu - lambda * S1 - mu * S1
		dI1 <- lambda * S1 - (gamma1 + mu) * I1
		dR1 <- gamma1 * I1 - (1-sigma) * lambda * R1 - mu * R1
		dI2 <- (1-sigma) * lambda * R1 - (gamma2 + mu) * I2
		dR2 <- gamma2 * I2 - (1-sigma) * lambda * R2 - mu * R2
		dI3 <- (1-sigma) * lambda * R2 - (gamma3 + mu) * I3
		dR3 <- gamma3 * I3 - (1-sigma) * lambda * R3 - mu * R3
		dI4 <- (1-sigma) * lambda * R3 - (gamma4 + mu) * I4
		dR4 <- gamma4 * I4 - (1-sigma) * lambda * R4 - mu * R4
		dI5 <- (1-sigma) * lambda * R4 - (gamma5 + mu) * I5
		dR5 <- gamma5 * I5 - (1-sigma) * lambda * R5 - mu * R5
		dI6 <- (1-sigma) * lambda * R5 - (gamma6 + mu) * I6
		dR6 <- gamma6 * I6 - (1-sigma) * lambda * R6 - mu * R6
		dI7 <- (1-sigma) * lambda * R6 - (gamma7 + mu) * I7
		dR7 <- gamma7 * I7 - (1-sigma) * lambda * R7 - mu * R7
		dI8 <- (1-sigma) * lambda * R7 - (gamma8 + mu) * I8
		dR8 <- gamma8 * I8 - (1-sigma) * lambda * R8 - mu * R8
		dI9 <- (1-sigma) * lambda * R8 - (gamma9 + mu) * I9
		dR9 <- gamma9 * I9 - (1-sigma) * lambda * R9 - mu * R9 + gamma10 * I10
		dI10 <- (1-sigma) * lambda * R9 - (gamma10 + mu) * I10
		
		list(c(dS1, dI1, dR1, dI2, dR2, dI3, dR3, dI4, dR4, dI5, dR5, dI6, dR6, dI7, dR7, dI8, dR8,
					 dI9, dR9, dI10),
				 prevalence=prevalence)
	})
}

simulate_cobey <- function(beta=1.1/(25/7),
													 kappa=25/7,
													 nu=220/7,
													 epsilon=0.4,
													 sigma=0.3,
													 mu=1/80/52,
													 yini,
													 tmin=0,
													 tmax=30) {
	parms <- c(
		beta=beta,
		kappa=kappa,
		nu=nu,
		epsilon=epsilon,
		sigma=sigma,
		mu=mu
	)
	
	R0 <- beta/(1/nu+mu)
	
	if (missing(yini)) {
		yini <- c(S1=1/R0-1e-6, I1=1e-6, R1=1-1/R0, I2=0, R2=0, I3=0, R3=0, I4=0, R4=0, I5=0, R5=0, 
							I6=0, R6=0, I7=0, R7=0, I8=0, R8=0, I9=0, R9=0, I10=0)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(ode(yini, times, model_cobey, parms))
	
	out
}

simulate_cobey_niche <- function(beta=1.1/(25/7),
																 kappa=25/7,
																 numin=25/7,
																 numax=220/7,
																 epsilon=0.4,
																 sigma=0.3,
																 mu=1/80/52,
																 mu_max=0.25,
																 tmin=0,
																 tmax=100,
																 Z=25) {
	nuvec <- seq(numax, numin, length.out=Z)
	
	simlist <- vector('list', length=Z-1)
	
	for (i in 1:(Z-1)) {
		tmplist <- vector('list', length=Z-i)
		
		for (j in (i+1):(Z)) {
			print(paste0(i, ", ", j))
			out1 <- tail(simulate_cobey(beta=beta, nu=nuvec[i], kappa=kappa, epsilon=epsilon, sigma=sigma, mu=mu,
																	tmin=tmin, tmax=tmax),1)
			
			nu1 <- (kappa + (nuvec[i] - kappa) * exp(-epsilon*(0:9)))
			
			out2 <- tail(simulate_cobey(beta=beta, nu=nuvec[j], kappa=kappa, epsilon=epsilon, sigma=sigma, mu=mu,
																	tmin=tmin, tmax=tmax),1)
			
			nu2 <- (kappa + (nuvec[j] - kappa) * exp(-epsilon*(0:9)))
			
			R01 <- beta/(1/nuvec[i]+mu)
			R02 <- beta/(1/nuvec[j]+mu)
			
			omega21 <- mu_max * (1-(i-1)/(Z-1))
			omega12 <- mu_max * (1-(j-1)/(Z-1))
			
			S_11 <- sum(c(out1$S1, 
										(1-sigma) * c(out1$R1, out1$R2, out1$R3, out1$R4, out1$R5,
																	out1$R6, out1$R7, out1$R8, out1$R9)) * (1/nu1[1]+mu)/(1/nu1+mu))
			S_21 <- sum(
				c(out1$S1, 
					c(out1$R1, out1$R2, out1$R3, out1$R4, out1$R5,
						out1$R6, out1$R7, out1$R8, out1$R9) * (1/nu2[1]+mu)/(1/nu2[-1]+mu),
					(1-omega21) * c(out1$I1, out1$I2, out1$I3, out1$I4, out1$I5,
													out1$I6, out1$I7, out1$I8, out1$I9, out1$I10) * (1/nu2[1]+mu)/(1/c(nu2[-1], nu2[10])+mu))
			)
			
			S_22 <- sum(c(out2$S1, 
										(1-sigma) * c(out2$R1, out2$R2, out2$R3, out2$R4, out2$R5,
																	out2$R6, out2$R7, out2$R8, out2$R9)) * (1/nu2[1]+mu)/(1/nu2+mu))
			S_12 <- sum(
				c(out2$S1, 
					c(out2$R1, out2$R2, out2$R3, out2$R4, out2$R5,
						out2$R6, out2$R7, out2$R8, out2$R9) * (1/nu1[1]+mu)/(1/nu1[-1]+mu),
					(1-omega12) * c(out2$I1, out2$I2, out2$I3, out2$I4, out2$I5,
													out2$I6, out2$I7, out2$I8, out2$I9, out2$I10) * (1/nu1[1]+mu)/(1/c(nu1[-1], nu1[10])+mu))
			)
			
			niche <- sqrt(S_11*S_22/(S_12*S_21))
			
			fitness <- R02/R01*sqrt(S_22*S_21/(S_11*S_12))
			
			tmplist[[j-i]] <- data.frame(
				nichediff=1/niche,
				fitnessdiff=1/fitness,
				system=paste0(i, ", ", j),
				meanR0=sqrt(R01*R02),
				nulldiff=sqrt(R01/R02)
			)
		}
		
		simlist[[i]] <- do.call("rbind", tmplist)
	}
	
	out <- do.call("rbind", simlist)
	
	out
}

model_cobey_vacc <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		prevalence <- I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I10
		lambda <- beta * prevalence
		
		gamma1 <- 1/(kappa + (nu - kappa) * exp(-0))
		gamma2 <- 1/(kappa + (nu - kappa) * exp(-epsilon))
		gamma3 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*2)))
		gamma4 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*3)))
		gamma5 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*4)))
		gamma6 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*5)))
		gamma7 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*6)))
		gamma8 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*7)))
		gamma9 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*8)))
		gamma10 <- 1/(kappa + (nu - kappa) * exp(-(epsilon*9)))
		
		dS1 <- mu - (1-p) * lambda * S1 - mu * S1
		dI1 <- (1-p) *lambda * S1 - (gamma1 + mu) * I1
		dR1 <- gamma1 * I1 - (1-p) * lambda * R1 - mu * R1
		dI2 <- (1-p) * lambda * R1 - (gamma2 + mu) * I2
		dR2 <- gamma2 * I2 - (1-p) * lambda * R2 - mu * R2
		dI3 <- (1-p) * lambda * R2 - (gamma3 + mu) * I3
		dR3 <- gamma3 * I3 - (1-p) * lambda * R3 - mu * R3
		dI4 <- (1-p) * lambda * R3 - (gamma4 + mu) * I4
		dR4 <- gamma4 * I4 - (1-p) * lambda * R4 - mu * R4
		dI5 <- (1-p) * lambda * R4 - (gamma5 + mu) * I5
		dR5 <- gamma5 * I5 - (1-p) * lambda * R5 - mu * R5
		dI6 <- (1-p) * lambda * R5 - (gamma6 + mu) * I6
		dR6 <- gamma6 * I6 - (1-p) * lambda * R6 - mu * R6
		dI7 <- (1-p) * lambda * R6 - (gamma7 + mu) * I7
		dR7 <- gamma7 * I7 - (1-p) * lambda * R7 - mu * R7
		dI8 <- (1-p) * lambda * R7 - (gamma8 + mu) * I8
		dR8 <- gamma8 * I8 - (1-p) * lambda * R8 - mu * R8
		dI9 <- (1-p) * lambda * R8 - (gamma9 + mu) * I9
		dR9 <- gamma9 * I9 - (1-p) * lambda * R9 - mu * R9 + gamma10 * I10
		dI10 <- (1-p) * lambda * R9 - (gamma10 + mu) * I10
		
		list(c(dS1, dI1, dR1, dI2, dR2, dI3, dR3, dI4, dR4, dI5, dR5, dI6, dR6, dI7, dR7, dI8, dR8,
					 dI9, dR9, dI10),
				 prevalence=prevalence)
	})
}

simulate_cobey_vacc <- function(beta=1.1/(25/7),
																kappa=25/7,
																nu=220/7,
																epsilon=0.4,
																p=0.6,
																mu=1/80/52,
																yini,
																tmin=0,
																tmax=30) {
	parms <- c(
		beta=beta,
		kappa=kappa,
		nu=nu,
		epsilon=epsilon,
		p=p,
		mu=mu
	)
	
	R0 <- beta/(1/nu+mu)
	
	if (missing(yini)) {
		yini <- c(S1=1/R0-1e-6, I1=1e-6, R1=1-1/R0, I2=0, R2=0, I3=0, R3=0, I4=0, R4=0, I5=0, R5=0, 
							I6=0, R6=0, I7=0, R7=0, I8=0, R8=0, I9=0, R9=0, I10=0)
	}
	
	times <- seq(52*tmin, 52*tmax, by=1)
	
	out <- as.data.frame(ode(yini, times, model_cobey_vacc, parms))
	
	out
}

simulate_cobey_vacc_niche <- function(beta=1.1/(25/7),
																			kappa=25/7,
																			numin=25/7,
																			numax=220/7,
																			epsilon=0.4,
																			sigma=0.3,
																			p=0.6,
																			mu=1/80/52,
																			mu_max=0.25,
																			tmin=0,
																			tmax=30,
																			Z=25,
																			vacc=c(1:7)) {
	if (sigma > p) {
		stop("sigma needs to be smaller than p")
	}
	
	nuvec <- seq(numax, numin, length.out=Z)
	
	simlist <- vector('list', length=Z-1)
	
	for (i in 1:(Z-1)) {
		tmplist <- vector('list', length=Z-i)
		
		for (j in (i+1):(Z)) {
			print(paste0(i, ", ", j))
			
			nu1 <- (kappa + (nuvec[i] - kappa) * exp(-epsilon*(0:9)))
			nu2 <- (kappa + (nuvec[j] - kappa) * exp(-epsilon*(0:9)))
			
			R01 <- beta/(1/nuvec[i]+mu)
			R02 <- beta/(1/nuvec[j]+mu)
			
			omega21 <- mu_max * (1-(i-1)/(Z-1))
			omega12 <- mu_max * (1-(j-1)/(Z-1))
			
			if (i %in% vacc) {
				out1 <- tail(simulate_cobey_vacc(beta=beta, nu=nuvec[i], kappa=kappa, epsilon=epsilon, mu=mu,
																				 tmin=tmin, tmax=tmax),1)	
				S_11 <- sum((1-p) * c(out1$S1, c(out1$R1, out1$R2, out1$R3, out1$R4, out1$R5,
																		out1$R6, out1$R7, out1$R8, out1$R9)) * (1/nu1[1]+mu)/(1/nu1+mu))
				
			} else {
				out1 <- tail(simulate_cobey(beta=beta, nu=nuvec[i], kappa=kappa, epsilon=epsilon, sigma=sigma, mu=mu,
																		tmin=tmin, tmax=tmax),1)
				
				S_11 <- sum(c(out1$S1, 
											(1-sigma) * c(out1$R1, out1$R2, out1$R3, out1$R4, out1$R5,
																		out1$R6, out1$R7, out1$R8, out1$R9)) * (1/nu1[1]+mu)/(1/nu1+mu))
			}
			
			if (j %in% vacc) {
				out2 <- tail(simulate_cobey_vacc(beta=beta, nu=nuvec[j], kappa=kappa, epsilon=epsilon, mu=mu,
																		tmin=tmin, tmax=tmax),1) 
				
				
				S_22 <- sum((1-p) * c(out2$S1, 
											c(out2$R1, out2$R2, out2$R3, out2$R4, out2$R5,
																		out2$R6, out2$R7, out2$R8, out2$R9)) * (1/nu2[1]+mu)/(1/nu2+mu))
				
				S_21 <- sum(
					c(out1$S1, 
						c(out1$R1, out1$R2, out1$R3, out1$R4, out1$R5,
							out1$R6, out1$R7, out1$R8, out1$R9) * (1/nu2[1]+mu)/(1/nu2[-1]+mu),
						(1-omega21) * c(out1$I1, out1$I2, out1$I3, out1$I4, out1$I5,
														out1$I6, out1$I7, out1$I8, out1$I9, out1$I10) * (1/nu2[1]+mu)/(1/c(nu2[-1], nu2[10])+mu))
				) * (1-p)
			} else {
				out2 <- tail(simulate_cobey(beta=beta, nu=nuvec[j], kappa=kappa, epsilon=epsilon, sigma=sigma, mu=mu,
																		tmin=tmin, tmax=tmax),1)
				
				
				S_22 <- sum(c(out2$S1, 
											(1-sigma) * c(out2$R1, out2$R2, out2$R3, out2$R4, out2$R5,
																		out2$R6, out2$R7, out2$R8, out2$R9)) * (1/nu2[1]+mu)/(1/nu2+mu))
				
				S_21 <- sum(
					c(out1$S1, 
						c(out1$R1, out1$R2, out1$R3, out1$R4, out1$R5,
							out1$R6, out1$R7, out1$R8, out1$R9) * (1/nu2[1]+mu)/(1/nu2[-1]+mu),
						(1-omega21) * c(out1$I1, out1$I2, out1$I3, out1$I4, out1$I5,
														out1$I6, out1$I7, out1$I8, out1$I9, out1$I10) * (1/nu2[1]+mu)/(1/c(nu2[-1], nu2[10])+mu))
				)
			}
			
			if (i %in% vacc) {
				S_12 <- sum(
					c(out2$S1, 
						c(out2$R1, out2$R2, out2$R3, out2$R4, out2$R5,
							out2$R6, out2$R7, out2$R8, out2$R9) * (1/nu1[1]+mu)/(1/nu1[-1]+mu),
						(1-omega12) * c(out2$I1, out2$I2, out2$I3, out2$I4, out2$I5,
														out2$I6, out2$I7, out2$I8, out2$I9, out2$I10) * (1/nu1[1]+mu)/(1/c(nu1[-1], nu1[10])+mu))
				) * (1-p)
			} else {
				S_12 <- sum(
					c(out2$S1, 
						c(out2$R1, out2$R2, out2$R3, out2$R4, out2$R5,
							out2$R6, out2$R7, out2$R8, out2$R9) * (1/nu1[1]+mu)/(1/nu1[-1]+mu),
						(1-omega12) * c(out2$I1, out2$I2, out2$I3, out2$I4, out2$I5,
														out2$I6, out2$I7, out2$I8, out2$I9, out2$I10) * (1/nu1[1]+mu)/(1/c(nu1[-1], nu1[10])+mu))
				)
			}
			
			niche <- sqrt(S_11*S_22/(S_12*S_21))
			
			fitness <- R02/R01*sqrt(S_22*S_21/(S_11*S_12))
			
			tmplist[[j-i]] <- data.frame(
				nichediff=1/niche,
				fitnessdiff=1/fitness,
				system=paste0(i, ", ", j),
				meanR0=sqrt(R01*R02),
				nulldiff=sqrt(R01/R02)
			)
		}
		
		simlist[[i]] <- do.call("rbind", tmplist)
	}
	
	out <- do.call("rbind", simlist)
	
	out
}



# R_s = R_p * (1-reduction)

## reduction in duration of infection
## R_p (I_p + I_s) * S_p = I_p 
## R_s (1-sigma) * (I_p + I_s) * S_s = I_s
## R_p * (1-reduction) (1-sigma) * (I_p + I_s) * S_s = I_s
## R_p * S_p + R_p * (1-reduction) * (1-sigma) * S_s = 1
## R_0 (S_p + (1-reduction) * (1-sigma) * S_s)

## reduction in transmission rate
## R_0 (I_p +  (1-epsilon) I_s) * S_p = I_p 
## R_0 * (I_p + (1-epsilon) I_s) * (1-sigma) * S_s = I_s
## R_0 * (I_p + (1-epsilon) I_s)/(I_p + I_s) * (S_p + (1-sigma) * S_s)

## reduction in all
## beta (sum_j (1-epsilon_j) I_j) * (1-sigma_i) * S_p = (gamma_i + mu) I_i
## R_0 * (sum_j (1-epsilon_j) I_j) * (1-reduciton_i) * (1-sigma_i) * S_p = I_i

## S_i -> I_i
## dI_i/dt = (1-sigma_i) lambda S_i - (gamma + mu) I_i
## During invasion, 
## I_i = I_i(0) exp(r)
## meaning
##  (1-sigma_i) lambda S_i - (gamma + mu) I_i = r I_i

simulate_koelle_niche <- function(beta=2.9*(7/9.5),
																	sigma=7/9.5,
																	mu=1/38.8/52,
																	gamma=0.956,
																	yini,
																	tmin=0,
																	tmax=100) {
	## dYI/dt = beta * XN * YI - (sigma+mu) Y_I
	## dXN/dt = mu - beta * XN * YI - mu * XN
	## dXI/dt = beta * XN * YI - mu * XI
	
	## XN* = (sigma+mu)/beta = 1/R0
	## mu (1/(sigma+mu)-1) = YI
	
	R0 <- beta/(mu+sigma)
	
	S_11 <- 1/R0
	S_12 <- 1/R0 + (1-gamma) * (1 - 1/R0)
	S_22 <- 1/R0
	S_21 <- 1/R0 + (1-gamma) * (1 - 1/R0)
	
	niche <- mean(sqrt(S_11*S_22/(S_12*S_21)))
	
	fitness <- mean(sqrt(S_22*S_21/(S_11*S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		system="Cholera Inaba, Ogawa",
		meanR0=sqrt(R0*R0),
		nulldiff=1
	)
}

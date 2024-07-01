## actually same model as pox!! not sure why I need a separate file
## but good to keep one I guess...

simulate_ferguson_niche <- function(R01=5,
																		sigma21=-0.5,
																		R02=2,
																		sigma12=-0.5,
																		system="DENV1, DENV2") {
	S_11 <- 1/R01
	S_21 <- 1/R01 + (1-1/R01) * (1-sigma21)
	
	S_22 <- 1/R02
	S_12 <- 1/R02 + (1-1/R02) * (1-sigma12)
	
	R_1 <- R01
	R_2 <- R02
	
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

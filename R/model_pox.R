# https://academic.oup.com/ije/article-abstract/17/3/643/729853
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7463189/
# https://www.bmj.com/content/379/bmj-2022-073153.full

simulate_pox_niche <- function(R0_small=5,
															 sigma_small=0.85,
															 R0_mpox=2,
															 sigma_mpox=0.5) {
	S_11 <- 1/R0_small
	S_21 <- 1/R0_small + (1-1/R0_small) * (1-sigma_small)
	
	S_22 <- 1/R0_mpox
	S_12 <- 1/R0_mpox + (1-1/R0_mpox) * (1-sigma_mpox)
	
	R_1 <- R0_small
	R_2 <- R0_mpox
	
	niche <- mean(sqrt(S_11*S_22/(S_12*S_21)))
	
	fitness <- mean(R_2/R_1*sqrt(S_22*S_21/(S_11*S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=pmax(fitness, 1/fitness),
		system="Smallpox, mpox",
		meanR0=mean(sqrt(R_1*R_2)),
		nulldiff=pmax(mean(sqrt(R_1/R_2)), mean(sqrt(R_2/R_1)))
	)
}

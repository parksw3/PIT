library(dplyr)
source("../R/model_pox.R")

R0_small_vec <- seq(3.5, 6, length.out=11)
sigma_small_vec <- seq(0.8, 0.95, length.out=11)
R0_mpox_vec <- seq(1.5, 2.7, length.out=11)
sigma_mpox_vec <- seq(0.1, 0.9, length.out=11)

paramdata <- expand.grid(R0_small_vec,
												 sigma_small_vec,
												 R0_mpox_vec,
												 sigma_mpox_vec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
	reslist[[i]] <- simulate_pox_niche(
		R0_small = paramdata[i,1],
		sigma_small = paramdata[i,2],
		R0_mpox = paramdata[i,3],
		sigma_mpox = paramdata[i,4]
	)
}

simulate_niche_pox <- reslist %>%
	bind_rows

save("simulate_niche_pox", file="simulate_niche_pox.rda")

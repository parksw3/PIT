library(deSolve)
source("../R/model_cobey.R")

data_cobey_niche_base <- simulate_cobey_niche(sigma=0, epsilon=0, beta=1.05/(25/7))
data_cobey_niche_anti <- simulate_cobey_niche(sigma=0.3, epsilon=0, beta=1.05/(25/7))
data_cobey_niche_nons <- simulate_cobey_niche(sigma=0, epsilon=0.4, beta=1.05/(25/7))
data_cobey_niche_both <- simulate_cobey_niche(sigma=0.3, epsilon=0.4, beta=1.05/(25/7))

save("data_cobey_niche_base", "data_cobey_niche_anti",
		 "data_cobey_niche_nons", "data_cobey_niche_both",
		 file="simulate_niche_cobey.rda")

library(dplyr)
source("../R/model_ferguson.R")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1692557/pdf/10365401.pdf
R0 <- c(4.29, 5.12, 4.61, 5.50)

reslist <- vector('list', 6)

count <- 1

for (i in 1:3) {
	for (j in (i+1):4) {
		system <- paste0("DENV", i, ", DENV", j)
		
		reslist[[count]] <- simulate_ferguson_niche(R01=R0[i], R02=R0[j],
																								system=system)
		
		count <- count + 1
	}
}

simulate_niche_dengue <- reslist %>%
	bind_rows

save("simulate_niche_dengue", file="simulate_niche_dengue.rda")

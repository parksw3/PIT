source("../R/model_meijers.R")

# https://www.nature.com/articles/s41591-021-01377-8
# https://github.com/m-meijers/Population_Immunity_SARSCOV2/blob/main/Create_Immune_Trajectory_Data.py
# https://www.medrxiv.org/content/10.1101/2021.12.19.21268038v1.full.pdf

lambda <- 3
T50 <- log10(0.2*94)
T0 <- log10(94)

R0delta <- c(3, 4, 5, 6)
R0ratio <- c(1.5, 2, 3)
cross <- c(0.2, 0.4, 0.6)

paramdata <- expand.grid(R0delta, R0ratio, cross)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
	print(i)
	tmplist <- vector('list', 4)
	
	pp <- paramdata[i,]
	
	TT <- 10^(-(T50-log(1/pp[[3]]-1)/lambda-T0))

	tmplist[[1]] <- simulate_meijers_vacc_niche(R_1=pp[[1]],
																							T11=1,
																							T21=TT,
																							tau1=3.5,
																							R_2=pp[[1]]*pp[[2]],
																							T22=1,
																							T12=TT,
																							tau2=3,
																							vprop=0,
																							bprop=0,
																							TV1=3.2,
																							TV2=47,
																							TB1=2.8,
																							TB2=6.7,
																							system="SARS-CoV-2 Delta, BA.1") %>%
		mutate(
			type="unvaccinated",
			R0delta=pp[[1]],
			Rratio=pp[[2]],
			cross=pp[[3]]
		)
	
	tmplist[[2]] <- simulate_meijers_vacc_niche(R_1=pp[[1]],
																							T11=1,
																							T21=TT,
																							tau1=3.5,
																							R_2=pp[[1]]*pp[[2]],
																							T22=1,
																							T12=TT,
																							tau2=3,
																							vprop=0.3,
																							bprop=0,
																							TV1=3.2,
																							TV2=47,
																							TB1=2.8,
																							TB2=6.7,
																							system="SARS-CoV-2 Delta, BA.1") %>%
		mutate(
			type="30% vac",
			R0delta=pp[[1]],
			Rratio=pp[[2]],
			cross=pp[[3]]
		)
		
	tmplist[[3]] <- simulate_meijers_vacc_niche(R_1=pp[[1]],
																							T11=1,
																							T21=TT,
																							tau1=3.5,
																							R_2=pp[[1]]*pp[[2]],
																							T22=1,
																							T12=TT,
																							tau2=3,
																							vprop=0.4,
																							bprop=0.3,
																							TV1=3.2,
																							TV2=47,
																							TB1=2.8,
																							TB2=6.7,
																							system="SARS-CoV-2 Delta, BA.1") %>%
		mutate(
			type="40% vac, 30% bst",
			R0delta=pp[[1]],
			Rratio=pp[[2]],
			cross=pp[[3]]
		)
	
	tmplist[[4]] <- simulate_meijers_vacc_niche(R_1=pp[[1]],
																							T11=1,
																							T21=TT,
																							tau1=3.5,
																							R_2=pp[[1]]*pp[[2]],
																							T22=1,
																							T12=TT,
																							tau2=3,
																							vprop=0.2,
																							bprop=0.5,
																							TV1=3.2,
																							TV2=47,
																							TB1=2.8,
																							TB2=6.7,
																							system="SARS-CoV-2 Delta, BA.1") %>%
		mutate(
			type="20% vac, 50% bst",
			R0delta=pp[[1]],
			Rratio=pp[[2]],
			cross=pp[[3]]
		)
	
	reslist[[i]] <- tmplist %>%
		bind_rows
}

simulate_niche_covid_omicron_sensitivity <- reslist %>%
	bind_rows

save("simulate_niche_covid_omicron_sensitivity", file="simulate_niche_covid_omicron_sensitivity.rda")

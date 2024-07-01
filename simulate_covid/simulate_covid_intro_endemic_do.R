source("../R/model_meijers.R")

lambda <- 3
T50 <- log10(0.2*94)
T0 <- log10(94)

R0 <- 2.5
R02 <- 2.5 * 1.7
R03 <- 2.5 * 1.7 * 1.5
R04 <- 2.5 * 1.7 * 1.5 * 1.5
cross <- 0.4
TT <- 10^(-(T50-log(1/cross-1)/lambda-T0))

1/(exp(-((T0-log10(27)-T50)*lambda))+1)

simulate_covid_delta_BA1_sens <- simulate_meijers_nstrain(R0=c(R03, R04),
																													Tmat=matrix(c(1, TT, TT, 1), 2, 2),
																													tau=c(3.5, 3),
																													i0=c(1e-6, 1e-6),
																													intro=c(1, 365*50),
																													tmax=365*55)

simulate_covid_delta_BA1_vacc_sens <- simulate_meijers_nstrain_vacc(R0=c(R03, R04),
																																		Tmat=matrix(c(1, TT, TT, 1), 2, 2),
																																		tau=c(3.5, 3),
																																		i0=c(1e-6, 1e-6),
																																		intro=c(1, 365*50),
																																		tmax=365*55,
																																		TV=c(3.2, 47),
																																		TB=c(2.8, 6.7),
																																		vprop=0.3,
																																		bprop=0)

simulate_covid_delta_BA1_boos_sens <- simulate_meijers_nstrain_vacc(R0=c(R03, R04),
																																		Tmat=matrix(c(1, TT, TT, 1), 2, 2),
																																		tau=c(3.5, 3),
																																		i0=c(1e-6, 1e-6),
																																		intro=c(1, 365*50),
																																		tmax=365*55,
																																		TV=c(3.2, 47),
																																		TB=c(2.8, 6.7),
																																		vprop=0.4,
																																		bprop=0.3)

simulate_covid_delta_BA1_more_sens <- simulate_meijers_nstrain_vacc(R0=c(R03, R04),
																																		Tmat=matrix(c(1, TT, TT, 1), 2, 2),
																																		tau=c(3.5, 3),
																																		i0=c(1e-6, 1e-6),
																																		intro=c(1, 365*50),
																																		tmax=365*55,
																																		TV=c(3.2, 47),
																																		TB=c(2.8, 6.7),
																																		vprop=0.2,
																																		bprop=0.5)

save("simulate_covid_delta_BA1_sens", 
		 "simulate_covid_delta_BA1_vacc_sens",
		 "simulate_covid_delta_BA1_boos_sens",
		 "simulate_covid_delta_BA1_more_sens",
		 file="simulate_covid_intro_endemic_do.rda")

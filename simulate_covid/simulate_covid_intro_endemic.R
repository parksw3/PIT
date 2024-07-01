source("../R/model_meijers.R")

R0 <- 2.5
F0 <- c(0.05212319, 0.02772249, 0.05731528, 0.07480703, 0.0, 0.0)

simulate_covid_wt_alpha <- simulate_meijers_nstrain(R0=c(R0, R0*exp(F0[1]*5.5)),
																										Tmat=matrix(c(1, 1.8, 1.8, 1), 2, 2),
																										tau=c(7.8, 5.5),
																										i0=c(1e-6, 1e-6),
																										intro=c(1, 365*50),
																										tmax=365*55)

simulate_covid_alpha_wt <- simulate_meijers_nstrain(R0=c(R0, R0*exp(F0[1]*5.5)),
																										Tmat=matrix(c(1, 1.8, 1.8, 1), 2, 2),
																										tau=c(7.8, 5.5),
																										i0=c(1e-6, 1e-6),
																										intro=c(365*50, 1),
																										tmax=365*55)

simulate_covid_alpha_delta <- simulate_meijers_nstrain(R0=c(R0*exp(F0[1]*5.5), R0*exp((F0[1]+F0[2])*3.5)),
																											 Tmat=matrix(c(1, 2.8, 3.5, 1), 2, 2),
																											 tau=c(5.5, 3.5),
																											 i0=c(1e-6, 1e-6),
																											 intro=c(1, 365*50),
																											 tmax=365*55)

simulate_covid_delta_alpha <- simulate_meijers_nstrain(R0=c(R0*exp(F0[1]*5.5), R0*exp((F0[1]+F0[2])*3.5)),
																											 Tmat=matrix(c(1, 2.8, 3.5, 1), 2, 2),
																											 tau=c(5.5, 3.5),
																											 i0=c(1e-6, 1e-6),
																											 intro=c(365*50, 1),
																											 tmax=365*55)

simulate_covid_delta_BA1 <- simulate_meijers_nstrain(R0=c(R0*exp((F0[1]+F0[2])*3.5), R0*exp((F0[1]+F0[2]+F0[3])*3)),
																										 Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																										 tau=c(3.5, 3),
																										 i0=c(1e-6, 1e-6),
																										 intro=c(1, 365*50),
																										 tmax=365*55)

simulate_covid_BA1_delta <- simulate_meijers_nstrain(R0=c(R0*exp((F0[1]+F0[2])*3.5), R0*exp((F0[1]+F0[2]+F0[3])*3)),
																										 Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																										 tau=c(3.5, 3),
																										 i0=c(1e-6, 1e-6),
																										 intro=c(365*50, 1),
																										 tmax=365*55)

simulate_covid_BA1_BA2 <- simulate_meijers_nstrain(R0=c(R0*exp((F0[1]+F0[2]+F0[3])*3), R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3)),
																									 Tmat=matrix(c(1, 2.4, 3.6, 1), 2, 2),
																									 tau=c(3, 3),
																									 i0=c(1e-6, 1e-6),
																									 intro=c(1, 365*50),
																									 tmax=365*55)

simulate_covid_BA2_BA1 <- simulate_meijers_nstrain(R0=c(R0*exp((F0[1]+F0[2]+F0[3])*3), R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3)),
																									 Tmat=matrix(c(1, 2.4, 3.6, 1), 2, 2),
																									 tau=c(3, 3),
																									 i0=c(1e-6, 1e-6),
																									 intro=c(365*50, 1),
																									 tmax=365*55)

save("simulate_covid_wt_alpha", "simulate_covid_alpha_wt",
		 "simulate_covid_alpha_delta", "simulate_covid_delta_alpha",
		 "simulate_covid_delta_BA1", "simulate_covid_BA1_delta",
		 "simulate_covid_BA1_BA2", "simulate_covid_BA2_BA1",
		 file="simulate_covid_intro_endemic.rda")

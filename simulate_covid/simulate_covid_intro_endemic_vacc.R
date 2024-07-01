source("../R/model_meijers.R")

R0 <- 2.5
F0 <- c(0.05212319, 0.02772249, 0.05731528, 0.07480703, 0.0, 0.0)

simulate_covid_wt_alpha_vacc <- simulate_meijers_nstrain_vacc(R0=c(R0, R0*exp(F0[1]*5.5)),
																															Tmat=matrix(c(1, 1.8, 1.8, 1), 2, 2),
																															tau=c(7.8, 5.5),
																															i0=c(1e-6, 1e-6),
																															intro=c(1, 365*50),
																															tmax=365*55,
																															TV=c(1, 1.8),
																															TB=c(1, 1.8),
																															vprop=0.3,
																															bprop=0)

simulate_covid_alpha_delta_vacc <- simulate_meijers_nstrain_vacc(R0=c(R0*exp(F0[1]*5.5), R0*exp((F0[1]+F0[2])*3.5)),
																																 Tmat=matrix(c(1, 2.8, 3.5, 1), 2, 2),
																																 tau=c(5.5, 3.5),
																																 i0=c(1e-6, 1e-6),
																																 intro=c(1, 365*50),
																																 tmax=365*55,
																																 TV=c(1.8, 3.2),
																																 TB=c(1.8, 2.8),
																																 vprop=0.3,
																																 bprop=0)

simulate_covid_alpha_delta_vacc2 <- simulate_meijers_nstrain_vacc(R0=c(R0*1.7, R0*1.7*1.5),
																																	Tmat=matrix(c(1, 2.8, 3.5, 1), 2, 2),
																																	tau=c(5.5, 3.5),
																																	i0=c(1e-6, 1e-6),
																																	intro=c(1, 365*50),
																																	tmax=365*55,
																																	TV=c(1.8, 3.2),
																																	TB=c(1.8, 2.8),
																																	vprop=0.3,
																																	bprop=0)

simulate_covid_alpha_delta_boos <- simulate_meijers_nstrain_vacc(R0=c(R0*exp(F0[1]*5.5), R0*exp((F0[1]+F0[2])*3.5)),
																																 Tmat=matrix(c(1, 2.8, 3.5, 1), 2, 2),
																																 tau=c(5.5, 3.5),
																																 i0=c(1e-6, 1e-6),
																																 intro=c(1, 365*50),
																																 tmax=365*55,
																																 TV=c(1.8, 3.2),
																																 TB=c(1.8, 2.8),
																																 vprop=0.4,
																																 bprop=0.3)

simulate_covid_alpha_delta_boos2 <- simulate_meijers_nstrain_vacc(R0=c(R0*1.7, R0*1.7*1.5),
																																	Tmat=matrix(c(1, 2.8, 3.5, 1), 2, 2),
																																	tau=c(5.5, 3.5),
																																	i0=c(1e-6, 1e-6),
																																	intro=c(1, 365*50),
																																	tmax=365*55,
																																	TV=c(1.8, 3.2),
																																	TB=c(1.8, 2.8),
																																	vprop=0.4,
																																	bprop=0.3)

simulate_covid_delta_BA1_vacc <- simulate_meijers_nstrain_vacc(R0=c(R0*exp((F0[1]+F0[2])*3.5), R0*exp((F0[1]+F0[2]+F0[3])*3)),
																															 Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																															 tau=c(3.5, 3),
																															 i0=c(1e-6, 1e-6),
																															 intro=c(1, 365*50),
																															 tmax=365*55,
																															 TV=c(3.2, 47),
																															 TB=c(2.8, 6.7),
																															 vprop=0.3,
																															 bprop=0)

simulate_covid_delta_BA1_vacc2 <- simulate_meijers_nstrain_vacc(R0=c(R0*1.7*1.5, R0*1.7*1.5*1.5),
																																Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																																tau=c(3.5, 3),
																																i0=c(1e-6, 1e-6),
																																intro=c(1, 365*50),
																																tmax=365*55,
																																TV=c(3.2, 47),
																																TB=c(2.8, 6.7),
																																vprop=0.3,
																																bprop=0)


simulate_covid_delta_BA1_boos <- simulate_meijers_nstrain_vacc(R0=c(R0*exp((F0[1]+F0[2])*3.5), R0*exp((F0[1]+F0[2]+F0[3])*3)),
																															 Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																															 tau=c(3.5, 3),
																															 i0=c(1e-6, 1e-6),
																															 intro=c(1, 365*50),
																															 tmax=365*55,
																															 TV=c(3.2, 47),
																															 TB=c(2.8, 6.7),
																															 vprop=0.4,
																															 bprop=0.3)

simulate_covid_delta_BA1_boos2 <- simulate_meijers_nstrain_vacc(R0=c(R0*1.7*1.5, R0*1.7*1.5*1.5),
																																Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																																tau=c(3.5, 3),
																																i0=c(1e-6, 1e-6),
																																intro=c(1, 365*50),
																																tmax=365*55,
																																TV=c(3.2, 47),
																																TB=c(2.8, 6.7),
																																vprop=0.4,
																																bprop=0.3)

simulate_covid_delta_BA1_more <- simulate_meijers_nstrain_vacc(R0=c(R0*exp((F0[1]+F0[2])*3.5), R0*exp((F0[1]+F0[2]+F0[3])*3)),
																															 Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																															 tau=c(3.5, 3),
																															 i0=c(1e-6, 1e-6),
																															 intro=c(1, 365*50),
																															 tmax=365*55,
																															 TV=c(3.2, 47),
																															 TB=c(2.8, 6.7),
																															 vprop=0.2,
																															 bprop=0.5)

simulate_covid_delta_BA1_more2 <- simulate_meijers_nstrain_vacc(R0=c(R0*1.7*1.5, R0*1.7*1.5*1.5),
																															  Tmat=matrix(c(1, 27, 27, 1), 2, 2),
																															  tau=c(3.5, 3),
																															  i0=c(1e-6, 1e-6),
																															  intro=c(1, 365*50),
																															  tmax=365*55,
																															  TV=c(3.2, 47),
																															  TB=c(2.8, 6.7),
																															  vprop=0.2,
																															  bprop=0.5)

simulate_covid_BA1_BA2_boos <- simulate_meijers_nstrain_vacc(R0=c(R0*exp((F0[1]+F0[2]+F0[3])*3), R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3)),
																														 Tmat=matrix(c(1, 2.4, 3.6, 1), 2, 2),
																														 tau=c(3, 3),
																														 i0=c(1e-6, 1e-6),
																														 intro=c(1, 365*50),
																														 tmax=365*55,
																														 TV=c(47, 47),
																														 TB=c(6.7, 5.2),
																														 vprop=0.4,
																														 bprop=0.3)

simulate_covid_BA1_BA2_more <- simulate_meijers_nstrain_vacc(R0=c(R0*exp((F0[1]+F0[2]+F0[3])*3), R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3)),
																														 Tmat=matrix(c(1, 2.4, 3.6, 1), 2, 2),
																														 tau=c(3, 3),
																														 i0=c(1e-6, 1e-6),
																														 intro=c(1, 365*50),
																														 tmax=365*55,
																														 TV=c(47, 47),
																														 TB=c(6.7, 5.2),
																														 vprop=0.2,
																														 bprop=0.5)

save("simulate_covid_wt_alpha_vacc", 
		 "simulate_covid_alpha_delta_vacc",
		 "simulate_covid_alpha_delta_vacc",
		 "simulate_covid_alpha_delta_vacc2",
		 "simulate_covid_alpha_delta_boos",
		 "simulate_covid_alpha_delta_boos2",
		 "simulate_covid_delta_BA1_vacc",
		 "simulate_covid_delta_BA1_vacc2",
		 "simulate_covid_delta_BA1_boos",
		 "simulate_covid_delta_BA1_boos2",
		 "simulate_covid_delta_BA1_more",
		 "simulate_covid_delta_BA1_more2",
		 "simulate_covid_BA1_BA2_boos",
		 "simulate_covid_BA1_BA2_more",
		 file="simulate_covid_intro_endemic_vacc.rda")

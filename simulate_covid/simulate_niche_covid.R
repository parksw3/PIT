source("../R/model_meijers.R")

# https://www.nature.com/articles/s41591-021-01377-8
# https://github.com/m-meijers/Population_Immunity_SARSCOV2/blob/main/Create_Immune_Trajectory_Data.py

R0 <- 2.5
F0 <- c(0.05212319, 0.02772249, 0.05731528, 0.07480703, 0.0, 0.0)

data_meijers_niche1 <- simulate_meijers_niche(R_1=R0,
																							T11=1,
																							T21=1.8,
																							tau1=7.8,
																							R_2=R0*exp(F0[1]*5.5),
																							T22=1,
																							T12=1.8,
																							tau2=5.5,
																							system="SARS-CoV-2 WT, Alpha")

data_meijers_niche2 <- simulate_meijers_niche(R_1=R0*exp(F0[1]*5.5),
																							T11=1,
																							T21=2.8,
																							tau1=5.5,
																							R_2=R0*exp((F0[1]+F0[2])*3.5),
																							T22=1,
																							T12=3.5,
																							tau2=3.5,
																							system="SARS-CoV-2 Alpha, Delta")

data_meijers_niche3 <- simulate_meijers_niche(R_1=R0*exp((F0[1]+F0[2])*3.5),
																							T11=1,
																							T21=27,
																							tau1=3.5,
																							R_2=R0*exp((F0[1]+F0[2]+F0[3])*3),
																							T22=1,
																							T12=27,
																							tau2=3,
																							system="SARS-CoV-2 Delta, BA.1")

data_meijers_niche4 <- simulate_meijers_niche(R_1=R0*exp((F0[1]+F0[2]+F0[3])*3),
																							T11=1,
																							T21=2.4,
																							tau1=3,
																							R_2=R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3),
																							T22=1,
																							T12=3.6,
																							tau2=3,
																							system="SARS-CoV-2 BA.1, BA.2")

save("data_meijers_niche1", "data_meijers_niche2", "data_meijers_niche3",
		 "data_meijers_niche4",
		 file="simulate_niche_covid.rda")

source("../R/model_meijers.R")

# https://www.nature.com/articles/s41591-021-01377-8
# https://github.com/m-meijers/Population_Immunity_SARSCOV2/blob/main/Create_Immune_Trajectory_Data.py

R0 <- 2.5
F0 <- c(0.05212319, 0.02772249, 0.05731528, 0.07480703, 0.0, 0.0)

data_meijers_vacc_niche1 <- simulate_meijers_vacc_niche(R_1=R0,
																												T11=1,
																												T21=1.8,
																												tau1=7.8,
																												R_2=R0*exp(F0[1]*5.5),
																												T22=1,
																												T12=1.8,
																												tau2=5.5,
																												vprop=0.3,
																												bprop=0,
																												TV1=1,
																												TV2=1.8,
																												TB1=1,
																												TB2=1.8,
																												system="SARS-CoV-2 WT, Alpha") %>%
	mutate(
		type="30% vac"
	)

data_meijers_vacc_niche2 <- simulate_meijers_vacc_niche(R_1=R0*exp(F0[1]*5.5),
																												T11=1,
																												T21=2.8,
																												tau1=5.5,
																												R_2=R0*exp((F0[1]+F0[2])*3.5),
																												T22=1,
																												T12=3.5,
																												tau2=3.5,
																												vprop=0.3,
																												bprop=0,
																												TV1=1.8,
																												TV2=3.2,
																												TB1=1.8,
																												TB2=2.8,
																												system="SARS-CoV-2 Alpha, Delta") %>%
	mutate(
		type="30% vac"
	)

data_meijers_vacc_niche3 <- simulate_meijers_vacc_niche(R_1=R0*exp(F0[1]*5.5),
																												T11=1,
																												T21=2.8,
																												tau1=5.5,
																												R_2=R0*exp((F0[1]+F0[2])*3.5),
																												T22=1,
																												T12=3.5,
																												tau2=3.5,
																												vprop=0.4,
																												bprop=0.3,
																												TV1=1.8,
																												TV2=3.2,
																												TB1=1.8,
																												TB2=2.8,
																												system="SARS-CoV-2 Alpha, Delta") %>%
	mutate(
		type="40% vac, 30% bst"
	)

data_meijers_vacc_niche4 <- simulate_meijers_vacc_niche(R_1=R0*exp((F0[1]+F0[2])*3.5),
																												T11=1,
																												T21=27,
																												tau1=3.5,
																												R_2=R0*exp((F0[1]+F0[2]+F0[3])*3),
																												T22=1,
																												T12=27,
																												tau2=3,
																												vprop=0.3,
																												bprop=0,
																												TV1=3.2,
																												TV2=47,
																												TB1=2.8,
																												TB2=6.7,
																												system="SARS-CoV-2 Delta, BA.1") %>%
	mutate(
		type="30% vac"
	)

data_meijers_vacc_niche5 <- simulate_meijers_vacc_niche(R_1=R0*exp((F0[1]+F0[2])*3.5),
																												T11=1,
																												T21=27,
																												tau1=3.5,
																												R_2=R0*exp((F0[1]+F0[2]+F0[3])*3),
																												T22=1,
																												T12=27,
																												tau2=3,
																												vprop=0.4,
																												bprop=0.3,
																												TV1=3.2,
																												TV2=47,
																												TB1=2.8,
																												TB2=6.7,
																												system="SARS-CoV-2 Delta, BA.1") %>%
	mutate(
		type="40% vac, 30% bst"
	)

data_meijers_vacc_niche6 <- simulate_meijers_vacc_niche(R_1=R0*exp((F0[1]+F0[2])*3.5),
																												T11=1,
																												T21=27,
																												tau1=3.5,
																												R_2=R0*exp((F0[1]+F0[2]+F0[3])*3),
																												T22=1,
																												T12=27,
																												tau2=3,
																												vprop=0.2,
																												bprop=0.5,
																												TV1=3.2,
																												TV2=47,
																												TB1=2.8,
																												TB2=6.7,
																												system="SARS-CoV-2 Delta, BA.1") %>%
	mutate(
		type="20% vac, 50% bst"
	)

data_meijers_vacc_niche7 <- simulate_meijers_vacc_niche(R_1=R0*exp((F0[1]+F0[2]+F0[3])*3),
																												T11=1,
																												T21=2.4,
																												tau1=3,
																												R_2=R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3),
																												T22=1,
																												T12=3.6,
																												tau2=3,
																												vprop=0.4,
																												bprop=0.3,
																												TV1=47,
																												TV2=47,
																												TB1=6.7,
																												TB2=5.2,
																												system="SARS-CoV-2 BA.1, BA.2") %>%
	mutate(
		type="40% vac, 30% bst"
	)

data_meijers_vacc_niche8 <- simulate_meijers_vacc_niche(R_1=R0*exp((F0[1]+F0[2]+F0[3])*3),
																												T11=1,
																												T21=2.4,
																												tau1=3,
																												R_2=R0*exp((F0[1]+F0[2]+F0[3]+F0[4])*3),
																												T22=1,
																												T12=3.6,
																												tau2=3,
																												vprop=0.2,
																												bprop=0.5,
																												TV1=47,
																												TV2=47,
																												TB1=6.7,
																												TB2=5.2,
																												system="SARS-CoV-2 BA.1, BA.2") %>%
	mutate(
		type="20% vac, 50% bst"
	)

save("data_meijers_vacc_niche1", "data_meijers_vacc_niche2", "data_meijers_vacc_niche3",
		 "data_meijers_vacc_niche4", "data_meijers_vacc_niche5", "data_meijers_vacc_niche6",
		 "data_meijers_vacc_niche7", "data_meijers_vacc_niche8",
		 file="simulate_niche_covid_vacc.rda")

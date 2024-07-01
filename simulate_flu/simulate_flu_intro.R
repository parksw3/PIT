library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../R/model_gog.R")

load("simulate_flu_base.rda")

x <- 166
y <- 182

yini_base <- tail(head(unlist(base_sim[base_sim$time==(12*365+185),]), -2), -1)

yini_intro1 <- yini_base
yini_intro1[300+x] <- 1e-6
simulate_flu_intro1 <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																		gamma=1/5, mu=1/30/365,
																		s=0.07, n=300, nscale=10,
																		yini=yini_intro1, tvec=(13*365-180):(13*365+185))

yini_intro1_intro2 <- tail(head(unlist(tail(simulate_flu_intro1,1)), -2), -1)
yini_intro1_intro2[300+y] <- 1e-6
simulate_flu_intro1_intro2 <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																					 gamma=1/5, mu=1/30/365,
																					 s=0.07, n=300, nscale=10,
																					 yini=yini_intro1_intro2, tvec=(14*365-180):(17*365+185))

yini_intro2 <- yini_base
yini_intro2[300+x] <- 0
yini_intro2[300+y] <- 1e-6
simulate_flu_intro2 <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																		gamma=1/5, mu=1/30/365,
																		s=0.07, n=300, nscale=10,
																		yini=yini_intro2, tvec=(13*365-180):(13*365+185))

yini_intro2_intro1 <- tail(head(unlist(tail(simulate_flu_intro2,1)), -2), -1)
yini_intro2_intro1[300+x] <- 1e-6
simulate_flu_intro2_intro1 <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																					 gamma=1/5, mu=1/30/365,
																					 s=0.07, n=300, nscale=10,
																					 yini=yini_intro2_intro1, tvec=(14*365-180):(17*365+185))

simulate_flu_intro_x <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																		 gamma=1/5, mu=1/30/365,
																		 s=0.07, n=300, nscale=10,
																		 yini=yini_intro1, tvec=(13*365-180):(17*365+185))

simulate_flu_intro_y <- simulate_gog(R0=1.8, Reff=1.28, phi=0.15,
																		 gamma=1/5, mu=1/30/365,
																		 s=0.07, n=300, nscale=10,
																		 yini=yini_intro2, tvec=(13*365-180):(17*365+185))

simulate_flu_intro_xy <- bind_rows(
	head(simulate_flu_intro1,-1),
	simulate_flu_intro1_intro2
)

simulate_flu_intro_yx <- bind_rows(
	head(simulate_flu_intro2,-1),
	simulate_flu_intro2_intro1
)

save("simulate_flu_intro_xy", "simulate_flu_intro_yx",
		 "simulate_flu_intro_x", "simulate_flu_intro_y",
		 file="simulate_flu_intro.rda")

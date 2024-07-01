library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/model_gog.R")

load("../simulate_flu/simulate_flu_base.rda")
load("../simulate_flu/simulate_flu_intro_endemic.rda")
load("../simulate_flu/simulate_flu_intro.rda")
load("../simulate_flu/simulate_flu_intro_factorial.rda")
load("../simulate_flu/simulate_flu_intro_factorial_endemic.rda")

which_strain <- base_sim[,c(1, 302:601)] %>%
	gather(key, value, -time) %>%
	group_by(key) %>%
	filter(value > 0) %>%
	filter(time==min(time)) %>%
	mutate(
		key=gsub("I", "", key),
		key=as.numeric(key)
	)

base_sim_sus <- base_sim[,1:301] %>%
	gather(key, value, -time) %>%
	mutate(
		key=gsub("S", "", key),
		key=as.numeric(key)
	) %>%
	filter(key %in% which_strain$key) %>%
	mutate(
		antigenic=key/10,
		group=cut(antigenic, antigenic_cut)
	)

base_sim_sus_mean <- base_sim_sus %>%
	group_by(group, time) %>%
	summarize(
		value=mean(value)
	)

base_prev_gather2 <- base_prev_gather %>%
	mutate(
		antigenic=key/10,
		group=cut(antigenic, antigenic_cut)
	) %>%
	group_by(time, group) %>%
	summarize(
		total=sum(value)
	)

g1 <- ggplot(filter(base_prev_gather2, total > 1e-6)) +
	geom_line(aes(time/365, y=total, group=group), lwd=1) +
	geom_ribbon(ymin=-Inf, aes(time/365, ymax=total, fill=group)) +
	annotate("text", x=13.3, y=0.008, label="strain x", family="Times") +
	annotate("text", x=15.1, y=0.011, label="strain y", family="Times") +
	scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#0072B2", "#CC79A7")) +
	scale_x_continuous("Year", limits=c(-1, 21), expan=c(0,0)) +
	scale_y_continuous("Prevalence", limits=c(0, 0.047), expand=c(0, 0)) +
	ggtitle("A") +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g2 <- ggplot(simulate_flu_intro_endemic_xy) +
	geom_line(aes(time/365, I166), col="#D55E00", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_endemic_x, time <= 550), aes(time/365, I166), col="#D55E00", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_endemic_xy, I182 > 0), aes(time/365, I182), col="#0072B2", lwd=1)  +
	annotate("point", x=1.5, y=1e-6, shape=21, size=3, stroke=1.5, col="#0072B2", fill="white") +
	scale_x_continuous("Year", expand=c(0,0), limits=c(NA, 7+185/365),
										 breaks=1:10) +
	scale_y_log10("Prevalence") +
	ggtitle("B") +
	coord_cartesian(ylim=c(1e-8, 0.1)) +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g3 <- ggplot(simulate_flu_intro_endemic_yx) +
	geom_line(aes(time/365, I182), col="#0072B2", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_endemic_y, time <= 550), aes(time/365, I182), col="#0072B2", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_endemic_yx, I166 > 0), aes(time/365, I166), col="#D55E00", lwd=1)  +
	annotate("point", x=1.5, y=1e-6, shape=21, size=3, stroke=1.5, col="#D55E00", fill="white") +
	scale_x_continuous("Year", expand=c(0,0), limits=c(NA, 7+185/365),
										 breaks=1:10) +
	scale_y_log10("Prevalence") +
	ggtitle("C") +
	coord_cartesian(ylim=c(1e-8, 0.1)) +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g6 <- ggplot(base_sim_sus) +
	geom_line(aes(time/365, value, group=key, col=group), alpha=0.2) +
	geom_line(data=base_sim_sus_mean, aes(time/365, value, col=group), lwd=1) +
	geom_hline(yintercept=1/1.8, lty=2) +
	scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#0072B2", "#CC79A7")) +
	scale_x_continuous("Year", limits=c(-1, 21), expand=c(0,0)) +
	scale_y_continuous("Proportion susceptible", limits=c(0.3, 1), expand=c(0, 0)) +
	ggtitle("D") +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	) 

g7 <- ggplot(simulate_flu_intro_endemic_xy) +
	geom_line(data=filter(simulate_flu_intro_endemic_x, time <= 550), aes(time/365, S166), col="#D55E00", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_endemic_x, time <= 550), aes(time/365, S182), col="#0072B2", lwd=1) +
	geom_line(aes(time/365, S166), col="#D55E00", lwd=1) +
	geom_line(aes(time/365, S182), col="#0072B2", lwd=1) +
	geom_hline(yintercept=1/1.8, lty=2) +
	scale_x_continuous("Year", expand=c(0,0), limits=c(NA, 7+185/365),
										 breaks=1:10) +
	scale_y_continuous("Proportion susceptible", limits=c(0.5, 0.6),
										 breaks=c(0.5, 0.55, 0.6),
										 expand=c(0, 0)) +
	ggtitle("E") +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g8 <- ggplot(simulate_flu_intro_endemic_yx) +
	geom_line(data=filter(simulate_flu_intro_endemic_y, time <= 550), aes(time/365, S166), col="#D55E00", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_endemic_y, time <= 550), aes(time/365, S182), col="#0072B2", lwd=1) +
	geom_line(aes(time/365, S182), col="#0072B2", lwd=1) +
	geom_line(aes(time/365, S166), col="#D55E00", lwd=1) +
	geom_hline(yintercept=1/1.8, lty=2) +
	scale_x_continuous("Year", expand=c(0,0), limits=c(NA, 7+185/365),
										 breaks=1:10) +
	scale_y_continuous("Proportion susceptible", limits=c(0.5, 0.6),
										 breaks=c(0.5, 0.55, 0.6),
										 expand=c(0, 0)) +
	ggtitle("F") +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

gcomb <- ggarrange(g1, g2, g3,
					g6, g7, g8, nrow=2, 
					widths=c(1, 1, 1),
					draw=FALSE)

simulate_flu_intro_factorial2 <- simulate_flu_intro_factorial %>%
	mutate(
		R01=R0,
		R02=R0*R0ratio,
		s=1-cross,
		S11=1/R01,
		S21=1/R01+(1-1/R01)*s,
		S22=1/R02,
		S12=1/R02+(1-1/R02)*s,
		nichediff_new=1/sqrt(S11*S22/(S12*S21)),
		fitnessdiff_new=R02/R01*sqrt(S22*S21/(S11*S12)),
		fitnessdiff_new=pmax(fitnessdiff_new, 1/fitnessdiff_new)
	) %>%
	mutate(
		invasion=nichediff > fitnessdiff,
		invasion_new=nichediff_new > fitnessdiff_new,
		group=ifelse(invasion_new, "Mutual invasion\n(Standard models)", "Competitive exclusion"),
		group=ifelse(invasion, "Mutual invasion\n(Gog & Grenfell, 2002)", group),
		group=ifelse(xy_trough_x >= 1e-7 & invasion, "Resident persistence", group),
		group=ifelse(xy_trough_y >= 1e-7 & invasion, "Invader persistence", group),
		group=ifelse(xy_trough_x >= 1e-7 & xy_trough_y >= 1e-7 & invasion, "Coexistence", group),
		R0=paste0("Resident~R[0]==",R0),
		R0=factor(R0, 
							levels=c("Resident~R[0]==1.5", "Resident~R[0]==3", "Resident~R[0]==6", "Resident~R[0]==12")),
		group=factor(group,
								 levels=c("Competitive exclusion",
								 				 "Mutual invasion\n(Standard models)",
								 				 "Mutual invasion\n(Gog & Grenfell, 2002)",
								 				 "Resident persistence",
								 				 "Invader persistence",
								 				 "Coexistence"),
								 labels=c("Mutual invasion not possible",
								 				 "Mutual invasion\n(Standard models)",
								 				 "Mutual invasion\n(Gog & Grenfell, 2002)",
								 				 "Resident persistence",
								 				 "Invader persistence",
								 				 "Co-circulation"))
	)

g11 <- ggplot(simulate_flu_intro_factorial2) +
	geom_raster(aes(cross, R0ratio, fill=group)) +
	geom_hline(yintercept=1, lty=2, col="white") +
	scale_x_continuous("Cross immunity", expand=c(0, 0),
										 breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
	scale_y_log10(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_manual(values=c("gray90", viridisLite::viridis(5))) +
	facet_wrap(~R0, nrow=1, labeller = label_parsed) +
	ggtitle("G") +
	theme(
		strip.background = element_blank(),
		legend.position = "bottom",
		legend.title = element_blank(),
		panel.spacing = unit(5, "mm"),
		strip.text = element_text(size=12)
	)

simulate_flu_intro_factorial_endemic2 <- simulate_flu_intro_factorial_endemic %>%
	mutate(
		R01=R0,
		R02=R0*R0ratio,
		s=1-cross,
		S11=1/R01,
		S21=1/(1 + (R01-1)*(1-s)),
		S22=1/R02,
		S12=1/(1 + (R02-1)*(1-s)),
		nichediff_new=1/sqrt(S11*S22/(S12*S21)),
		fitnessdiff_new=R02/R01*sqrt(S22*S21/(S11*S12)),
		fitnessdiff_new=pmax(fitnessdiff_new, 1/fitnessdiff_new)
	) %>%
	mutate(
		invasion=nichediff > fitnessdiff,
		invasion_new=nichediff_new > fitnessdiff_new,
		group=ifelse(invasion, "Mutual invasion\n(Standard models)", "Competitive exclusion"),
		group=ifelse(invasion_new, "Mutual invasion\n(Gog & Grenfell, 2002)", group),
		group=ifelse(xy_trough_x >= 1e-7 & invasion_new, "Resident persistence", group),
		group=ifelse(xy_trough_y >= 1e-7 & invasion_new, "Invader persistence", group),
		group=ifelse(xy_trough_x >= 1e-7 & xy_trough_y >= 1e-7 & invasion_new, "Coexistence", group),
		R0=paste0("Resident~R[0]==",R0),
		R0=factor(R0, 
							levels=c("Resident~R[0]==1.5", "Resident~R[0]==3", "Resident~R[0]==6", "Resident~R[0]==12")),
		group=factor(group,
								 levels=c("Competitive exclusion",
								 				 "Mutual invasion\n(Standard models)",
								 				 "Mutual invasion\n(Gog & Grenfell, 2002)",
								 				 "Resident persistence",
								 				 "Invader persistence",
								 				 "Coexistence"),
								 labels=c("Mutual invasion not possible",
								 				 "Mutual invasion\n(Standard models)",
								 				 "Mutual invasion\n(Gog & Grenfell, 2002)",
								 				 "Resident persistence",
								 				 "Invader persistence",
								 				 "Co-circulation"))
	)

g12 <- ggplot(simulate_flu_intro_factorial_endemic2) +
	geom_raster(aes(cross, R0ratio, fill=group)) +
	geom_hline(yintercept=1, lty=2, col="white") +
	scale_x_continuous("Cross immunity", expand=c(0, 0),
										 breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
	scale_y_log10(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_manual(values=c("gray90", viridisLite::viridis(5)[1:5])) +
	facet_wrap(~R0, nrow=1, labeller = label_parsed) +
	# ggtitle("L. Resident at equilibrium, invader near equilibrium") +
	theme(
		strip.background = element_blank(),
		legend.title = element_blank(),
		panel.spacing = unit(5, "mm"),
		strip.text = element_text(size=12),
		legend.position = "bottom"
	)

gfinal <- arrangeGrob(gcomb, g11, nrow=2, heights=c(1.5, 1.2))

ggsave("figure_flu.pdf", gfinal, width=12, height=8)

ggsave("figure_flu_equi.pdf", g12, width=12, height=4)

g4 <- ggplot(simulate_flu_intro_xy) +
	geom_line(aes(time/365, I166), col="#D55E00", lwd=1) +
	# geom_line(data=simulate_flu_intro_x, aes(time/365, I166), col="#D55E00", lty=2) +
	geom_line(data=filter(simulate_flu_intro_xy, I182 > 0), aes(time/365, I182), col="#0072B2", lwd=1) +
	annotate("point", x=13.5, y=1e-6, shape=21, size=3, stroke=1.5, col="#0072B2", fill="white") +
	scale_x_continuous("Year", expand=c(0,0), breaks=13:17) +
	scale_y_log10("Prevalence") +
	ggtitle("A") +
	coord_cartesian(ylim=c(1e-8, 0.1)) +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g5 <- ggplot(simulate_flu_intro_yx) +
	geom_line(aes(time/365, I182), col="#0072B2", lwd=1) +
	geom_line(data=filter(simulate_flu_intro_yx, I166 > 0), aes(time/365, I166), col="#D55E00", lwd=1) +
	annotate("point", x=13.5, y=1e-6, shape=21, size=3, stroke=1.5, col="#D55E00", fill="white") +
	scale_x_continuous("Year", expand=c(0,0), breaks=13:17) +
	scale_y_log10("Prevalence") +
	ggtitle("B") +
	coord_cartesian(ylim=c(1e-8, 0.1)) +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g9 <- ggplot(simulate_flu_intro_xy) +
	geom_line(aes(time/365, S166), col="#D55E00", lwd=1) +
	geom_line(aes(time/365, S182), col="#0072B2", lwd=1) +
	geom_hline(yintercept=1/1.8, lty=2) +
	scale_x_continuous("Year", expand=c(0,0), breaks=13:17) +
	scale_y_continuous("Proportion susceptible", limits=c(0.3, 0.71), expand=c(0, 0)) +
	ggtitle("C") +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

g10 <- ggplot(simulate_flu_intro_yx) +
	geom_line(aes(time/365, S182), col="#0072B2", lwd=1) +
	geom_line(aes(time/365, S166), col="#D55E00", lwd=1) +
	geom_hline(yintercept=1/1.8, lty=2) +
	scale_x_continuous("Year", expand=c(0,0), breaks=13:17) +
	scale_y_continuous("Proportion susceptible", limits=c(0.3, 0.71), expand=c(0, 0)) +
	ggtitle("D") +
	theme(
		panel.grid = element_blank(),
		legend.position = "none"
	)

gcomb_supp <- ggarrange(g4, g5,
												g9, g10, 
												nrow=2, 
												draw=FALSE)

ggsave("figure_flu_nonendemic.pdf", gcomb_supp, width=6, height=4)

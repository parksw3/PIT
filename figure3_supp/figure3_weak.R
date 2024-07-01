library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)

load("simulate_rsv_factorial_weak.rda")

simulate_rsv_factorial2 <- simulate_rsv_factorial_weak %>%
	mutate(
		invasion=nichediff > fitnessdiff,
		group=ifelse(invasion, "Mutual invasion possible", "Mutual invasion impossible"),
		group=ifelse(B_trough >= 1e-7 & invasion, "Resident persistence", group),
		group=ifelse(A_trough >= 1e-7 & invasion, "Invader persistence", group),
		group=ifelse(A_trough >= 1e-7 & B_trough >= 1e-7 & invasion, "Co-circulation", group),
		group=factor(group,
								 levels=c("Mutual invasion impossible",
								 				 "Mutual invasion possible",
								 				 "Resident persistence",
								 				 "Invader persistence",
								 				 "Co-circulation"),
								 labels=c("Strain replacement & mutual invasion impossible",
								 				 "Strain replacement & extinction",
								 				 "Strain replacement & resident persistence",
								 				 "Strain replacement & invader persistence",
								 				 "Co-circulation")),
		returntime=ifelse(!invasion, NA, returntime)
	)

g1 <- ggplot(simulate_rsv_factorial2) +
	geom_raster(aes(immune_dur, alpha, fill=group)) +
	scale_x_log10("Duration of immunity (years)", expand=c(0, 0)) +
	scale_y_reverse(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_viridis_d("Return time\n(years)") +
	guides(fill=guide_legend(nrow=2)) +
	theme(
		strip.background = element_blank(),
		legend.title = element_blank(),
		legend.position = "top"
	)

g3 <- ggplot(simulate_rsv_factorial2) +
	geom_raster(aes(immune_dur, alpha, fill=returntime/52)) +
	scale_x_log10("Duration of immunity (years)", expand=c(0, 0))  +
	scale_y_reverse(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_viridis_c("Return time\n(years)", option="A",
											 breaks=c(1:5*2), labels=c(2, 4, 6, 8, ">10"),
											 limits=c(0, 10))

g2 <- ggplot(simulate_rsv_factorial2) +
	geom_raster(aes(immune_dur, alpha, fill=overcomp)) +
	scale_x_log10("Duration of immunity (years)", expand=c(0, 0))  +
	scale_y_reverse(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_viridis_c("Degree of\novercompensation", option="A")

g1a <- ggarrange(g1, labels="A")

g2a <- ggarrange(g2, g3, nrow=2, labels=c("B", "C"))

gfinal <- arrangeGrob(g1a, g2a, nrow=1)

ggsave("figure3_weak.pdf", gfinal, width=12, height=6)

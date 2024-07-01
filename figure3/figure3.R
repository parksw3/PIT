library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)

load("simulate_rsv_factorial.rda")

simulate_rsv_factorial2 <- simulate_rsv_factorial %>%
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
	geom_point(x=log10(1/0.52), y=-1/0.9159, size=10, shape=21, stroke=2) +
	scale_x_log10("Duration of immunity (years)", expand=c(0, 0)) +
	scale_y_reverse(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_viridis_d("Return time\n(years)") +
	guides(fill=guide_legend(nrow=2)) +
	theme(
		strip.background = element_blank(),
		legend.title = element_blank(),
		legend.position = "top"
	)

pointdata <- data.frame(
	immune_dur=c((1/0.52), (1/0.52), (10)),
	alpha=c(1/0.9159, 1.4, 1),
	label=c("A", "B", "C")
)

g2 <- ggplot(simulate_rsv_factorial2) +
	geom_raster(aes(immune_dur, alpha, fill=returntime/52)) +
	geom_point(data=pointdata, aes(immune_dur, alpha), col="white", shape=1, size=7, stroke=2) +
	geom_text(data=pointdata, aes(immune_dur, alpha, label=label), col="white", shape=1, size=7, stroke=2, hjust=-1, vjust=1.5) +
	scale_x_log10("Duration of immunity (years)", expand=c(0, 0))  +
	scale_y_reverse(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_viridis_c("Return time\n(years)", option="A",
											 breaks=c(1:5*2), labels=c(2, 4, 6, 8, ">10"),
											 limits=c(0, 10))

g3 <- ggplot(simulate_rsv_factorial2) +
	geom_raster(aes(immune_dur, alpha, fill=overcomp)) +
	geom_point(data=pointdata, aes(immune_dur, alpha), col="white", shape=1, size=7, stroke=2) +
	geom_text(data=pointdata, aes(immune_dur, alpha, label=label), col="white", shape=1, size=7, stroke=2, hjust=-1, vjust=1.5) +
	scale_x_log10("Duration of immunity (years)", expand=c(0, 0))  +
	scale_y_reverse(expression(Invader~R[0]/Resident~R[0]), expand=c(0, 0)) +
	scale_fill_viridis_c("Degree of\novercompensation", option="A")

ggsave("figure3a.pdf", g1, width=6, height=6)
ggsave("figure3b.pdf", g2, width=5, height=4)
ggsave("figure3c.pdf", g3, width=5, height=4)

ggsave("figure3a.png", g1, width=6, height=6)
ggsave("figure3b.png", g2, width=5, height=4)
ggsave("figure3c.png", g3, width=5, height=4)

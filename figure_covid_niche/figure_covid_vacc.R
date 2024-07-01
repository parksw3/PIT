library(tidyr)
library(dplyr)
library(ggrepel)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
load("../simulate_covid/simulate_niche_covid.rda")
load("../simulate_covid/simulate_niche_covid_vacc.rda")

data_meijiers_baseline <- bind_rows(
	data_meijers_niche1,
	data_meijers_niche2,
	data_meijers_niche3,
	data_meijers_niche4
) %>%
	mutate(
		type="naive"
	)

data_meijiers_vacc <- bind_rows(
	data_meijers_vacc_niche1,
	data_meijers_vacc_niche2,
	data_meijers_vacc_niche3,
	data_meijers_vacc_niche4,
	data_meijers_vacc_niche5,
	data_meijers_vacc_niche6,
	data_meijers_vacc_niche7,
	data_meijers_vacc_niche8
)

data_meijiers_comb <-
	bind_rows(
		data_meijiers_baseline, 
		data_meijiers_vacc
	) %>%
	mutate(
		fitnessdiff=pmax(fitnessdiff, 1/fitnessdiff),
		type=factor(type,
								levels=c("naive", "30% vac", "40% vac, 30% bst", "20% vac, 50% bst"),
								labels=c("unvaccinated", "30% vac", "40% vac, 30% bst", "20% vac, 50% bst")),
		system=gsub("SARS-CoV-2 ", "", system),
		system=factor(system,
									levels=c("WT, Alpha",
													 "Alpha, Delta",
													 "Delta, BA.1",
													 "BA.1, BA.2"))
	)

xmin <- 0
xend <- 1
xlength <- 21

data_below <- data.frame(
	x=seq(xmin, xend, length.out=xlength),
	y=1-seq(xmin, xend, length.out=xlength)
)

data_above <- data.frame(
	x=seq(xmin, xend, length.out=xlength),
	y=1/(1-seq(xmin, xend, length.out=xlength))
)

g1 <- ggplot(data_meijiers_comb) +
	geom_hline(yintercept=1, lty=2, col="gray") +
	geom_ribbon(data=data_below, ymin=-Inf,  aes(x=x,ymax=y), fill="gray40", alpha=0.3) +
	geom_ribbon( data=data_above, ymax=Inf,  aes(x=x,ymin=y), fill="gray40", alpha=0.3) +
	geom_function(fun=function(x) 1-x) +
	geom_function(fun=function(x) 1/(1-x)) +
	geom_point(aes(1-1/nichediff, pmax(fitnessdiff, fitnessdiff), col=system, shape=type), size=2, stroke=1) +
	annotate("text", x=0.69, y=2.7, label="Mutual invasion\npossible", family='Times', hjust=0) +
	annotate("text", x=1-0.5, y=2.7, label="Mutual invasion\nnot possible", family='Times', hjust=1) +
	# annotate("segment", x=0.31, y=1.12, xend=0.325, yend=1.17, arrow=arrow(length=unit(0.04, "inches")), col="#35B779FF", alpha=0.8) +
	# annotate("segment", x=0.36, y=1.21, xend=0.45, yend=1.24, arrow=arrow(length=unit(0.04, "inches")), col="#35B779FF", alpha=0.8) +
	# annotate("segment", x=0.485, y=1.23, xend=0.505, yend=1.18, arrow=arrow(length=unit(0.04, "inches")), col="#35B779FF", alpha=0.8) +
	# geom_label_repel(data=data_meijiers_baseline, aes(nichediff, fitnessdiff, label=system), family="Times", size=3, box.padding = 1,
	# 								 max.overlaps = 20,
	# 								 fill="white") +
	scale_x_continuous("Immunological niche difference", limits=c(0, 1), expand=c(0, 0),
										 breaks=0:5*0.2) +
	scale_y_log10("Fitness difference", expand=c(0, 0), 
								breaks=c(0.5, 1, 2)) +
	scale_shape_manual("Vaccination", values=1:4) +
	scale_color_viridis_d("Variants") +
	coord_cartesian(xlim=c(0, 1), ylim=c(0.5, 3)) +
	guides(color=guide_legend(ncol=1, title.position = "top"), shape=guide_legend(ncol=1, title.position = "top")) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		legend.position = "right",
		legend.box="horizontal"
	)

ggsave("figure_covid_vacc.pdf", g1, width=7, height=4.5)

library(tidyr)
library(dplyr)
library(ggrepel)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
load("../simulate_covid/simulate_niche_covid_omicron_sensitivity.rda")

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

simulate_niche_covid_omicron_sensitivity2 <- simulate_niche_covid_omicron_sensitivity %>%
	mutate(
		cross=paste0(cross*100, "% cross immunity"),
		type=factor(type,
								levels=c("unvaccinated", "30% vac", "40% vac, 30% bst", "20% vac, 50% bst"),
								labels=c("unvaccinated", "30% vac", "40% vac, 30% bst", "20% vac, 50% bst")),
		Rratio=paste0(Rratio, "x transmission advantage")
	)

g1 <- ggplot(simulate_niche_covid_omicron_sensitivity2) +
	geom_hline(yintercept=1, lty=2, col="gray") +
	geom_ribbon( data=data_above, ymax=Inf,  aes(x=x,ymin=y), fill="gray40", alpha=0.3) +
	geom_function(fun=function(x) 1/(1-x)) +
	geom_point(aes(1-1/nichediff, fitnessdiff, col=factor(R0delta), shape=type), size=2, stroke=1) +
	scale_x_continuous("Immunological niche difference", limits=c(0, 1), expand=c(0, 0),
										 breaks=0:5*0.2) +
	scale_y_log10("Fitness difference", expand=c(0, 0),
								breaks=c(1, 2)) +
	scale_color_viridis_d(expression("Delta"~R[0])) +
	scale_shape_manual("Scenario", values=1:4) +
	coord_cartesian(xlim=c(0, 1.05), ylim=c(1, 3)) +
	facet_grid(Rratio~cross) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		strip.background = element_blank()
	)

ggsave("figure_niche_covid_omicron_sensitivity.pdf", g1, width=8, height=6)

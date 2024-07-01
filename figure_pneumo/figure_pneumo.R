library(tidytext)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(gridExtra)
load("../simulate_niche/simulate_niche_cobey.rda")

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

data_cobey_niche_comb <- bind_rows(
	data_cobey_niche_base %>% mutate(type="Baseline"),
	data_cobey_niche_anti %>% mutate(type="Serotype-specific immunity"),
	data_cobey_niche_nons %>% mutate(type="Non-specific immunity"),
	data_cobey_niche_both %>% mutate(type="Serotype-specific+non-specific immunity")
) %>%
	mutate(
		type=factor(type, levels=c("Baseline",
															 "Serotype-specific immunity",
															 "Non-specific immunity",
															 "Serotype-specific+non-specific immunity")),
		coexist=1/nichediff < fitnessdiff & fitnessdiff<nichediff
	)

dd1 <- data_cobey_niche_comb %>%
	filter(type=="Baseline", system==c("1, 25")) %>%
	mutate(
		xend=0.5, yend=3.5
	)

dd2 <- data_cobey_niche_comb %>%
	filter(type=="Baseline", system==c("13, 25")) %>%
	mutate(
		xend=0.6, yend=2.5
	)

dd3 <- data_cobey_niche_comb %>%
	filter(type=="Baseline", system==c("24, 25")) %>%
	mutate(
		xend=0.9, yend=1.8
	)

dd4 <- data_cobey_niche_comb %>%
	filter(type=="Serotype-specific+non-specific immunity", system==c("1, 25")) %>%
	mutate(
		xend=0.5, yend=3
	)

dd5 <- data_cobey_niche_comb %>%
	filter(type=="Serotype-specific+non-specific immunity", system==c("13, 25")) %>%
	mutate(
		xend=0.6, yend=2.3
	)

g3 <- ggplot(data_cobey_niche_comb) +
	geom_hline(yintercept=1, lty=2, col="gray") +
	geom_ribbon(data=data_above, ymax=Inf,  aes(x=x,ymin=y), fill="gray40", alpha=0.3) +
	geom_ribbon(data=data_below, ymin=-Inf,  aes(x=x,ymax=y), fill="gray40", alpha=0.3) +
	geom_function(fun=function(x) 1-x) +
	geom_function(fun=function(x) 1/(1-x)) +
	geom_segment(data=dd1, aes(1-1/nichediff, fitnessdiff, xend=1-xend, yend=yend), col="orange") +
	geom_text(data=dd1, aes(1-xend, yend), col="orange", label=c("Serotypes 1, 25"), hjust=1, vjust=0, family="Times") +
	geom_segment(data=dd2, aes(1-1/nichediff, fitnessdiff, xend=1-xend, yend=yend), col="orange") +
	geom_text(data=dd2, aes(1-xend, yend), col="orange", label=c("Serotypes 13, 25"), hjust=1, vjust=-0.5, family="Times") +
	geom_segment(data=dd3, aes(1-1/nichediff, fitnessdiff, xend=1-xend, yend=yend), col="black") +
	geom_text(data=dd3, aes(1-xend, yend), col="black", label=c("Serotypes 24, 25"), hjust=0.2, vjust=-0.5, family="Times") +
	geom_segment(data=dd4, aes(1-1/nichediff, fitnessdiff, xend=1-xend, yend=yend), col="black") +
	geom_text(data=dd4, aes(1-xend, yend), col="black", label=c("Serotypes 1, 25"), hjust=1, vjust=-0.5, family="Times") +
	geom_segment(data=dd5, aes(1-1/nichediff, fitnessdiff, xend=1-xend, yend=yend), col="black") +
	geom_text(data=dd5, aes(1-xend, yend), col="black", label=c("Serotypes 13, 25"), hjust=1, vjust=-0.5, family="Times") +
	geom_point(aes(1-1/nichediff, fitnessdiff, col=coexist, shape=coexist), size=1) +
	scale_x_continuous("Immunological niche difference", limits=c(0, 1), expand=c(0, 0),
										 breaks=0:5*0.2) +
	scale_y_log10("Fitness difference", expand=c(0, 0)) +
	scale_color_manual(values=c("orange", "black")) +
	scale_shape_manual(values=c(2, 1)) +
	coord_cartesian(xlim=c(0,1), ylim=c(0.5, 4)) +
	facet_wrap(~type, nrow=2, scale="free") +
	ggtitle("C") +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		strip.background = element_blank(),
		legend.position = "none"
	)

niche_data <- data.frame(
	effect=-c(1/data_cobey_niche_anti$nichediff-1/data_cobey_niche_base$nichediff,
					 1/data_cobey_niche_nons$nichediff-1/data_cobey_niche_base$nichediff,
					 1/data_cobey_niche_both$nichediff-1/data_cobey_niche_base$nichediff-(1/data_cobey_niche_anti$nichediff-1/data_cobey_niche_base$nichediff)-(1/data_cobey_niche_nons$nichediff-1/data_cobey_niche_base$nichediff)),
	group=c(rep("Serotype-specific immunity", nrow(data_cobey_niche_anti)),
					rep("Non-specific immunity", nrow(data_cobey_niche_anti)),
					rep("Interaction effect", nrow(data_cobey_niche_anti)))
)

fitness_data <- data.frame(
	effect=c(data_cobey_niche_anti$fitnessdiff-data_cobey_niche_base$fitnessdiff,
					 data_cobey_niche_nons$fitnessdiff-data_cobey_niche_base$fitnessdiff,
					 data_cobey_niche_both$fitnessdiff-data_cobey_niche_base$fitnessdiff-(data_cobey_niche_anti$fitnessdiff-data_cobey_niche_base$fitnessdiff)-(data_cobey_niche_nons$fitnessdiff-data_cobey_niche_base$fitnessdiff)),
	group=c(rep("Serotype-specific immunity", nrow(data_cobey_niche_anti)),
					rep("Non-specific immunity", nrow(data_cobey_niche_anti)),
					rep("Interaction effect", nrow(data_cobey_niche_anti)))
)

g1 <- ggplot(niche_data) +
	geom_vline(xintercept=0, lty=2) +
	geom_boxplot(aes(effect, group)) +
	scale_x_continuous("Changes in niche difference", limits=c(-0.62, 0.62),
										 breaks=(-3):3*2/10) +
	scale_y_discrete("") +
	ggtitle("A") +
	theme(
		panel.border = element_rect(linewidth=1)
	)

g2 <- ggplot(fitness_data) +
	geom_vline(xintercept=0, lty=2) +
	geom_boxplot(aes(effect, group)) +
	scale_x_continuous("Changes in fitness difference", limits=c(-2.3, 2.3)) +
	scale_y_discrete("") +
	ggtitle("B") +
	theme(
		panel.border = element_rect(linewidth=1)
	)

gcomb1 <- ggarrange(g1, g2, nrow=2)

gfinal <- arrangeGrob(gcomb1, g3, nrow=1, widths=c(1, 1.5))

ggsave("figure_pneumo.pdf", gfinal, width=10, height=6)

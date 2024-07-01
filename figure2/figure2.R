library(tidyr)
library(dplyr)
library(ggrepel)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
load("../simulate_niche/simulate_niche.rda")
load("../simulate_niche/simulate_niche_pox.rda")
load("../simulate_niche/simulate_niche_dengue.rda")
load("../simulate_covid/simulate_niche_covid.rda")
load("../smith_niche/smith_niche.rda")

data_all_niche <- bind_rows(
	data_kissler_niche %>% mutate(interaction="HCoV HKU1, OC43"),
	data_white_niche %>% mutate(interaction="RSV A, B"),
	data_rohani_niche %>% mutate(interaction="Prevaccination childhood infections"),
	data_rohani_niche2 %>% mutate(interaction="Prevaccination childhood infections"),
	data_rohani_niche3 %>% mutate(interaction="Prevaccination childhood infections"),
	data_koelle_niche %>% mutate(interaction="Cholera serotypes"),
	data_yang_niche %>% mutate(interaction="Seasonal influenza (between subtype)"),
	data_yang_niche2 %>% mutate(interaction="Seasonal influenza (between subtype)"),
	data_yang_niche3 %>% mutate(interaction="Seasonal influenza (between subtype)"),
	data_bhattacharyya_niche1 %>% mutate(interaction="Paramyxoviruses family"),
	data_bhattacharyya_niche2 %>% mutate(interaction="Paramyxoviruses family"),
	data_bhattacharyya_niche3 %>% mutate(interaction="Paramyxoviruses family"),
	data_bhattacharyya_niche4 %>% mutate(interaction="Paramyxoviruses family"),
	data_meijers_niche1 %>% mutate(interaction="SARS-CoV-2 variants"),
	data_meijers_niche2 %>% mutate(interaction="SARS-CoV-2 variants"),
	data_meijers_niche3 %>% mutate(interaction="SARS-CoV-2 variants"),
	data_meijers_niche4 %>% mutate(interaction="SARS-CoV-2 variants"),
	simulate_niche_dengue %>% mutate(interaction="Dengue serotypes"),
	simulate_niche_pox %>% summarize(nichediff=mean(nichediff), fitnessdiff=mean(fitnessdiff)) %>% mutate(interaction = "Smallpox, mpox"),
	smith_niche %>% mutate(interaction = "Seasonal influenza (within subtype)")
) %>%
	mutate(
		fitnessdiff=pmax(fitnessdiff, 1/fitnessdiff),
		interaction=factor(interaction,
											 levels=c("SARS-CoV-2 variants",
											 				 "Seasonal influenza (within subtype)",
											 				 "Cholera serotypes",
											 				 "Dengue serotypes",
											 				 "HCoV HKU1, OC43",
											 				 "Paramyxoviruses family",
											 				 "Prevaccination childhood infections",
											 				 "RSV A, B",
											 				 "Seasonal influenza (between subtype)",
											 				 "Smallpox, mpox"))
	)

simulate_niche_pox2 <- simulate_niche_pox %>%
	mutate(
		nicchange=100*(log(nichediff)-log(meanR0))/log(meanR0),
		fitchange=100*(fitnessdiff-nulldiff)/nulldiff
	) %>%
	summarize(
		system=unique(system),
		nicchange_lwr=min(nicchange),
		nicchange_upr=max(nicchange),
		fitchange_lwr=min(fitchange),
		fitchange_upr=max(fitchange),
		nicchange=mean(nicchange),
		fitchange=mean(fitchange),
		interaction="Smallpox, mpox"
	)

data_all_niche2 <- data_all_niche %>%
	mutate(
		nicchange=100*(log(nichediff)-log(meanR0))/log(meanR0),
		fitchange=100*(fitnessdiff-nulldiff)/nulldiff,
		nicchange_lwr=NA,
		nicchange_upr=NA,
		fitchange_lwr=NA,
		fitchange_upr=NA
	) %>%
	bind_rows(simulate_niche_pox2)

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

# simulate_niche_pox_hull <- simulate_niche_pox[chull(simulate_niche_pox$nichediff, simulate_niche_pox$fitnessdiff),]

## jitter fitnessdiff to draw ellipse
set.seed(101)
ellipsedata1 <- data_all_niche %>% 
	filter(interaction=="Seasonal influenza (within subtype)") %>%
	mutate(
		fitnessdiff=c(rep(1, 4), 1.12, 0.88, rep(1, 4)),
		nichediff=1/c(0.66, 0.66, 0.66, 0.66, 0.75, 0.75, 0.9, 0.9, 0.9, 0.9)
	)

ellipsedata2 <- data_all_niche %>% 
	filter(interaction=="SARS-CoV-2 variants") 

ellipsedata3 <- bind_rows(
	ellipsedata2,
	ellipsedata2,
	ellipsedata2,
	ellipsedata2,
	ellipsedata2
)

varcol <- c("#E02938", viridis::viridis(9)[1])

g1 <- ggplot(data_all_niche) +
	geom_hline(yintercept=1, lty=2, col="gray") +
	geom_ribbon(data=data_below, ymin=-Inf,  aes(x=x,ymax=y), fill="gray40", alpha=0.3) +
	geom_ribbon(data=data_above, ymax=Inf,  aes(x=x,ymin=y), fill="gray40", alpha=0.3) +
	geom_function(fun=function(x) 1-x) +
	geom_function(fun=function(x) 1/(1-x)) +
	# geom_polygon(data=simulate_niche_pox_hull, aes(1/nichediff, fitnessdiff), fill="darkred", col="darkred", alpha=0.2) +
	# annotate("label", x=1.7, y=3.5, label="Smallpox, mpox", family="Times", col="darkred") +
	annotate("text", x=0.69, y=2.7, label="Mutual invasion\npossible", family='Times', hjust=0) +
	annotate("text", x=1-0.5, y=2.7, label="Mutual invasion\nnot possible", family='Times', hjust=1) +
	geom_point(aes(1-1/nichediff, fitnessdiff,
								 col=interaction, shape=interaction, size=interaction), stroke=1) +
	stat_ellipse(data=ellipsedata1, aes(1-1/nichediff, fitnessdiff), lty=2, lwd=0.7, col=varcol[2]) +
	stat_ellipse(data=ellipsedata3, aes(1-1/nichediff, fitnessdiff), lty=2, lwd=0.7, col=varcol[1]) +
	annotate("text", x=1-0.35, y=1.34, label="SARS-CoV-2 variants", col=varcol[1], family="Times") +
	annotate("text", x=0.25, y=0.83, label="Seasonal influenza\n(within subtype)", col=varcol[2], family="Times", hjust=0) +
	# geom_label_repel(aes(nichediff, fitnessdiff, label=system), family="Times", size=3, box.padding = 0.9,
	#  								 min.segment.length = 0.3,
	#  								 max.overlaps = 25,
	#  								 fill="white") +
	scale_x_continuous("Immunological niche difference", expand=c(0, 0),
										 breaks=0:5*0.2) +
	scale_y_log10("Fitness difference", expand=c(0, 0), 
								breaks=c(0.5, 1, 2)) +
	scale_shape_manual(values=1:10) +
	scale_size_manual(values=c(4, 4, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5)) +
	scale_color_manual(values=c("#E02938", viridis::viridis(9))) +
	coord_cartesian(ylim=c(0.5, 3), xlim=c(0, 1)) +
	guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=1)) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		legend.position = "right",
		legend.direction = "vertical",
		legend.title = element_blank()
	)

ggsave("figure2.pdf", g1, width=7, height=4.5)

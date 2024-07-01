library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(ggtree)

tree <- read.tree("../data/nextstrain_flu_seasonal_h3n2_ha_12y_timetree.nwk")

g_flu_phylo <- ggplot(tree) + 
	geom_tree() +
	scale_x_continuous("Year", breaks=c(0, 5, 10, 15),
										 labels=c("'05", "'10", "'15", "'20")) +
	ggtitle("E. Strain replacement, Seasonal influenza") +
	theme(
		panel.border = element_blank(),
		panel.grid = element_blank(),
		axis.line.y = element_blank(),
		axis.title = element_blank(),
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank()
	)

save("g_flu_phylo", file="figure1_flu_phylo.rda")

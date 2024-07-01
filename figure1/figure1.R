library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(ggtree)
library(gridExtra)
load("figure1_rsv_korea.rda")
load("figure1_dengue_thailand.rda")
load("figure1_hfmd.rda")
load("figure1_covid.rda")
load("figure1_flu_phylo.rda")

gfinal <- arrangeGrob(g_rsv_korea, g_hfmd, g_dengue_thailand,
											g_covid, g_flu_phylo, nrow=2)

ggsave("figure1.pdf", gfinal, width=14, height=8)

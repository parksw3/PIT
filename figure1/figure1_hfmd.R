library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(vroom)

covid_color <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

data_hfmd <- vroom("../data/saki_hfmd.csv")

data_hfmd_raw <- data_hfmd %>%
	group_by(PROV_ID) %>%
	arrange(PROV_ID, YEAR, WEEK) %>%
	mutate(
		time=2009+1:n()/52
	) %>%
	group_by(time) %>%
	summarize(
		EVA71=sum(EVA71),
		CVA16=sum(CVA16)
	) %>%
	gather(key, value, -time)

data_hfmd_raw2 <- data_hfmd_raw %>%
	group_by(time) %>%
	mutate(
		value=value/sum(value)
	)
 
g1 <- ggplot(data_hfmd_raw) +
	geom_bar(aes(time, value, fill=key), stat="identity") +
	scale_x_continuous(expand=c(0, 0),
										 limits=c(2009, 2014)) +
	scale_y_continuous("Cases", expand=c(0, 0), limits=c(0, 70000)) +
	scale_fill_manual(values=covid_color[1:6]) +
	guides(fill=guide_legend(nrow=1)) +
	ggtitle("B. Stable coexistence, Enterovirus") +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		# axis.line = element_line(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank(),
		legend.position = "none",
		legend.title = element_blank()
	)

g2 <- ggplot(data_hfmd_raw2) +
	geom_bar(aes(time, value, fill=key), stat="identity") +
	scale_x_continuous(expand=c(0, 0),
										 limits=c(2009, 2014),
										 breaks=c(2009, 2010, 2011, 2012, 2013, 2014),
										 labels=c("'09", "'10", "'11", "'12", "'13", "'14")) +
	scale_y_continuous("Proportion", expand=c(0, 0), limits=c(0, 1)) +
	scale_fill_manual(values=covid_color[1:6]) +
	guides(fill=guide_legend(nrow=1)) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		# axis.line = element_line(),
		axis.title.x = element_blank(),
		legend.position = "bottom",
		legend.title = element_blank()
	)

g_hfmd <- ggarrange(g1, g2, nrow=2, heights=c(1, 1), draw=FALSE)

ggsave("figure1_hfmd.pdf", g_hfmd, width=8, height=4)
save("g_hfmd", file="figure1_hfmd.rda")

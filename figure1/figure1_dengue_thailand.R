library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(vroom)
library(lubridate)
library(readxl)
# https://www.science.org/doi/10.1126/science.abk0058?intcmp=trendmd-sci

covid_color <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

data_dengue <- read_xlsx("../data/dengue_thailand.xlsx")

data_dengue_raw <- data_dengue %>%
	gather(key, value, -year)

data_dengue_prop <- data_dengue_raw %>%
	group_by(year) %>%
	mutate(
		value=value/sum(value)
	)

g1 <- ggplot(data_dengue_raw) +
	geom_bar(aes(year, value, fill=key), stat="identity") +
	scale_x_continuous(limits=c(1993.5, 2014.5), expand=c(0, 0)) +
	scale_y_continuous("Cases", expand=c(0, 0), limits=c(0, 1400)) +
	scale_fill_manual(values=covid_color[1:6]) +
	ggtitle("C. Stable coexistence, Dengue") +
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

g2 <- ggplot(data_dengue_prop) +
	geom_bar(aes(year, value, fill=key), stat="identity") +
	scale_x_continuous(limits=c(1993.5, 2014.5), expand=c(0, 0),
										 breaks=c(1995, 2000, 2005, 2010),
										 labels=c("'95", "'00", "'05", "'10")) +
	scale_y_continuous("Proportion", expand=c(0, 0), limits=c(0, 1)) +
	scale_fill_manual(values=covid_color[1:6]) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		# axis.line = element_line(),
		axis.title.x = element_blank(),
		legend.position = "bottom",
		legend.title = element_blank()
	)

g_dengue_thailand <- ggarrange(g1, g2, nrow=2, heights=c(1, 1))

ggsave("figure1_dengue_thailand.pdf", g_dengue_thailand, width=8, height=4)
save("g_dengue_thailand", file="figure1_dengue_thailand.rda")

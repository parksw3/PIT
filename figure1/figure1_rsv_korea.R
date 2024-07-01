library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(vroom)

covid_color <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

data_rsv <- vroom("../data/RSV_Korea.csv") 

data_rsv_gather <- data_rsv %>%
	gather(key, value, -year, -month) %>%
	group_by(year, month) %>%
	mutate(
		key=paste0("RSV ", key)
	)

data_rsv_gather2 <- data_rsv_gather %>%
	mutate(value=value/sum(value))
	
g1 <- ggplot(data_rsv_gather) +
	geom_bar(aes(year+(month-1)/12, value, fill=key), lwd=0.8, stat="identity") +
	scale_x_continuous(expand=c(0, 0),
										 breaks=c(2012, 2014, 2016, 2018),
										 labels=c("Jan. '12", "Jan. '14", "Jan. '16", "Jan. '18")) +
	scale_y_continuous("Detection rate (%)", expand=c(0, 0), limits=c(0, 59)) +
	scale_fill_manual(values=covid_color[1:6]) +
	guides(fill=guide_legend(nrow=1)) +
	ggtitle("A. Stable coexistence, RSV") +
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

g2 <- ggplot(data_rsv_gather2) +
	geom_bar(aes(year+(month-1)/12, value, fill=key), lwd=0.8, stat="identity") +
	scale_x_continuous(expand=c(0, 0),
										 breaks=c(2012, 2014, 2016, 2018),
										 labels=c("Jan. '12", "Jan. '14", "Jan. '16", "Jan. '18")) +
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

g_rsv_korea <- ggarrange(g1, g2, nrow=2, heights=c(1, 1), draw=FALSE)

ggsave("figure1_rsv_korea.pdf", g_rsv_korea, width=8, height=4)
save("g_rsv_korea", file="figure1_rsv_korea.rda")

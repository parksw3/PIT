library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(vroom)

data_covid <- vroom("../data/R_average.txt") 

covid_color <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

data_covid_raw <- data_covid %>%
	mutate(
		date=time+as.Date("1900-01-01"),
		WT=x_wt*tot_cases_rate,
		Alpha=x_ALPHA*tot_cases_rate,
		Delta=x_DELTA*tot_cases_rate,
		BA.1=x_BA.1*tot_cases_rate,
		BA.2=x_BA.2*tot_cases_rate,
		BA.45=(x_BA.4+x_BA.5)*tot_cases_rate
	) %>%
	select(date, WT, Alpha, Delta, BA.1, BA.2, BA.45) %>%
	gather(key, value, -date) %>%
	mutate(
		key=factor(key, levels=c("WT", "Alpha", "Delta", "BA.1", "BA.2", "BA.45"),
							 labels=c("WT", "Alpha", "Delta", "BA.1", "BA.2", "BA.4/5"))
	) %>%
	filter(
		value > 1e-6
	)

data_covid_raw2 <- data_covid_raw %>%
	mutate(
		year=year(date),
		month=month(date)
	) %>%
	group_by(year, month, key) %>%
	summarize(
		value=sum(value)
	) %>%
	filter(!(year + month/12 <= 2020+9/12 & key != "WT"))

data_covid_raw3 <- data_covid_raw2 %>%
	mutate(
		value=value/sum(value)
	)

g1 <- ggplot(data_covid_raw2) +
	geom_bar(aes(year+(month-1)/12, value/7, fill=key), lwd=0.8, stat="identity") +
	scale_x_continuous(expand=c(0, 0),
										 limits=c(2020+1.5/12, 2022+7.5/12),
										 breaks=c(2020+6/12, 2021, 2021+6/12, 2022, 2022+6/12),
										 labels=c("Jul. '20", "Jan. '21", "Jul. '21", "Jan. '22", "Jul. '22")) +
	scale_y_continuous("Cases per capita", expand=c(0, 0), limits=c(0, 0.04)) +
	scale_fill_manual(values=covid_color[1:6]) +
	guides(fill=guide_legend(nrow=1)) +
	ggtitle("D. Strain replacement, SARS-CoV-2") +
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

g2 <- ggplot(data_covid_raw3) +
	geom_bar(aes(year+(month-1)/12, value, fill=key), lwd=0.8, stat="identity") +
	scale_x_continuous(expand=c(0, 0),
										 limits=c(2020+1.5/12, 2022+7.5/12),
										 breaks=c(2020+6/12, 2021, 2021+6/12, 2022, 2022+6/12),
										 labels=c("Jul. '20", "Jan. '21", "Jul. '21", "Jan. '22", "Jul. '22")) +
	scale_y_continuous("Proportion", expand=c(0, 0), limits=c(0, NA)) +
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

g_covid <- ggarrange(g1, g2, nrow=2, heights=c(1, 1), draw=FALSE)

ggsave("figure1_covid.pdf", g_covid, width=8, height=4)
save("g_covid", file="figure1_covid.rda")

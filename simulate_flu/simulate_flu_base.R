library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../R/model_gog.R")

base_sim <- simulate_gog_stochastic(seed=913)

base_prev <- base_sim[,c(1, which(grepl("I", colnames(base_sim))))]

base_prev_gather <- base_prev %>%
	gather(key, value, -time) %>%
	mutate(
		key=as.numeric(gsub("I", "", key))
	)

base_prev_gather2 <- base_prev_gather %>% 
	filter(value > 0) %>% 
	group_by(key) %>% 
	summarize(
		key=unique(key)
	)

antigenic_cut <- c(-1, 1.5, 4, 8, 14.5, 17.5, 24, 30)

base_prev_summ <- base_prev_gather %>%
	mutate(
		year=floor((time+180)/365)
	) %>%
	group_by(year, key) %>%
	summarize(
		total=sum(value)
	) %>%
	mutate(
		antigenic=key/10,
		group=cut(antigenic, antigenic_cut)
	)

save("base_sim",
		 "base_prev", 
		 "base_prev_gather",
		 "antigenic_cut",
		 "base_prev_summ",
		 file="simulate_flu_base.rda")

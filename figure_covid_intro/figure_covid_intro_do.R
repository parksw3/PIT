library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)

load("../simulate_covid/simulate_covid_intro_endemic.rda")
load("../simulate_covid/simulate_covid_intro_endemic_vacc.rda")
load("../simulate_covid/simulate_covid_intro_endemic_do.rda")

summfun <- function(x, s1="WT", s2="Alpha", vacc="naive") {
	x %>%
		as.data.frame %>%
		mutate(
			time=1:n()
		) %>%
		rename(
			resident=imat.1,
			invader=imat.2
		) %>%
		mutate(
			variant=paste0(s1, ", ", s2),
			vacc=vacc
		)
}

alldata <- list(
	summfun(simulate_covid_delta_BA1_vacc, s1="Delta", s2="BA.1", vacc="30% vac") %>% mutate(type="Original parameterization"),
	summfun(simulate_covid_delta_BA1_boos, s1="Delta", s2="BA.1", vacc="40% vac, 30% bst") %>% mutate(type="Original parameterization"),
	summfun(simulate_covid_delta_BA1_more, s1="Delta", s2="BA.1", vacc="20% vac, 50% bst") %>% mutate(type="Original parameterization"),
	summfun(simulate_covid_delta_BA1_vacc2, s1="Delta", s2="BA.1", vacc="30% vac") %>% mutate(type="Increased transmissibility"),
	summfun(simulate_covid_delta_BA1_boos2, s1="Delta", s2="BA.1", vacc="40% vac, 30% bst") %>% mutate(type="Increased transmissibility"),
	summfun(simulate_covid_delta_BA1_more2, s1="Delta", s2="BA.1", vacc="20% vac, 50% bst") %>% mutate(type="Increased transmissibility"),
	summfun(simulate_covid_delta_BA1_vacc_sens, s1="Delta", s2="BA.1", vacc="30% vac") %>% mutate(type="Increased transmissibility\nIncreased cross immunity"),
	summfun(simulate_covid_delta_BA1_boos_sens, s1="Delta", s2="BA.1", vacc="40% vac, 30% bst") %>% mutate(type="Increased transmissibility\nIncreased cross immunity"),
	summfun(simulate_covid_delta_BA1_more_sens, s1="Delta", s2="BA.1", vacc="20% vac, 50% bst") %>% mutate(type="Increased transmissibility\nIncreased cross immunity")
) %>%
	bind_rows %>%
	mutate(
		variant=factor(variant, levels=c("WT, Alpha", "Alpha, Delta", "Delta, BA.1", "BA.1, BA.2")),
		vacc=factor(vacc, levels=c("naive", "30% vac", "40% vac, 30% bst", "20% vac, 50% bst")),
		type=factor(type, levels=c("Original parameterization", "Increased transmissibility", "Increased transmissibility\nIncreased cross immunity"))
	)

g1 <- ggplot(alldata) +
	geom_line(aes(time/365-50, resident, col="Resident"), lwd=1) +
	geom_line(aes(time/365-50, invader, col="Invader"), lwd=1) +
	scale_x_continuous("Time since invasion (years)", limits=c(49.5, 52)-50, breaks=0:2, expand=c(0, 0)) +
	scale_y_log10("Prevalence", expand=c(0, 0)) +
	coord_cartesian(ylim=c(1e-7, 1)) +
	facet_grid(type~vacc) +
	scale_color_manual(values=c("orange", "black")) +
	theme(
		panel.grid = element_blank(),
		strip.background = element_blank(),
		legend.position = "top",
		legend.title = element_blank()
	)

ggsave("figure_covid_intro_do.pdf", g1, width=8, height=6)

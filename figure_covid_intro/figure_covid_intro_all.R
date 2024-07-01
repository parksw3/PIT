library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)

load("../simulate_covid/simulate_covid_intro_endemic.rda")
load("../simulate_covid/simulate_covid_intro_endemic_vacc.rda")

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
	summfun(simulate_covid_wt_alpha),
	summfun(simulate_covid_wt_alpha_vacc, vacc="30% vac"),
	summfun(simulate_covid_alpha_delta, s1="Alpha", s2="Delta", vacc="naive"),
	summfun(simulate_covid_alpha_delta_vacc, s1="Alpha", s2="Delta", vacc="30% vac"),
	summfun(simulate_covid_alpha_delta_boos, s1="Alpha", s2="Delta", vacc="40% vac, 30% bst"),
	summfun(simulate_covid_delta_BA1_vacc, s1="Delta", s2="BA.1", vacc="30% vac"),
	summfun(simulate_covid_delta_BA1_boos, s1="Delta", s2="BA.1", vacc="40% vac, 30% bst"),
	summfun(simulate_covid_delta_BA1_more, s1="Delta", s2="BA.1", vacc="20% vac, 50% bst"),
	summfun(simulate_covid_BA1_BA2_boos, s1="BA.1", s2="BA.2", vacc="40% vac, 30% bst"),
	summfun(simulate_covid_BA1_BA2_more, s1="BA.1", s2="BA.2", vacc="20% vac, 50% bst")
) %>%
	bind_rows %>%
	mutate(
		variant=factor(variant, levels=c("WT, Alpha", "Alpha, Delta", "Delta, BA.1", "BA.1, BA.2")),
		vacc=factor(vacc, levels=c("naive", "30% vac", "40% vac, 30% bst", "20% vac, 50% bst"))
	)

g1 <- ggplot(alldata) +
	geom_line(aes(time/365-50, resident, col="Resident"), lwd=1) +
	geom_line(aes(time/365-50, invader, col="Invader"), lwd=1) +
	scale_x_continuous("Time since invasion (years)", limits=c(49.5, 52)-50, breaks=0:5, expand=c(0, 0)) +
	scale_y_log10("Incidence", expand=c(0, 0)) +
	coord_cartesian(ylim=c(1e-7, 1)) +
	facet_grid(variant~vacc) +
	scale_color_manual(values=c("orange", "black")) +
	theme(
		panel.grid = element_blank(),
		strip.background = element_blank(),
		legend.position = "top",
		legend.title = element_blank()
	)

ggsave("figure_covid_intro_all.pdf", g1, width=8, height=6)

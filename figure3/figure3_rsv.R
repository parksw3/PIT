library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../R/model_white.R")

tmax <- 25

yini_A <- c(S=1-1e-7, Ip1=1e-7, Ip2=0,
						R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)

ss_A <- simulate_white(yini=yini_A, tmax=tmax)

yini_AB <- unlist(tail(ss_A,1)[2:13])
yini_AB["Ip2"] <- 1e-7

ss_AB <- simulate_white(yini=yini_AB, tmax=20)

yini_B <- c(S=1-1e-7, Ip1=0, Ip2=1e-7,
						R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)

ss_B <- simulate_white(yini=yini_B, tmax=tmax)

yini_BA <- unlist(tail(ss_B,1)[2:13])
yini_BA["Ip1"] <- 1e-7

ss_BA <- simulate_white(yini=yini_BA, tmax=20)

g1 <- ggplot(ss_AB) +
	geom_line(data=ss_A, aes(time/52-tmax, prevalence1, col="RSV A"), lwd=1) +
	geom_line(aes(time/52, prevalence1, col="RSV A"), lwd=1) +
	geom_line(aes(time/52, prevalence2, col="RSV B"), lwd=1) +
	annotate("point", x=0, y=1e-7, shape=21, size=3, stroke=1.5, col="#0072B2", fill="white") +
	scale_x_continuous("Time since introduction (years)", expand=c(0, 0)) +
	scale_y_log10("Prevalence", breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
	scale_color_manual(values=c("#D55E00", "#0072B2")) +
	coord_cartesian(xlim=c(-1, 5), ylim=c(5e-8, 1)) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		# axis.line = element_line(),
		legend.position = c(0.4, 0.92),
		legend.direction = "horizontal",
		legend.background = element_rect(fill=NA),
		legend.title = element_blank()
	)

g2 <- ggplot(ss_BA) +
	geom_line(data=ss_B, aes(time/52-tmax, prevalence2, col="RSV B"), lwd=1) +
	geom_line(aes(time/52, prevalence1, col="RSV A"), lwd=1) +
	geom_line(aes(time/52, prevalence2, col="RSV B"), lwd=1) +
	annotate("point", x=0, y=1e-7, shape=21, size=3, stroke=1.5, col="#D55E00", fill="white") +
	scale_x_continuous("Time since introduction (years)", expand=c(0, 0)) +
	scale_y_log10("Prevalence", breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
	scale_color_manual(values=c("#D55E00", "#0072B2")) +
	coord_cartesian(xlim=c(-1, 5), ylim=c(5e-8, 1)) +
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(size=1),
		# axis.line = element_line(),
		legend.position = "none",
		legend.title = element_blank(),
		axis.title.y = element_blank()
	)

figure4_rsv <- ggarrange(g1, g2, nrow=1)

ggsave("figure3_rsv.pdf", figure4_rsv, width=5, height=3)

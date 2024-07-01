# https://www.cell.com/cell/pdf/S0092-8674(23)01076-0.pdf

simulate_meijers <- function(R0=2.5,
														 T11=1,
														 tau=5.5,
														 lambda=3,
														 T50=log10(0.2*94),
														 T0=log10(94),
														 tau_decay = 90,
														 i0=1e-6,
														 tmax=365*30) {
	# C11 <- 1/(1+exp(- lambda * (T0 - log10(T11) - T50 - 1:tmax/tau_decay)))
	
	ivec <- rep(0, tmax)
	ivec[1] <- i0
	Fvec <- rep(0, tmax)
	Fvec[1] <- log(R0)/tau
	
	for (i in 2:tmax) {
		C11 <- sum(1/(1+exp(- lambda * (T0 - log10(T11) - T50 - log10(exp(1:(i-1)/tau_decay))))) * rev(ivec[1:(i-1)]))
		
		Fvec[i] <- (log(R0) - C11)/tau
		
		ivec[i] <- ivec[i-1] + Fvec[i] * ivec[i-1]
	}
	
	data.frame(
		time=1:tmax,
		ivec=ivec,
		Fvec=Fvec,
		Rvec=exp(Fvec*tau)
	)
}

simulate_meijers_niche <- function(R_1=2.5,
																	 T11=1,
																	 T21=1.8,
																	 tau1=5.5,
																	 R_2=2.5*1.7,
																	 T22=1,
																	 T12=1.8,
																	 tau2=5.5,
																	 lambda=3,
																	 T50=log10(0.2*94),
																	 T0=log10(94),
																	 tau_decay = 90,
																	 tmax=365*50,
																	 invmin=365*45,
																	 invmax=365*50,
																	 system="SARS-CoV-2 WT, Alpha") {
	out1 <- simulate_meijers(R0=R_1, T11=T11, tau=tau1, tmax=tmax,
													 lambda=lambda, T50=T50, T0=T0, tau_decay=tau_decay)
	out2 <- simulate_meijers(R0=R_2, T11=T22, tau=tau2, tmax=tmax,
													 lambda=lambda, T50=T50, T0=T0, tau_decay=tau_decay)
	
	out1_filter <- out1[out1$time > invmin & out1$time < invmax,]
	out2_filter <- out2[out2$time > invmin & out2$time < invmax,]
	
	invt <- invmax-invmin-1
	
	C11 <- sum(1/(1+exp(- lambda * (T0 - log10(T11) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out1_filter$ivec))
	C21 <- sum(1/(1+exp(- lambda * (T0 - log10(T21) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out1_filter$ivec))
	C22 <- sum(1/(1+exp(- lambda * (T0 - log10(T22) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out2_filter$ivec))
	C12 <- sum(1/(1+exp(- lambda * (T0 - log10(T12) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out2_filter$ivec))
	
	S_11 <- exp(-C11)
	S_21 <- exp(-C21)
	
	S_22 <- exp(-C22)
	S_12 <- exp(-C12)
	
	R_1 <- R_1
	R_2 <- R_2
	
	niche <- mean(sqrt(S_11*S_22/(S_12*S_21)))
	
	fitness <- mean(R_2/R_1*sqrt(S_22*S_21/(S_11*S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=fitness,
		system=system,
		meanR0=mean(sqrt(R_1*R_2)),
		nulldiff=pmax(mean(sqrt(R_1/R_2)), mean(sqrt(R_2/R_1)))
	)	
}

simulate_meijers_vacc <- function(R0=2.5,
																	T11=1,
																	vprop=0.3,
																	bprop=0,
																	TV=1,
																	TB=1,
																	tau=5.5,
																	lambda=3,
																	T50=log10(0.2*94),
																	T0=log10(94),
																	T0_vacc=log10(223),
																	T0_boos=log10(223*4),
																	tau_decay = 90,
																	i0=1e-6,
																	tmax=365*30) {
	ivec <- rep(0, tmax)
	ivec[1] <- i0
	Fvec <- rep(0, tmax)
	Fvec[1] <- log(R0)/tau
	
	for (i in 2:tmax) {
		C11 <- sum(1/(1+exp(- lambda * (T0 - log10(T11) - T50 - log10(exp(1:(i-1)/tau_decay))))) * rev(ivec[1:(i-1)]))
		Cv <- 1/(1+exp(- lambda * (T0_vacc - log10(TV) - T50))) * vprop
		Cb <- 1/(1+exp(- lambda * (T0_boos - log10(TB) - T50))) * bprop
		
		Fvec[i] <- (log(R0) - C11 - Cv - Cb)/tau
		
		ivec[i] <- ivec[i-1] + Fvec[i] * ivec[i-1]
	}
	
	data.frame(
		time=1:tmax,
		ivec=ivec,
		Fvec=Fvec,
		Rvec=exp(Fvec*tau)
	)
}

simulate_meijers_vacc_niche <- function(R_1=2.5,
																	 T11=1,
																	 T21=1.8,
																	 tau1=5.5,
																	 R_2=2.5*1.7,
																	 T22=1,
																	 T12=1.8,
																	 tau2=5.5,
																	 vprop=0.3,
																	 bprop=0,
																	 TV1=1,
																	 TV2=1.8,
																	 TB1=1,
																	 TB2=1.2,
																	 lambda=3,
																	 T50=log10(0.2*94),
																	 T0=log10(94),
																	 T0_vacc=log10(223),
																	 T0_boos=log10(223*4),
																	 tau_decay = 90,
																	 tmax=365*50,
																	 invmin=365*45,
																	 invmax=365*50,
																	 system="SARS-CoV-2 WT, Alpha") {
	out1 <- simulate_meijers_vacc(R0=R_1, T11=T11, tau=tau1, tmax=tmax,
																vprop=vprop, bprop=bprop, TV=TV1, TB=TB1,
																lambda=lambda, T50=T50, T0=T0, 
																T0_vacc=T0_vacc, T0_boos=T0_boos, tau_decay=tau_decay)
	
	out2 <- simulate_meijers_vacc(R0=R_2, T11=T22, tau=tau2, tmax=tmax,
																vprop=vprop, bprop=bprop, TV=TV2, TB=TB2,
																lambda=lambda, T50=T50, T0=T0, 
																T0_vacc=T0_vacc, T0_boos=T0_boos, tau_decay=tau_decay)
	
	out1_filter <- out1[out1$time > invmin & out1$time < invmax,]
	out2_filter <- out2[out2$time > invmin & out2$time < invmax,]
	
	invt <- invmax-invmin-1
	
	C11 <- sum(1/(1+exp(- lambda * (T0 - log10(T11) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out1_filter$ivec))
	C21 <- sum(1/(1+exp(- lambda * (T0 - log10(T21) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out1_filter$ivec))
	C22 <- sum(1/(1+exp(- lambda * (T0 - log10(T22) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out2_filter$ivec))
	C12 <- sum(1/(1+exp(- lambda * (T0 - log10(T12) - T50 - log10(exp((1:invt)/tau_decay))))) * rev(out2_filter$ivec))
	
	Cv1 <- 1/(1+exp(- lambda * (T0_vacc - log10(TV1) - T50))) * vprop 
	Cv2 <- 1/(1+exp(- lambda * (T0_vacc - log10(TV2) - T50))) * vprop 
	
	Cb1 <- 1/(1+exp(- lambda * (T0_boos - log10(TB1) - T50))) * bprop 
	Cb2 <- 1/(1+exp(- lambda * (T0_boos - log10(TB2) - T50))) * bprop 
	
	S_11 <- exp(-C11-Cv1-Cb1)
	S_21 <- exp(-C21-Cv2-Cb2)
	
	S_22 <- exp(-C22-Cv2-Cb2)
	S_12 <- exp(-C12-Cv1-Cb1)
	
	R_1 <- R_1
	R_2 <- R_2
	
	niche <- mean(sqrt(S_11*S_22/(S_12*S_21)))
	
	fitness <- mean(R_2/R_1*sqrt(S_22*S_21/(S_11*S_12)))
	
	data.frame(
		nichediff=1/niche,
		fitnessdiff=fitness,
		system=system,
		meanR0=mean(sqrt(R_1*R_2)),
		nulldiff=pmax(mean(sqrt(R_1/R_2)), mean(sqrt(R_2/R_1)))
	)	
}

simulate_meijers_nstrain <- function(R0=c(2.5, 3.33, 3.3),
																		 Tmat=matrix(c(1, 1.8, 3.2, 1.8, 1, 2.8, 3.2, 3.5, 1), nrow=3, ncol=3),
																		 tau=c(7.8, 5.5, 3.5),
																		 lambda=3,
																		 T50=log10(0.2*94),
																		 T0=log10(94),
																		 tau_decay = 90,
																		 N=1,
																		 i0=c(1e-6, 1e-6, 1e-6),
																		 intro=c(1, 150, 300),
																		 tmax=600) {
	# C11 <- 1/(1+exp(- lambda * (T0 - log10(T11) - T50 - 1:tmax/tau_decay)))
	
	nstrain <- length(R0)
	
	imat <- matrix(0, nrow=tmax, ncol=nstrain)
	Fmat <- matrix(0, nrow=tmax, ncol=nstrain)
	
	imat[1,intro==1] <- i0[intro==1]
		
	for (t in 2:tmax) {
		Cmat <- matrix(0, nstrain, nstrain)
		
		for (x in 1:nstrain) {
			for (y in 1:nstrain) {
				Cmat[x, y] <- sum(1/(1+exp(- lambda * (T0 - log10(Tmat[x,y]) - T50 - log10(exp(1:(t-1)/tau_decay))))) * rev(imat[1:(t-1),y]))
				
			}
		}
		
		Fmat[t,] <- (log(R0)-rowSums(Cmat))/tau
		imat[t,] <- imat[t-1,] + Fmat[t,] * imat[t-1,]
		
		if (any(intro==t)) {
			imat[t,intro==t] <- i0[intro==t]
		}
	}
	
	list(
		time=1:tmax,
		imat=imat,
		Fmat=Fmat
	)
}

simulate_meijers_nstrain_vacc <- function(R0=c(2.5, 3.33, 3.3),
																					Tmat=matrix(c(1, 1.8, 3.2, 1.8, 1, 2.8, 3.2, 3.5, 1), nrow=3, ncol=3),
																					tau=c(7.8, 5.5, 3.5),
																					lambda=3,
																					T50=log10(0.2*94),
																					T0=log10(94),
																					tau_decay = 90,
																					N=1,
																					vprop=0.3,
																					bprop=0,
																					TV=1,
																					TB=1,
																					T0_vacc=log10(223),
																					T0_boos=log10(223*4),
																					i0=c(1e-6, 1e-6, 1e-6),
																					intro=c(1, 150, 300),
																					tmax=600) {
	Cv <- 1/(1+exp(- lambda * (T0_vacc - log10(TV) - T50))) * vprop
	Cb <- 1/(1+exp(- lambda * (T0_boos - log10(TB) - T50))) * bprop
	
	nstrain <- length(R0)
	
	imat <- matrix(0, nrow=tmax, ncol=nstrain)
	Fmat <- matrix(0, nrow=tmax, ncol=nstrain)
	
	imat[1,intro==1] <- i0[intro==1]
	
	for (t in 2:tmax) {
		Cmat <- matrix(0, nstrain, nstrain)
		
		for (x in 1:nstrain) {
			for (y in 1:nstrain) {
				Cmat[x, y] <- sum(1/(1+exp(- lambda * (T0 - log10(Tmat[x,y]) - T50 - log10(exp(1:(t-1)/tau_decay))))) * rev(imat[1:(t-1),y]))
				
			}
		}
		
		Fmat[t,] <- (log(R0)-rowSums(Cmat)-Cv-Cb)/tau
		imat[t,] <- imat[t-1,] + Fmat[t,] * imat[t-1,]
		
		if (any(intro==t)) {
			imat[t,intro==t] <- i0[intro==t]
		}
	}
	
	list(
		time=1:tmax,
		imat=imat,
		Fmat=Fmat
	)
}

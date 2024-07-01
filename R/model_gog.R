model_gog <- function(t, y, param) {
	with(param, {
		beta <- R0 * (1 + phi * cos(2*pi*t/365)) * (gamma+mu)
		
		S <- y[1:n]
		I <- y[(n+1):(2*n)]
		
		deriv <- c(
			mu - S * (sigmamat %*% (beta * I)) - mu * S,
			beta * S * I - (gamma+mu) * I
		)
		
		total <- sum(I)
		mean <- sum(1:n*I)/total
		
		list(deriv,
				 total=total,
				 mean=mean)
	})
}

simulate_gog <- function(R0=1.8,
												 Reff=1.28,
												 phi=0.15,
												 gamma=1/5,
												 mu=1/30/365,
												 s=0.07,
												 n=300,
												 nscale=10, ## each strain is 1/10 antigenic units apart
												 tvec=(-180):(185),
												 yini) {
	sigmamat <- pmax(1-outer(1:n, 1:n, function(x, y) abs(x-y)*s/nscale),0)
	
	if (missing(yini)) {
		yini <-  c(1-(1-Reff/R0) * sigmamat[,1], 1e-6, rep(0, n-1))
		names(yini)[1:n] <- paste0("S", 1:n)
		names(yini)[(n+1):(2*n)] <- paste0("I", 1:n)		
	}
	
	parms <- list(
		R0=R0,
		phi=phi,
		gamma=gamma, 
		mu=mu,
		n=n,
		sigmamat=sigmamat
	)
	
	out <- as.data.frame(deSolve::rk4(yini, tvec, model_gog, parms))
	
	out
}

simulate_gog_stochastic  <- function(R0=1.8,
																		 Reff=1.28,
																		 phi=0.15,
																		 gamma=1/5,
																		 mu=1/30/365,
																		 s=0.07,
																		 n=300,
																		 nscale=10, ## each strain is 1/10 antigenic units apart
																		 delta_mean=1.3,
																		 delta_sd=0.2,
																		 tmax=20,
																		 yini,
																		 seed=101) {
	set.seed(seed)
	sigmamat <- pmax(1-outer(1:n, 1:n, function(x, y) abs(x-y)*s/nscale),0)
	
	yini <-  c(Reff/R0, rep(1, n-1), 1e-6, rep(0, n-1))
	names(yini)[1:n] <- paste0("S", 1:n)
	names(yini)[(n+1):(2*n)] <- paste0("I", 1:n)	
	
	currstrain <- 1
	simlist <- vector('list', tmax+1)
	
	simlist[[1]] <- simulate_gog(R0=R0, Reff=Reff, phi=phi, gamma=gamma, mu=mu, s=s, n=n, nscale=nscale,
															 yini=yini,
															 tvec=(-180):(185))
	
	for (i in 1:tmax) {
		print(i)
		# 1/shape if squared CV
		# shape = 1/CV
		# mean = shape * scale
		# scale = mean/shape = sd
		
		currstrain <- newstrain <- currstrain + round(rgamma(1, shape=delta_mean/delta_sd, scale=delta_sd)*nscale)
		
		yini <- tail(head(unlist(tail(simlist[[i]],1)), -2), -1)
		yini[n+newstrain] <- 1e-6
		
		simlist[[i]] <- head(simlist[[i]], -1)
		
		simlist[[i+1]] <- simulate_gog(R0=R0, Reff=Reff, phi=phi, gamma=gamma, mu=mu, s=s, n=n, nscale=nscale,
																	 yini=yini,
																	 tvec=(i*365-180):(i*365+185))
	}
	
	out <- do.call("rbind", simlist)
	
	out
}

model_gog_SI <- function(t, y, param) {
	with(param, {
		beta <- R0 * (1 + phi * cos(2*pi*t/365)) * (gamma+mu)
		
		S <- y[1:n]
		I <- y[(n+1):(2*n)]
		
		deriv <- c(
			mu - S * (sigmamat %*% (beta * I)) - mu * S + gamma * (sigmamat %*% I),
			beta * S * I - (gamma+mu) * I
		)
		
		total <- sum(I)
		mean <- sum(1:n*I)/total
		
		list(deriv,
				 total=total,
				 mean=mean)
	})
}

simulate_gog_SI <- function(R0=1.8,
												 Reff=1.28,
												 phi=0.15,
												 gamma=1/5,
												 mu=1/30/365,
												 s=0.07,
												 n=300,
												 nscale=10, ## each strain is 1/10 antigenic units apart
												 tvec=(-180):(185),
												 yini) {
	sigmamat <- pmax(1-outer(1:n, 1:n, function(x, y) abs(x-y)*s/nscale),0)
	
	if (missing(yini)) {
		yini <-  c(1-(1-Reff/R0) * sigmamat[,1], 1e-6, rep(0, n-1))
		names(yini)[1:n] <- paste0("S", 1:n)
		names(yini)[(n+1):(2*n)] <- paste0("I", 1:n)		
	}
	
	parms <- list(
		R0=R0,
		phi=phi,
		gamma=gamma, 
		mu=mu,
		n=n,
		sigmamat=sigmamat
	)
	
	out <- as.data.frame(deSolve::rk4(yini, tvec, model_gog_SI, parms))
	
	out
}

model_gog_SIRS <- function(t, y, param) {
	with(param, {
		beta <- R0 * (1 + phi * cos(2*pi*t/365)) * (gamma+mu)
		
		S <- y[1:n]
		I <- y[(n+1):(2*n)]
		
		deriv <- c(
			mu - S * (sigmamat %*% (beta * I)) - mu * S + delta * (sigmamat %*% (1-S-I)),
			beta * S * I - (gamma+mu) * I
		)
		
		total <- sum(I)
		mean <- sum(1:n*I)/total
		
		list(deriv,
				 total=total,
				 mean=mean)
	})
}

simulate_gog_SIRS <- function(R0=1.8,
															Reff=1.28,
															phi=0.15,
															gamma=1/5,
															delta=1/365,
															mu=1/30/365,
															s=0.07,
															n=300,
															nscale=10, ## each strain is 1/10 antigenic units apart
															tvec=(-180):(185),
															yini) {
	sigmamat <- pmax(1-outer(1:n, 1:n, function(x, y) abs(x-y)*s/nscale),0)
	
	if (missing(yini)) {
		yini <-  c(1-(1-Reff/R0) * sigmamat[,1], 1e-6, rep(0, n-1))
		names(yini)[1:n] <- paste0("S", 1:n)
		names(yini)[(n+1):(2*n)] <- paste0("I", 1:n)		
	}
	
	parms <- list(
		R0=R0,
		phi=phi,
		gamma=gamma, 
		delta=delta,
		mu=mu,
		n=n,
		sigmamat=sigmamat
	)
	
	out <- as.data.frame(deSolve::rk4(yini, tvec, model_gog_SIRS, parms))
	
	out
}
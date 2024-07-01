library(deSolve)
source("../R/model_kissler.R")
source("../R/model_white.R")
source("../R/model_rohani.R")
source("../R/model_koelle.R")
source("../R/model_yang.R")
source("../R/model_bhattacharyya.R")

data_kissler_niche <- simulate_kissler_niche()
data_white_niche <- simulate_white_niche()
data_rohani_niche <- simulate_rohani_niche(tmax=200,
																					 invmin=190, invmax=200)
# https://royalsocietypublishing.org/doi/pdf/10.1098/rspb.2005.3454
data_rohani_niche2 <- simulate_rohani_niche(b01=1250/52,
																						sigma1=7/8,
																						gamma1=7/5,
																						b02=49/11,
																						sigma2=7/9,
																						gamma2=7/11,
																						tmax=200,
																						invmin=190, invmax=200,
																						system="Measles, Rubella")
data_rohani_niche3 <- simulate_rohani_niche(b01=1250/52,
																						sigma1=7/8,
																						gamma1=7/5,
																						b02=77/10,
																						sigma2=7/10,
																						gamma2=7/10,
																						tmax=200,
																						invmin=190, invmax=200,
																						system="Measles, Chickenpox")
data_koelle_niche <- simulate_koelle_niche()
data_yang_niche <- simulate_yang_niche(system="Flu A(H1N1), A(H3N2)")
data_yang_niche2 <- simulate_yang_niche(R01=1.44,
																				R02=1.43,
																				D1=2.64/7,
																				D2=3.09/7,
																				L1=3.12*52,
																				L2=3.08*52,
																				c12=0.32,
																				c21=0.23,
																				system="Flu A(H1N1), B")
data_yang_niche3 <- simulate_yang_niche(R01=1.60,
																				R02=1.43,
																				D1=3.0/7,
																				D2=3.09/7,
																				L1=2.28*52,
																				L2=3.08*52,
																				c12=0.24,
																				c21=0.11,
																				system="Flu A(H3N2), B")

data_bhattacharyya_niche1 <- simulate_bhattacharyya_niche(b01=3.4,
																													b02=2.9,
																													theta1=0.4,
																													theta2=0.33,
																													phi1=-17.5,
																													phi2=-29.5,
																													epsilon12=0.9,
																													epsilon21=0.54,
																													system="RSV, HPIV-1")

data_bhattacharyya_niche2 <- simulate_bhattacharyya_niche(b01=3.4,
																													b02=2.4,
																													theta1=0.4,
																													theta2=0.28,
																													phi1=-17.99,
																													phi2=-20.001,
																													epsilon12=0.9,
																													epsilon21=0.796,
																													system="RSV, HPIV-2")

data_bhattacharyya_niche3 <- simulate_bhattacharyya_niche(b01=3.4,
																													b02=2.0,
																													theta1=0.4,
																													theta2=0.31,
																													phi1=-18,
																													phi2=-11,
																													epsilon12=0.87,
																													epsilon21=0.63,
																													system="RSV, HPIV-3")

data_bhattacharyya_niche4 <- simulate_bhattacharyya_niche(b01=3.4,
																													b02=3.9,
																													theta1=0.4,
																													theta2=0.3,
																													phi1=0.005,
																													phi2=4.99,
																													epsilon12=0.92,
																													epsilon21=0.45,
																													system="RSV, hMPV")

save("data_kissler_niche", "data_white_niche", "data_rohani_niche",
		 "data_rohani_niche2", "data_rohani_niche3",
		 "data_koelle_niche",
		 "data_yang_niche", "data_yang_niche2", "data_yang_niche3",
		 "data_bhattacharyya_niche1", "data_bhattacharyya_niche2", "data_bhattacharyya_niche3", "data_bhattacharyya_niche4",
		 file="simulate_niche.rda")

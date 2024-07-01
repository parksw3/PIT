library(Racmacs)
options(RacOptimizer.num_cores = 1)

path_to_titer_file <- system.file("extdata/h3map2004_hitable.csv", package = "Racmacs")
titer_table        <- read.titerTable(path_to_titer_file)

target_name <- c("HK/1/68", "EN/42/72", "VI/3/75", "TE/1/77",
                 "BA/1/79", "SI/2/87", "BE/353/89", "BE/32/92",
                 "WU/359/95", "SY/5/97", "FU/411/02")

map <- acmap(
  titer_table = titer_table
)

map <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 2,
  number_of_optimizations = 500,
  minimum_column_basis    = "none"
)

ag <- agCoords(map, 1)

data <- ag[match(target_name, rownames(ag)),]

s <- 0.07
R0 <- 1.8

reslist <- vector('list', nrow(data)-1)

for (i in 1:(nrow(data)-1)) {
  x_ag1 <- data[i,1]
  x_ag2 <- data[i,2]
  
  y_ag1 <- data[i+1,1]
  y_ag2 <- data[i+1,1]
  
  dist <- sqrt((y_ag1-x_ag1)^2+(y_ag2-x_ag2)^2)
  
  sigma <- 1-s*dist
  
  S_11 <- 1/R0
  S_21 <- 1/R0 + (1-1/R0) * (1-sigma)
  
  S_22 <- 1/R0
  S_12 <- 1/R0 + (1-1/R0) * (1-sigma)
  
  niche <- sqrt(S_11*S_22/(S_12*S_21))
  
  fitness <- sqrt(S_22*S_21/(S_11*S_12)) ## == 1
  
  reslist[[i]] <- data.frame(
    dist=dist,
    sigma=sigma,
    nichediff=1/niche,
    fitnessdiff=pmax(fitness, 1/fitness)
  )
}

smith_niche <- reslist %>%
  bind_rows

save("smith_niche", file="smith_niche.rda")

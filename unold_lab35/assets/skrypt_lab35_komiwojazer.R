# clean old data
rm(list=ls())
dev.off(dev.list()["RStudioGD"])

# load libraries
require("GA")
require("globalOptTests")
require("rgl")
require("TSP")
require("psoptim")

# TSP with GA ----

instances <- c("eil51", "eil76", "eil101")
best_solutions <- c(426, 538, 629)
found_solutions <- c()
solutions_quality <- c()

for (i in 1:length(instances)) {
  
  fileName = paste("examples/", instances[i], ".tsp", sep="")
  graphTitle = paste("TSPLIB: ", instances[i], sep="")
  
  drill <- read_TSPLIB(system.file(fileName, package = "TSP"))
  D <- as.matrix(dist(drill, method = "euclidean"))
  N <- max(dim(D))
  
  tourLength <- function(tour, distMatrix) {
    tour <- c(tour, tour[1])
    route <- embed(tour, 2)[,2:1]
    sum(distMatrix[route])
  }
  
  fit <- function(tour, distMatrix) 1/tourLength(tour, distMatrix)
  
  GA <- ga(type = "permutation", 
           fitness = fit, 
           distMatrix = D, 
           min = 1, 
           max = N, 
           maxiter=200, 
           pmutation=0.2, 
           run=300, 
           optim = FALSE)
  
  tour <- GA@solution[1, ]

  plot(drill, tour, cex=.6, col = "red", pch= 3, main = graphTitle)
  
  tl <- tourLength(tour, D)
  ol <- best_solutions[i]
  
  found_solutions <- c(found_solutions, tl)
  solutions_quality <- c(solutions_quality, (ol/tl) * 100)

}

png(file = "tsp_results.png", width=600, height=400, units="px")
plot(1:length(instances), xaxt = "n", solutions_quality, col="red", main="Wyniki dla ró¿nych instancji", type="l", xlab="instancje", ylab="jakoœæ rozwi¹zañ [%]")
axis(1, at=1:length(instances), labels=instances)
dev.off()


# PSO tests ----

n <- 50
m.l <- 50
w <- 0.95
c1 <- 0.2
c2 <- 0.2
xmin <- c(-5.12, -5.12)
xmax <- c(5.12, 5.12)
vmax <- c(4, 4)

g <- function(x){  
  -(20 + x[,1]^2 + x[,2]^2 - 10*(cos(2*pi*x[,1]) + cos(2*pi*x[,2])))
}

psoptim(FUN=g, n=n, max.loop=m.l, w=w, c1=c1, c2=c2,
        xmin=xmin, xmax=xmax, vmax=vmax, seed=5, anim=FALSE)

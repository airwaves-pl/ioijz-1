
# whole source code is located here:
# https://github.com/cran/GA/tree/master/R

# clean old data
rm(list=ls())
dev.off(dev.list()["RStudioGD"])

# load libraries
require("GA")
require("globalOptTests")
require("rgl")

# ATSP example ----

require(GA)

data("eurodist", package = "datasets")
D <- as.matrix(eurodist)
N <- max(dim(D))

tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}

fit <- function(tour, distMatrix) 1/tourLength(tour, distMatrix)

GA <- ga(type = "permutation", fitness = fit, distMatrix = D, min = 1, max = N, maxiter=2000, pmutation=0.2, run=500, optim = FALSE)

summary(GA)

apply(GA@solution, 1, tourLength, D)

mds <- cmdscale(eurodist)
x <- mds[, 1]
y <- -mds[, 2]
plot(x, y, type = "n", asp = 1, xlab = "", ylab = "")
tour <- GA@solution[1, ]
tour <- c(tour, tour[1])
n <- length(tour)
arrows(x[tour[-n]], y[tour[-n]], x[tour[-1]], y[tour[-1]],length = 0.15, angle = 25, col = "steelblue", lwd = 2)
text(x, y, labels(eurodist), cex=0.8)



# TSP tests ----

require(TSP)

instances <- c("u159", "u574", "u724", "u1060", "u1432", "u1817", "u2152", "u2319")
best_solutions <- c(42080, 36905, 41910, 224094, 152970, 57201, 64253, 234256)
found_solutions <- c()
solutions_quality <- c()

for (i in 1:length(instances)) {
  
  fileName = paste("examples/", instances[i], ".tsp", sep="")
  graphTitle = paste("TSPLIB: ", instances[i], sep="")
  
  ## Drilling problem from TSP
  drill <- read_TSPLIB(system.file(fileName, package = "TSP"))
  tour <- solve_TSP(drill, method = "nn", two_opt = TRUE)
  plot(drill, tour, cex=.6, col = "red", pch= 3, main = graphTitle)
  
  tl <- tour_length(tour)
  ol <- best_solutions[i]
  
  found_solutions <- c(found_solutions, tl)
  solutions_quality <- c(solutions_quality, ol/tl)

}

#png(file = "tsp_results.png", width=600, height=400, units="px")
plot(1:length(instances), xaxt = "n", solutions_quality, col="red", main="Summary", type="l", xlab="instancje", ylab="jakość rozwiązań")
axis(1, at=1:length(instances), labels=instances)
#dev.off()



# PSO tests ----

require(psoptim)

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

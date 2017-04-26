# clean old data
rm(list=ls())
dev.off(dev.list()["RStudioGD"])
 
# load libraries
require("GA")
require("globalOptTests")
require("rgl")
require("TSP")
require("psoptim")

numberOfMeasurements <- 2 #15


# TSP with GA ----

# instances to test and best known solutions
instances <- c("eil51", "eil76", "eil101")
best_solutions <- c(426, 538, 629)
colors <- c("red", "green", "blue")

tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}

fit <- function(tour, distMatrix) 1/tourLength(tour, distMatrix)

performTest <- function(testName, graphMain, graphXLab, 
                        sequenceType, sequence, 
                        popsize=50, pcrossover=0.8, 
                        pmutation=0.1, maxiter=100,
                        optim=FALSE) {
  
  solution_qualities <- c()
  
  # each instance as separate serie
  for (i in 1:length(instances)) {
    
    fileName = paste("examples/", instances[i], ".tsp", sep="")
    graphTitle = paste("TSPLIB: ", instances[i], sep="")
    
    drill <- read_TSPLIB(system.file(fileName, package = "TSP"))
    D <- as.matrix(dist(drill, method = "euclidean"))
    N <- max(dim(D))
    
    solution_quality <- c()
    
    bestTour <- NA
    bestTourLength <- .Machine$integer.max
    averageLength <- 0
    
    for (s in 1:length(sequence)) {
      
      for (n in 1:numberOfMeasurements) {
        
        message(paste("Instancja: ", i))
        message(paste("Sekwencja: ", s))
        message(paste("Pomiar: ", n))
        
        GA <- ga(type = "permutation", 
                 fitness = fit, 
                 distMatrix = D, 
                 min = 1, 
                 max = N, 
                 popSize = if (sequenceType == "popsize") sequence[s] else popsize, 
                 pcrossover = if (sequenceType == "pcrossover") sequence[s] else pcrossover, 
                 pmutation = if (sequenceType == "pmutation") sequence[s] else pmutation, 
                 maxiter = if (sequenceType == "maxiter") sequence[s] else maxiter,
                 optim = optim)
        
        tour <- GA@solution[1, ]
        tl <- tourLength(tour, D)
        
        if (tl < bestTourLength) {
          bestTourLength <- tl
          bestTour <- tour
        }
        
        averageLength <- averageLength + (tl - averageLength) / n
        
      }
      
      solution_quality <- c(solution_quality, 
                            (best_solutions[i]/averageLength) * 100)
      
    }

    plot(drill, bestTour, cex=.6, col = "red", pch=3, main = graphTitle)
    
    solution_qualities <- c(solution_qualities, solution_quality)
    
  }
  
  
  qualities = matrix(solution_qualities, 
                  nrow=length(instances), ncol=length(sequence), byrow = TRUE)
  
  # save graph with measurement series to file
  png(file = paste(testName, ".png", sep=""), width=600, height=400, units="px")
  plot(0, 0, main=graphMain, 
       ylim=c(0,100),
       xlim=c(min(sequence),max(sequence)),
       type="n", xlab=graphXLab, ylab="jakość rozwiązań [%]")
  for (i in 1:length(instances)) {
    lines(sequence, qualities[i,], col = colors[i], type = 'l')
  }
  legend("topright", instances, lwd=rep(2,length(instances)), lty=rep(1,length(instances)), col=colors)
  dev.off()
  
}

performTest(testName = "tsp_pop", 
            graphMain = "Pomiary dla różnych rozmiarów populacji", 
            graphXLab = "rozmiar populacji", 
            sequenceType = "popsize", sequence = seq(50, 500, 50))

performTest(testName = "tsp_pop_hyb", 
            graphMain = "Pomiary dla różnych rozmiarów populacji (hybrydowy)", 
            graphXLab = "rozmiar populacji", 
            sequenceType = "popsize", sequence = seq(50, 500, 50), optim = TRUE)

performTest(testName = "tsp_mut", 
            graphMain = "Pomiary dla różnych p. mutacji", 
            graphXLab = "p. mutacji", 
            sequenceType = "pmutation", sequence = seq(0, 1, 0.1))

performTest(testName = "tsp_cross", 
            graphMain = "Pomiary dla różnych p. krzyżowania", 
            graphXLab = "p. krzyżowania", 
            sequenceType = "pcrossover", sequence = seq(0, 1, 0.1))




# PSO tests ----

n <- 500
m.l <- 50
w <- 0.95
c1 <- 0.2
c2 <- 0.2
xmin <- c(-5.12, -5.12)
xmax <- c(5.12, 5.12)
vmax <- c(4, 4)

g <- function(x){  
  -(200 + x[,1]^2 + x[,2]^2 + cos(2*pi*x[,2]))
}

psoptim(FUN=g, n=n, max.loop=m.l, w=w, c1=c1, c2=c2,
        xmin=xmin, xmax=xmax, vmax=vmax, seed=5, anim=TRUE)

#determine permutation by number
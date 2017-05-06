# clean old data
rm(list=ls())
dev.off(dev.list()["RStudioGD"])
 
# load libraries
require("GA")
require("globalOptTests")
require("rgl")
require("TSP")
require("psoptim")

numberOfMeasurements <- 15

# instances for testing and best known solutions
instances <- c("eil51", "eil76", "eil101")
best_solutions <- c(426, 538, 629)
colors <- c("red", "green", "blue")

tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}

fit <- function(tour, distMatrix) 1/tourLength(tour, distMatrix)

customMutation <- function(object, parent, ...) {
  # insertion mutation
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sample(1:n, size = 1)
  pos <- sample(1:(n-1), size = 1)
  i <- c(setdiff(1:pos,m), m, setdiff((pos+1):n,m))
  mutate <- parent[i]

  # displacement mutation
  parent <- mutate
  m <- sort(sample(1:n, size = 2))
  m <- seq(m[1], m[2], by = 1)
  l <- max(m)-min(m)+1
  pos <- sample(1:max(1,(n-l)), size = 1)
  i <- c(setdiff(1:n,m)[1:pos], m, setdiff(1:n,m)[-(1:pos)])
  mutate <- parent[na.omit(i)]
  
  # scramble mutation
  parent <- mutate
  m <- sort(sample(1:n, size = 2))
  m <- seq(min(m), max(m), by = 1)
  m <- sample(m, replace = FALSE)
  i <- c(setdiff(1:min(m),m), m, setdiff(max(m):n,m))
  mutate <- parent[i]
  return(mutate)
} 

performTest <- function(testName, graphMain, graphXLab, 
                        sequenceType, sequence, 
                        popsize=50, pcrossover=0.8, 
                        pmutation=0.1, maxiter=100, mutation = NULL) {
  
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
                 mutation = if (is.null(mutation)) gaControl("permutation")$mutation else mutation) 
        
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

    png(file = paste(testName, "_", instances[i], ".png", sep=""), width=600, height=400, units="px")
    plot(drill, bestTour, cex=.6, col = "red", pch=3, main = graphTitle)
    dev.off()
    
    solution_qualities <- c(solution_qualities, solution_quality)
  }
  
  qualities = matrix(solution_qualities, nrow=length(instances), ncol=length(sequence), byrow = TRUE)
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

performTest(testName = "tsp_pop", graphMain = "Pomiary dla różnych rozmiarów populacji", 
            graphXLab = "rozmiar populacji", 
            sequenceType = "popsize", sequence = seq(50, 500, 50))

performTest(testName = "tsp_mut", graphMain = "Pomiary dla różnych p. mutacji", 
            graphXLab = "p. mutacji", 
            sequenceType = "pmutation", sequence = seq(0, 1, 0.1))

performTest(testName = "tsp_mut_custom", graphMain = "Pomiary dla różnych p. mutacji (własny op. mutacji)", 
            graphXLab = "p. mutacji", 
            sequenceType = "pmutation", sequence = seq(0, 1, 0.1), mutation = customMutation)

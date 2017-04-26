# initialize ----
# clean old data
rm(list=ls())
dev.off(dev.list()["RStudioGD"])
 
# load libraries
require("GA")
require("globalOptTests")
require("rgl")
require("parallel")
require("doParallel")

# custom functions ----
# mutation function
myMutationFunction <- function(object, parent) {
  # get GA population
  population <- parent <- as.vector(object@population[parent, ])
  
  # calculate randoms
  rnd <- sample(1:length(population), 1)
  rndMinOrMax <- sample(1:10, 1)
  
  # get min and max from population vector
  max_value <- which.max(population)
  min_value <- which.min(population)
 
  # set random element to min value 
  population[rnd] = min_value
 
  return (population);
}
# crossover function
myCrossoverFunction <- function(object, parent) {
  
}

# Settings ----

nOfRuns <- 1 # 15 number of runs to calc avg scores

numOfCores <- FALSE # number of cores to use (FALSE, 1 - n)

# colors and titles for plot series
colors <- c("red", "blue", "purple", "black")
series <- c("Seria 1", "Seria 2", "Seria 3", "Seria 4")

# default parameters for measurements
# each row is a different serie
# [mutations,crossovers,populations,iterations,color]
'params = matrix(
  c(0, 0, 50, 100, 1,
    0, 0.8, 50, 100, 2,
    0.1, 0, 50, 100, 3,
    0.1, 0.8, 50, 100, 4),
  nrow=4, ncol=5, byrow = TRUE)'

params = matrix(
  c(0.1, 0.8, 50, 100, 4),
  nrow=1, ncol=5, byrow = TRUE)

# names of functions from globalOptTests package
functions <- c("Branin")#, "Gulf", "CosMix4", "EMichalewicz", 
	#"Hartman6", "PriceTransistor", "Schwefel", "Zeldasine20")

# graph settings
graphs <- TRUE #true if you want to print graphs
quality <- 100 #number of probes

# sequences of parameters for each serie
#mutationTests <- seq(0, 1, 0.1)
#crossoverTests <- seq(0, 1, 0.1)
#populationTests <- seq(10, 100, 5)
#iterationTests <- seq(10, 200, 10)
#elitismTests <- seq(0, 1, 0.1)

# test only mutation (precisely)
mutationTests <- seq(0, 1, 0.01)
crossoverTests <- c(0.2)
populationTests <- c(30)
iterationTests <- c(5)
elitismTests <- c(0.05)

# hybrid algorithm
hybrid = FALSE
poptim = 0.05 #a value [0,1] specifying the probability of performing a local search at each iteration of GA (def 0.1)
pressel = 0.5 #a value [0,1] specifying the pressure selection (def 0.5)

# Processing ----

customMeasure <- function(fileName, graphName, values, mType, xlab, main) {
  
  gMin <- .Machine$integer.max
  gBest <- NA
  
  # main measurement loop (for each serie and sequence calculate average results)
  temp <- c()
  for (defRow in 1:nrow(params)) {
    averages <- c()
    for (value in values) {
      sum <- 0
      for (i in 1:nOfRuns) {
        GAmin <- ga(type = "real-valued",
            mutation = myMutationFunction,
            #crossover = myCrossoverFunction,
            fitness =  function(xx) -f(xx),
            min = c(B[1,]), max = c(B[2,]),
            popSize = if (mType == "pop") value else params[defRow,3],
            maxiter = if (mType == "itr") value else params[defRow,4],
            pmutation = if (mType == "mut") value else params[defRow,1], 
            pcrossover = if (mType == "crs") value else params[defRow,2],
            elitism = if (mType == "elt") value else max(1, round(params[defRow,3] * 0.05)),
            parallel = numOfCores,
            #hybrid algorithm
            optim = hybrid,
            optimArgs = list (
              poptim = poptim, 
              pressel = pressel))
        solution <- matrix(unlist(GAmin@solution),ncol=dim,byrow=TRUE)
        eval <- f(solution[1,])
        if (eval < gMin) {
          gMin <- eval
          gBest <- GAmin
        }
        sum <- sum + eval
      }
      averages <- c(averages, (sum / nOfRuns))
    }
    temp <- c(temp, averages)
  }
  result <- matrix(c(temp),nrow = nrow(params),ncol = length(values))
  write.table(result, file = paste(funcName, fileName, sep=""), row.names=FALSE, 
              na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    
    # save graph with measurement series to file
    png(file = paste(funcName, graphName, ".png", sep=""), width=600, height=400, units="px")
    plot(0, 0, main=main,
         ylim=c(min(c(temp,globalOpt)),max(c(temp,globalOpt))),
         xlim=c(min(values),max(values)),
         type="n", xlab=xlab, ylab="wartosc")
    abline(globalOpt,0, col="green")
    colorNames <- c()
    seriesNames <- c()
    for (i in 1:nrow(params)) {
      color <- colors[params[i,5]]
      colorNames <- c(colorNames, color)
      seriesNames <- c(seriesNames, series[params[i,5]])
      lines(values, result[i,], col = color, type = 'l')
    }
    legend("topright", seriesNames, lwd=rep(2,nrow(params)), lty=rep(1,nrow(params)), col=colorNames)
    dev.off()
    
    summary(gBest)
    
    # save overview of best found minimum to file
    png(file = paste(funcName, graphName, mType, ".png", sep=""), width=600, height=400, units="px")
    filled.contour(x, y, z, color.palette = jet.colors, nlevels = 24, 
         plot.axes = { axis(1); axis(2);
           points(solution[1,1], solution[1,2], 
                  pch = 3, cex = 5, col = "black", lwd = 2) 
         }
    )
    dev.off()
    
    # save best fitness graph to file
    png(file = paste(funcName, graphName, mType, "fitness", ".png", sep=""), width=600, height=400, units="px")
    plot(gBest)
    dev.off()
  }
}

for (funcName in functions) {

  # get data from globalOptTests package
	dim <- getProblemDimen(funcName)
	B <- matrix(unlist(getDefaultBounds(funcName)),ncol=dim,byrow=TRUE)
	f <- function(xx) goTest(par=c(xx, rep(0, dim-length(xx))), 
							 fnName=funcName, checkDim = TRUE)
	globalOpt <- getGlobalOpt(funcName)

	if (graphs) {
	  # prepare two versions of graphs (interactive and static)
	  xprobes <- abs(B[2,1] - B[1,1]) / quality
	  yprobes <- abs(B[2,2] - B[1,2]) / quality
	  x <- seq(B[1,1], B[2,1], by = xprobes)
	  y <- seq(B[1,2], B[2,2], by = yprobes)
	  z <- outer(x, y, Vectorize(function(x,y) f(c(x,y))))
	  nbcol = 100
	  color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
	  zcol  = cut(z, nbcol)
	  persp3d(x, y, z, theta=50, phi=25, expand=0.75, col=color[zcol],
			  ticktype="detailed",axes=TRUE)
	  
	  png(file = paste(funcName, "1.png", sep=""), width=600, height=400, units="px")
	  persp3D(x, y, z, theta = -45, phi = 20, color.palette = jet.colors)
	  dev.off()
	}
	
	# for each function perform set of measurements
	customMeasure("resultsMutations.csv", "2", mutationTests, "mut", 
		"p. mutacji", "Znalezione minimum dla roznych prawdopodobienstw mutacji")
	#customMeasure("resultsCrossover.csv", "3", crossoverTests, "crs", 
	#	"p. krzyzowania", "Znalezione minimum dla roznych prawdopodobienstw krzyzowania")
	#customMeasure("resultsPopulation.csv", "4", populationTests, "pop", 
	#	"rozmiar populacji", "Znalezione minimum dla roznych rozmiarow populacji")
	#customMeasure("resultsIterations.csv", "5", iterationTests, "itr",
	#	"ilosc iteracji", "Znalezione minimum dla roznych ilosci iteracji")
	#customMeasure("resultsElitism.csv", "6", elitismTests, "elt", 
	#	"elityzm", "Znalezione minimum dla roznych wartosci elityzmu")
}

# whole source code is located here:
# https://github.com/cran/GA/tree/master/R

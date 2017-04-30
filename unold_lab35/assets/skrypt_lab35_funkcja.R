# initialize ----
# clean old data
rm(list=ls())
dev.off(dev.list()["RStudioGD"])

# load libraries
require("GA")
require("globalOptTests")
require("rgl")
require("psoptim")

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

# Settings ----

nOfRuns <- 1 # 30 number of runs to calc avg scores

# colors and titles for plot series
colors <- c("red", "purple")
series <- c("GA", "GA + własna mutacja")

GAWithHybridSeries <- c("GA", "GA + własna mutacja", "Mem", "Mem + własna mutacja")
GAWithHybridColors <- c("red", "purple", "blue", "orange")

# name of function from globalOptTests package
funcName <- "Hartman6"

# graph settings
graphs <- TRUE #true if you want to print graphs
quality <- 100 #number of probes

#hybrid algorithm settings
poptim = 0.05 #a value [0,1] specifying the probability of performing a local search at each iteration of GA (def 0.1)
pressel = 0.5 #a value [0,1] specifying the pressure selection (def 0.5)

# Processing ----

customGAMeasure <- function(values, mType, xlab, main) {
  
  # main measurement loop (for each serie and sequence calculate average results)
  temp <- c()
  for (serie in 1:length(series)) {
    averages <- c()
    for (value in values) {
      sum <- 0
      for (i in 1:nOfRuns) {
        
        message(paste("Seria: ", serie))
        message(paste("Sekwencja: ", value))
        message(paste("Przebieg: ", i))
        
        GAmin <- ga(type = "real-valued",
            mutation = if (serie == 2) myMutationFunction else gaControl("real-valued")$mutation,
            fitness =  function(xx) -f(xx),
            min = c(B[1,]), max = c(B[2,]),
            popSize = if (mType == "pop") value else 50,
            pmutation = if (mType == "mut") value else 0.1)
        solution <- matrix(unlist(GAmin@solution),ncol=dim,byrow=TRUE)
        eval <- f(solution[1,])
        sum <- sum + eval
      }
      averages <- c(averages, (sum / nOfRuns))
    }
    temp <- c(temp, averages)
  }
  result <- matrix(c(temp), nrow = length(series), ncol = length(values))
  
  if (graphs) {
    
    # save graph with measurement series to file
    png(file = paste(funcName, mType, ".png", sep=""), width=600, height=400, units="px")
    plot(0, 0, main=main,
         ylim=c(min(c(temp,globalOpt)),max(c(temp,globalOpt))),
         xlim=c(min(values),max(values)),
         type="n", xlab=xlab, ylab="wartosc")
    abline(globalOpt,0, col="green")
    colorNames <- c()
    seriesNames <- c()
    for (i in 1:length(series)) {
      color <- colors[i]
      colorNames <- c(colorNames, color)
      seriesNames <- c(seriesNames, series[i])
      lines(values, result[i,], col = color, type = 'l')
    }
    legend("topright", seriesNames, lwd=rep(2,length(series)), lty=rep(1,length(series)), col=colorNames)
    dev.off()
    
  }
}

customMeasureGAWithHybrid <- function(values, mType, xlab, main) {
  
  # main measurement loop (for each serie and sequence calculate average results)
  temp <- c()
  for (serie in 1:length(GAWithHybridSeries)) {
    averages <- c()
    for (value in values) {
      sum <- 0
      for (i in 1:nOfRuns) {
        
        message(paste("Seria: ", GAWithHybridSeries[serie]))
        message(paste("Sekwencja: ", value))
        message(paste("Przebieg: ", i))
        
        if(GAWithHybridSeries[serie] == "GA" || GAWithHybridSeries[serie] == "GA + własna funkcja")
        {
          GAmin <- ga(type = "real-valued",
                      mutation = if (serie == 2) myMutationFunction else gaControl("real-valued")$mutation,
                      fitness =  function(xx) -f(xx),
                      min = c(B[1,]), max = c(B[2,]),
                      popSize = if (mType == "pop") value else 50,
                      pmutation = if (mType == "mut") value else 0.1)
        }
        else
        {
          GAmin <- ga(type = "real-valued",
                      mutation = if (serie == 4) myMutationFunction else gaControl("real-valued")$mutation,
                      fitness =  function(xx) -f(xx),
                      min = c(B[1,]), max = c(B[2,]),
                      optim = TRUE,
                      optimArgs = list (
                        poptim = if (mType == "poptim") value else 0.05, 
                        pressel = if (mType == "pressel") value else 0.5))
        }
        
        solution <- matrix(unlist(GAmin@solution),ncol=dim,byrow=TRUE)
        eval <- f(solution[1,])
        sum <- sum + eval
      }
      averages <- c(averages, (sum / nOfRuns))
    }
    temp <- c(temp, averages)
  }
  result <- matrix(c(temp), nrow = length(GAWithHybridSeries), ncol = length(values))
  
  if (graphs) {
    # create standalone graph for each serie
    for (serie in 1:length(GAWithHybridSeries)) {
      legendColors <- rep("darkslateblue", length(GAWithHybridSeries))
      legendColors[serie] = "red"
      # save graph with measurement series to file
      png(file = paste(funcName, mType, serie, ".png", sep=""), width=600, height=400, units="px")
      plot(0, 0, main=main,
           ylim=c(min(c(temp,globalOpt)),max(c(temp,globalOpt))),
           xlim=c(min(values),max(values)),
           type="n", xlab=xlab, ylab="wartosc")
      abline(globalOpt,0, col="green")
      
      lastLine <- NA
      seriesNames <- c()
      for (i in 1:length(GAWithHybridSeries)) {
        seriesNames <- c(seriesNames, GAWithHybridSeries[i])
        if (i != serie)
        {
          lines(values, result[i,], col = "darkslateblue", type = 'l', lwd = 2)
        }
      }
      lines(values, result[serie,], col = "red", type = 'l', lwd = 2)
      
      legend("topright", seriesNames, lwd=rep(2,length(GAWithHybridSeries)), lty=rep(1,length(GAWithHybridSeries)), col = legendColors)
      dev.off()
    }
  }
}

{
# get data from globalOptTests package
dim <- getProblemDimen(funcName)
B <- matrix(unlist(getDefaultBounds(funcName)),ncol=dim,byrow=TRUE)
f <- function(xx) goTest(par=c(xx, rep(0, dim-length(xx))), 
						 fnName=funcName, checkDim = TRUE)
globalOpt <- getGlobalOpt(funcName)

if (graphs) {
  # prepare overview graph
  xprobes <- abs(B[2,1] - B[1,1]) / quality
  yprobes <- abs(B[2,2] - B[1,2]) / quality
  x <- seq(B[1,1], B[2,1], by = xprobes)
  y <- seq(B[1,2], B[2,2], by = yprobes)
  z <- outer(x, y, Vectorize(function(x,y) f(c(x,y))))
  png(file = paste(funcName, "_overview.png", sep=""), width=600, height=400, units="px")
  persp3D(x, y, z, theta = -45, phi = 20, color.palette = jet.colors)
  dev.off()
}
}


# perform set of measurements ----
customMeasureGAWithHybrid(seq(0, 1, 0.1), "mut", "p.mutacji", "Znalezione minimum dla różnych p. mutacji")

customGAMeasure(seq(0, 1, 0.1), "mut", 
	"p. mutacji", "Znalezione minimum dla różnych p. mutacji")
customGAMeasure(seq(10, 100, 10), "pop", 
	"rozmiar populacji", "Znalezione minimum dla różnych rozmiarów populacji")



# hybrid algorithm ----

customHybridMeasure <- function(values, mType, xlab, main) {

  averages <- c()
  for (value in values) {
    sum <- 0
    for (i in 1:nOfRuns) {
      
      message(paste("Sekwencja: ", value))
      message(paste("Przebieg: ", i))
      
      GAmin <- ga(type = "real-valued",
                  fitness =  function(xx) -f(xx),
                  min = c(B[1,]), max = c(B[2,]),
                  optim = TRUE,
                  optimArgs = list (
                    poptim = if (mType == "poptim") value else 0.05, 
                    pressel = if (mType == "pressel") value else 0.5))
      solution <- matrix(unlist(GAmin@solution),ncol=dim,byrow=TRUE)
      eval <- f(solution[1,])
      sum <- sum + eval
    }
    averages <- c(averages, (sum / nOfRuns))
  }

  if (graphs) {
    
    # save graph with measurement series to file
    png(file = paste(funcName, mType, ".png", sep=""), width=600, height=400, units="px")
    plot(0, 0, main=main,
         ylim=c(min(c(averages,globalOpt)),max(c(averages,globalOpt))),
         xlim=c(min(values),max(values)),
         type="n", xlab=xlab, ylab="wartość")
    abline(globalOpt,0, col="green")
    lines(values, averages, col = "red", type = 'l')
    legend("topright", c("memetyczny"), lwd=rep(2,1), lty=rep(1,1), col=c("red"))
    dev.off()
  }
}

customHybridMeasure(seq(0, 1, 0.05), "poptim", 
                "p. lokalnego searcha", "Znalezione minimum dla różnych poptimów")
customHybridMeasure(seq(0, 1, 0.1), "pressel", 
                    "ciśnienie", "Znalezione minimum dla różnych ciśnień")


# PSO tests ----


nOfRuns = 1 # zostaje bo niby nie można uśredniać?


#TODO
customPSOMeasure <- function(values, mType, xlab, main) {
  
  averages <- c()
  for (value in values) {
    sum <- 0
    for (i in 1:nOfRuns) {
      
      message(paste("Sekwencja: ", value))
      message(paste("Przebieg: ", i))
      
      GAmin <- ga(type = "real-valued",
                  fitness =  function(xx) -f(xx),
                  min = c(B[1,]), max = c(B[2,]),
                  optim = TRUE,
                  optimArgs = list (
                    poptim = if (mType == "poptim") value else 0.05, 
                    pressel = if (mType == "pressel") value else 0.5))
      
      solution <- matrix(unlist(GAmin@solution),ncol=dim,byrow=TRUE)
      eval <- f(solution[1,])
      sum <- sum + eval
    }
    averages <- c(averages, (sum / nOfRuns))
  }
  
  if (graphs) {
    
    # save graph with measurement series to file
    png(file = paste(funcName, mType, ".png", sep=""), width=600, height=400, units="px")
    plot(0, 0, main=main,
         ylim=c(min(c(averages,globalOpt)),max(c(averages,globalOpt))),
         xlim=c(min(values),max(values)),
         type="n", xlab=xlab, ylab="wartość")
    abline(globalOpt,0, col="green")
    lines(values, averages, col = "red", type = 'l')
    legend("topright", c("memetyczny"), lwd=rep(2,1), lty=rep(1,1), col=c("red"))
    dev.off()
  }
}


n <- 500 #ilosc czastek
m.l <- 50 #ilosc przebiegow
w <- 0.95
c1 <- 0.2
c2 <- 0.2
xmin <- c(-5.12, -5.12)
xmax <- c(5.12, 5.12)
vmax <- c(4, 4)

#inaczej są parametry podawane, trzeba zrobić dodatkowego wrappera na f()
g <- function(x) {  
  -(200 + x[,1]^2 + x[,2]^2 + cos(2*pi*x[,2]))
}

psoptim(FUN=g, n=n, max.loop=m.l, w=w, c1=c1, c2=c2,
        xmin=xmin, xmax=xmax, vmax=vmax, seed=NULL, anim=FALSE)



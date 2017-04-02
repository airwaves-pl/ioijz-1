
rm(list=ls())
dev.off(dev.list()["RStudioGD"])

require("GA")
require("globalOptTests")
require("rgl")

# Params ----

nOfRuns <- 10 # minimum 5

colors <- c("red", "blue", "orange")
series <- c("Seria 1", "Seria 2", "Seria 3")

baseParams = matrix(
  c(0, 0, 50, 100, 1,
    0.1, 0.8, 50, 100, 2,
    0.1, 0.8, 25, 50, 3), 
  # [mutations,crossovers,populations,iterations,color]*
  nrow=3, ncol=5, byrow = TRUE)

graphs <- TRUE
quality <- 100 #graph resolutions

mutationTests <- seq(0, 1, 0.1)
crossoverTests <- seq(0, 1, 0.1)
populationTests <- seq(10, 100, 5)
iterationTests <- seq(10, 200, 10)
elitismTests <- seq(0, 1, 0.1)

# Functions ----

functions <- c("Branin", "Gulf", "CosMix4", "EMichalewicz", "Hartman6", "PriceTransistor", "Schwefel", "Zeldasine20")
for (abc in functions)
{
  funcName <- abc

# for manual launch
#funcName <- "Branin" #2d
#funcName <- "Gulf" #3d
#funcName <- "CosMix4" #4d
#funcName <- "EMichalewicz" #5d
#funcName <- "Hartman6" #6d
#funcName <- "PriceTransistor" #9d
#funcName <- "Schwefel" #10d
#funcName <- "Zeldasine20" #20d

# Processing ----

customMeasure <- function(fileName, graphName, values, mType, xlab, main) {
  
  gMin <- .Machine$integer.max
  gBest <- NA
  
  temp <- c()
  for (defRow in 1:nrow(baseParams)) {
    averages <- c()
    for (value in values) {
      sum <- 0
      for (i in 1:nOfRuns) {
        GAmin <- ga(type = "real-valued",
            fitness =  function(xx) -f(xx),
            min = c(B[1,]), max = c(B[2,]),
            popSize = if (mType == "pop") value else baseParams[defRow,3],
            maxiter = if (mType == "itr") value else baseParams[defRow,4],
            pmutation = if (mType == "mut") value else baseParams[defRow,1], 
            pcrossover = if (mType == "crs") value else baseParams[defRow,2],
            elitism = if (mType == "elt") value else max(1, round(baseParams[defRow,3] * 0.05)))
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
  result <- matrix(c(temp),nrow = nrow(baseParams),ncol = length(values))
  write.table(result, file = fileName, row.names=FALSE, 
              na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    
    png(file = paste(funcName, graphName, ".png", sep=""), width=600, height=400, units="px")
    
    plot(0, 0, main=main,
         ylim=c(min(c(temp,globalOpt)),max(c(temp,globalOpt))),
         xlim=c(min(values),max(values)),
         type="n", xlab=xlab, ylab="wartosc")
    abline(globalOpt,0, col="green")
    colorNames <- c()
    seriesNames <- c()
    for (i in 1:nrow(baseParams)) {
      color <- colors[baseParams[i,5]]
      colorNames <- c(colorNames, color)
      seriesNames <- c(seriesNames, series[baseParams[i,5]])
      lines(values, result[i,], col = color, type = 'l')
    }
    legend("topright", seriesNames, lty=rep(1,nrow(baseParams)), col=colorNames)

    dev.off()
    
    summary(gBest)
    
    png(file = paste(funcName, graphName, mType, ".png", sep=""), width=600, height=400, units="px")
    filled.contour(x, y, z, color.palette = jet.colors, nlevels = 24, 
         plot.axes = { axis(1); axis(2);
           points(solution[1,1], solution[1,2], 
                  pch = 3, cex = 5, col = "black", lwd = 2) 
         }
    )
    dev.off()
    
    png(file = paste(funcName, graphName, mType, "fitness", ".png", sep=""), width=600, height=400, units="px")
    plot(gBest)
    dev.off()
  }
  
}

dim <- getProblemDimen(funcName)
B <- matrix(unlist(getDefaultBounds(funcName)),ncol=dim,byrow=TRUE)
f <- function(xx) goTest(par=c(xx, rep(0, dim-length(xx))), 
                         fnName=funcName, checkDim = TRUE)
globalOpt <- getGlobalOpt(funcName)

if (graphs) {

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

customMeasure("resultsMutations.csv", "2", mutationTests, "mut", "p. mutacji",
              "Znalezione minimum dla roznych prawdopodobienstw mutacji")

customMeasure("resultsCrossover.csv", "3", crossoverTests, "crs", "p. krzyzowania",
              "Znalezione minimum dla roznych prawdopodobienstw krzyzowania")

customMeasure("resultsPopulation.csv", "4", populationTests, "pop", "rozmiar populacji",
              "Znalezione minimum dla roznych rozmiarow populacji")

customMeasure("resultsIterations.csv", "5", iterationTests, "itr","ilosc iteracji",
              "Znalezione minimum dla roznych ilosci iteracji")

customMeasure("resultsElitism.csv", "6", elitismTests, "elt", "elityzm",
              "Znalezione minimum dla roznych wartosci elityzmu")

}
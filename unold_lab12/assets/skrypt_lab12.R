
rm(list=ls())

require("GA")
require("globalOptTests")
require("rgl")

graphs <- TRUE

# Function proxies ----
ackley <- function(xx)
{
  return(goTest(par=c(xx, rep(0, 10-length(xx))), fnName="Ackleys", checkDim = TRUE))
}

branin <- function(xx)
{
  return(goTest(par=c(xx), fnName="Branin", checkDim = TRUE))
}

schwef <- function(xx)
{
  return(goTest(par=c(xx,0,0,0,0,0,0,0,0), fnName="Schwefel", checkDim = TRUE))
}

# Params ----

isSingleTest = FALSE
n <- 7               # default 7
GAPopulation <- 50  # default 500
GAIterations <- 5   # default 50
GAMutations <- 0.1   # % (def 0.1)
GACrossovers <- 0.8  # % (def 0.8)

# Ackleys 2D function ----

# xFrom <- -32.768
# xTo <- 32.768
# f <- Vectorize(ackley)
# is3D = FALSE

# Ackleys 3D function ----

d <- 0.8
B <- matrix(c(-32.768, -32.768, 32.768, 32.768), nrow=2, ncol=2, byrow = TRUE)
f <- function(x, y) ackley(c(x,y))
is3D = TRUE

# Branin function ----

# d <- 0.5
# B <- matrix(c(-5, 0, 10, 15), nrow=2, ncol=2, byrow = TRUE)
# f <- function(x, y) branin(c(x,y))
# is3D = TRUE

# Schwefel function ----

# d <- 5
# B <- matrix(c(-500, -500, 500, 500), nrow=2, ncol=2, byrow = TRUE)
# f <- function(x, y) schwef(c(x,y))
# is3D = TRUE

# Processing ----

if (is3D) {
  x <- seq(B[1,1], B[2,1], by = d)
  y <- seq(B[1,2], B[2,2], by = d)
}

if (graphs) {
  if (is3D) {
    z <- outer(x, y, Vectorize(f))
    nbcol = 100
    color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
    zcol  = cut(z, nbcol)
    persp3d(x, y, z, theta=50, phi=25, expand=0.75, col=color[zcol],
            ticktype="detailed",axes=TRUE)
    persp3D(x, y, z, theta = -45, phi = 20, color.palette = jet.colors)
  } else {
    curve (f, xFrom, xTo)
  }
}

if (isSingleTest) {

  sum <- 0
  vector <- rep(NA,n)
  for (i in 1:n) {
    if (is3D) {
      GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x[1],x[2]), min = c(B[1,1],B[1,2]), max = c(B[2,1],B[2,2]), popSize = GAPopulation, maxiter = GAIterations, pmutation=GAMutations, pcrossover = GACrossovers)
      eval <- f(GAmin@solution[,1],GAmin@solution[,2])
    } else {
      GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x), min = xFrom, max = xTo, popSize = GAPopulation, maxiter = GAIterations, pmutation=GAMutations, pcrossover = GACrossovers)
      eval <- f(GAmin@solution[,1])
    }
    vector[i] <- eval
  }
  result <- matrix(c(vector),nrow = n,ncol = 1)
  write.table(result, file = "resultsSingle.csv", row.names=FALSE, na="", col.names=FALSE, sep=";")
  
} else {
  
  temp <- c()
  values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  averages <- c()
  for (mutation in values) {
    sum <- 0
    vector <- rep(NA,n)
    for (i in 1:n) {
      if (is3D) {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x[1],x[2]), min = c(B[1,1],B[1,2]), max = c(B[2,1],B[2,2]), popSize = GAPopulation, maxiter = GAIterations, pmutation=mutation, pcrossover = GACrossovers)
        eval <- f(GAmin@solution[,1],GAmin@solution[,2])
      } else {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x), min = xFrom, max = xTo, popSize = GAPopulation, maxiter = GAIterations, pmutation=mutation, pcrossover = GACrossovers)
        eval <- f(GAmin@solution[,1])
      }
      sum <- sum + eval
      vector[i] <- eval
    }
    temp <- c(temp, vector)
    avg <- sum / n
    averages <- c(averages, avg)
  }
  result <- matrix(c(temp),nrow = n,ncol = length(values))
  write.table(result, file = "resultsMutations.csv", row.names=FALSE, na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    plot(values, averages, 
         main="Wartoœæ funkcji celu dla ró¿nych wartoœci P. mutacji", 
         ylim=c(-0.05,max(averages)), 
         type="l", col="red", xlab="parametry", ylab="wartoœæ funkcji")
    abline(0,0, col="green")
  }
  
  temp <- c()
  values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  averages <- c()
  for (crossover in values) {
    sum <- 0
    vector <- rep(NA,n)
    for (i in 1:n) {
      if (is3D) {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x[1],x[2]), min = c(B[1,1],B[1,2]), max = c(B[2,1],B[2,2]), popSize = GAPopulation, maxiter = GAIterations, pmutation=GAMutations, pcrossover = crossover)
        eval <- f(GAmin@solution[,1],GAmin@solution[,2])
      } else {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x), min = xFrom, max = xTo, popSize = GAPopulation, maxiter = GAIterations, pmutation=GAMutations, pcrossover = crossover)
        eval <- f(GAmin@solution[,1])
      }
      sum <- sum + eval
      vector[i] <- eval
    }
    temp <- c(temp, vector)
    avg <- sum / n
    averages <- c(averages, avg)
  }
  result <- matrix(c(temp),nrow = n,ncol = length(values))
  write.table(result, file = "resultsCrossover.csv", row.names=FALSE, na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    plot(values, averages, 
         main="Wartoœæ funkcji celu dla ró¿nych wartoœci P. krzy¿owania", 
         ylim=c(-0.05,max(averages)), 
         type="l", col="red", xlab="parametry", ylab="wartoœæ funkcji")
    abline(0,0, col="green")
  }
  
  temp <- c()
  values <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  averages <- c()
  for (elitism in values) {
    sum <- 0
    vector <- rep(NA,n)
    for (i in 1:n) {
      if (is3D) {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x[1],x[2]), min = c(B[1,1],B[1,2]), max = c(B[2,1],B[2,2]), popSize = GAPopulation, maxiter = GAIterations, pmutation=GAMutations, pcrossover = GACrossovers, elitism = elitism)
        eval <- f(GAmin@solution[,1],GAmin@solution[,2])
      } else {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x), min = xFrom, max = xTo, popSize = GAPopulation, maxiter = GAIterations, pmutation=GAMutations, pcrossover = GACrossovers, elitism = elitism)
        eval <- f(GAmin@solution[,1])
      }
      sum <- sum + eval
      vector[i] <- eval
    }
    temp <- c(temp, vector)
    avg <- sum / n
    averages <- c(averages, avg)
  }
  result <- matrix(c(temp),nrow = n,ncol = length(values))
  write.table(result, file = "resultsElitism.csv", row.names=FALSE, na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    plot(values, averages, 
         main="Wartoœæ funkcji celu dla ró¿nych wartoœci elityzmu", 
         ylim=c(-0.05,max(averages)), 
         type="l", col="red", xlab="parametry", ylab="wartoœæ funkcji")
    abline(0,0, col="green")
  }
  
  temp <- c()
  values <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
  averages <- c()
  for (population in values) {
    sum <- 0
    vector <- rep(NA,n)
    for (i in 1:n) {
      if (is3D) {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x[1],x[2]), min = c(B[1,1],B[1,2]), max = c(B[2,1],B[2,2]), popSize = population, maxiter = GAIterations, pmutation=GAMutations, pcrossover = GACrossovers)
        eval <- f(GAmin@solution[,1],GAmin@solution[,2])
      } else {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x), min = xFrom, max = xTo, popSize = population, maxiter = GAIterations, pmutation=GAMutations, pcrossover = GACrossovers)
        eval <- f(GAmin@solution[,1])
      }
      sum <- sum + eval
      vector[i] <- eval
    }
    temp <- c(temp, vector)
    avg <- sum / n
    averages <- c(averages, avg)
  }
  result <- matrix(c(temp),nrow = n,ncol = length(values))
  write.table(result, file = "resultsPopulation.csv", row.names=FALSE, na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    plot(values, averages, 
         main="Wartoœæ funkcji celu dla ró¿nych rozmiarów populacji", 
         ylim=c(-0.05,max(averages)), 
         type="l", col="red", xlab="parametry", ylab="wartoœæ funkcji")
    abline(0,0, col="green")
  }
  
  temp <- c()
  values <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  averages <- c()
  for (iterations in values) {
    sum <- 0
    vector <- rep(NA,n)
    for (i in 1:n) {
      if (is3D) {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x[1],x[2]), min = c(B[1,1],B[1,2]), max = c(B[2,1],B[2,2]), popSize = GAPopulation, maxiter = iterations, pmutation=GAMutations, pcrossover = GACrossovers)
        eval <- f(GAmin@solution[,1],GAmin@solution[,2])
      } else {
        GAmin <- ga(type = "real-valued", fitness =  function(x) -f(x), min = xFrom, max = xTo, popSize = GAPopulation, maxiter = iterations, pmutation=GAMutations, pcrossover = GACrossovers)
        eval <- f(GAmin@solution[,1])
      }
      sum <- sum + eval
      vector[i] <- eval
    }
    temp <- c(temp, vector)
    avg <- sum / n
    averages <- c(averages, avg)
  }
  result <- matrix(c(temp),nrow = n,ncol = 10)
  write.table(result, file = "resultsIterations.csv", row.names=FALSE, na="", col.names=FALSE, sep=";")
  
  if (graphs) {
    plot(values, averages, 
         main="Wartoœæ funkcji celu dla ró¿nych iloœci iteracji", 
         ylim=c(-0.05,max(averages)), 
         type="l", col="red", xlab="parametry", ylab="wartoœæ funkcji")
    abline(0,0, col="green")
  }
  
}

if (graphs) {
  summary(GAmin)
  if (is3D) {
    filled.contour(x, y, z, color.palette = jet.colors, nlevels = 24, 
       plot.axes = { 
         axis(1); 
         axis(2);
         points(GAmin@solution[,1], GAmin@solution[,2], pch = 3, cex = 5, col = "black", lwd = 2) 
       }
    )
  } else {
    points(GAmin@solution, f(GAmin@solution), col="red", pch=19)
    abline(v=GAmin@solution, h=f(GAmin@solution), col="red", lty=2)
  }
  plot(GAmin)
}

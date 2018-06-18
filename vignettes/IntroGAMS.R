#---- Introducing GAMS

library(gamair)
library(ggplot2)
data("engine")
attach(engine)
x <- size - min(size); x <- x / max(x)
detach(engine)
engine$capacity <- x
p <- ggplot(engine, aes(capacity, wear)) + geom_point()

rk <- function(x, z) {
  #R(x,z) for cubic spline on [0,1]
  ( (z - 0.5)^2 - 1/12 ) *  (( x - 0.5)^2 - 1/12)/4 - ((abs(x-z) - 0.5)^4 - (abs(x-z)-0.5)^2/2+7/240)/24
}

spl.X <- function(x, xk){
  #set u-p model matrix for cubic penalised regression spline
  q <- length(xk) + 2 #number of parameters
  n <- length(x) # number of data
  X <- matrix(1,n,q) # initialised model matrix
  X[, 2] <- x #set second column to x
  X[, 3:q] <- outer(x, xk, FUN = rk)
  X
}

xk <- 1:4/5 #chose some knots
X <- spl.X(engine$capacity, xk) #generate model matrix
mod.1 <- lm(engine$wear~X-1) #fit model
xp <- 0:100/100 #x calues for prediction
Xp <- spl.X(xp, xk) #prediction matrix
data.frame( x = xp,
            y = Xp %*% coef(mod.1) )
p + geom_line(data = data.frame( x = xp,
                                 y = Xp %*% coef(mod.1) ), aes(x, y))


#------ GAMMs with R
data("sole")


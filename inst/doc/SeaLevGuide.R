## ----include=FALSE-------------------------------------------------------
library(knitr)

opts_chunk$set(cache = FALSE,
               fig.path = 'Rgraphics/fig', dev = 'pdf', tidy = FALSE,
               fig.width = 5.4, fig.height = 5.4)
library(SeaLev)
data(Brest.tide)

## ----options, echo = FALSE-----------------------------------------------
options(prompt = "R> ", continue = "   ", width = 80)
RVersion <- R.Version()$version.string

## ----label=BrestTide, fig.keep="last", fig.show="hide"-------------------
library(SeaLev)
data(Brest.tide)
class(Brest.tide)
str(Brest.tide)
plot(Brest.tide, col = "SeaGreen", type = "l",
     main = "Density of high-tide sea level in Brest")
grid(); abline(h = 0)

## ----label=BrestTide0----------------------------------------------------
Brest.tide$y[c(1L, length(Brest.tide$y))]
h <- diff(Brest.tide$x)[1]                     
h * sum(Brest.tide$y) 

## ----label=BrestSurgePar1------------------------------------------------
u <- 50 
theta.y <- c("scale" = 10, "shape" = -0.01)
lambda <- 1.6

## ----label=BrestConv0, fig.show="hide", warning=FALSE--------------------
conv.gpd0 <- convSL(dens.x = Brest.tide,
                   threshold.y = u, distname.y = "GPD",
                   lambda = lambda, par.y = theta.y,
                   main = "Sea-level with GPD surges: given parameters")

## ----label=BrestFitSurge, fig.show="hide"--------------------------------
library(Renext); data(Brest) 
fit.gpd1 <- Renouv(x = Brest$OTdata$Surge,
                   effDuration = as.numeric(Brest$OTinfo$effDuration),
                   threshold = 50, distname.y = "GPD",
                   main = "GPD surge")
coef(fit.gpd1)

## ----label=BrestSL-------------------------------------------------------
cov1 <- vcov(fit.gpd1)
cov1

## ----label=BrestConv1, include=FALSE, warning=FALSE----------------------
conv.gpd1 <- convSL(dens.x = Brest.tide,
                    threshold.y = 50,
                    distname.y = "GPD",
                    lambda = lambda, par.y = theta.y,
                    covpar.y = cov1,
                    main = "Sea-level for Brest with GPD surges")

## ----label=BrestConv15, fig.keep="all", fig.show="hide", warning=FALSE----
conv.gpd2a <- convSL(dens.x = Brest.tide,
                     threshold.y = 50,
                     distname.y = "GPD",
                     lambda = lambda, par.y = theta.y,
                     pct.conf = c(95, 90),
                     filled.conf = TRUE, mono = FALSE,
                     covpar.y = cov1,
                     main = "Sea-level for Brest with GPD surges (lambda known)")

## ----label=BrestSL1------------------------------------------------------
cov1[-1, -1]

## ----label=BrestConv2, fig.show="hide", fig.keep="all", warning=FALSE----
conv.gpd2 <- convSL(dens.x = Brest.tide,
                    threshold.y = 50,
                    distname.y = "GPD",
                    lambda = lambda, par.y = theta.y,
                    use.covlambda = FALSE,
                    pct.conf = c(95, 90),
                    filled.conf = TRUE, mono = FALSE,
                    covpar.y = cov1,
                    main = "Sea-level for Brest with GPD surges (lambda known)")

## ----label=BrestConvGEV, fig.show="hide", fig.keep="all", warning=FALSE----
par.y <- c(loc = -10.8, scale = 10)
res.gumbel <- convSL(dens.x = Brest.tide,
                     threshold.y = NA,           
                     distname.y = "gumbel",
                     lambda = 705.8,
                     par.y = par.y,
                     filled.conf = TRUE, mono = FALSE,
                     main = "Gumbel surges")

## ----label=BrestConv25---------------------------------------------------
head(conv.gpd2$pred, n = 3)

## ----label=BrestConv3, fig.show="hide", fig.keep="all", warning=FALSE----
conv.gpd3 <- convSL(dens.x = Brest.tide,
                    threshold.y = 50, distname.y = "GPD",
                    lambda = lambda, par.y = theta.y, covpar.y = cov1,
                    ylim = c(300, 600),
                    main = "Sea-level for Brest with GPD surges (lambda known)")

## ----label=BrestConv4, fig.show="hide", fig.keep="all", warning=FALSE----
conv.gpd3 <- convSL(dens.x = Brest.tide,
                    threshold.y = 50, distname.y = "GPD",
                    lambda = lambda, par.y = theta.y, covpar.y = cov1,
                    Tlim = c(100, 3000),
                    main = "Sea-level for Brest with GPD surges (lambda known)")

## ----empirical, fig.show="hide", fig.keep="all", warning=FALSE-----------
res.g2 <- convSL(dens.x = Brest.tide,
                 threshold.y = NA, distname.y = "gumbel",
                 lambda = 705.8, par.y = par.y,
                 filled.conf = TRUE, mono = FALSE,
                 main = "Artificial empirical points (1 set)",
                 z = c(500, 490, 480, 460),
                 duration = 200)

## ----empirical1, fig.show="hide", fig.keep="all", warning=FALSE----------
res.g3 <- convSL(dens.x = Brest.tide,
                 threshold.y = NA, distname.y = "gumbel",
                 lambda = 705.8, par.y = par.y,
                 filled.conf = TRUE, mono = FALSE,
                 main = "Artificial empirical points (2 sets)",
                 z = list(c(500, 490, 480), c(440, 420, 380, 350)),
                 duration = c(200, 170))

## ----SplineDens1, fig.show = "hide", fig.keep="all", warning=FALSE-------
SD <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y)
SD10 <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y, nKnots = 10)
plot(SD, main = "default")
plot(SD10, main = "10 knots")

## ----SplineDens2, fig.show = "hide", fig.keep="all", warning=FALSE-------
set.seed(1234)
u <- 50
par.y <- c("scale" = rgamma(1, shape = 2, scale = 30),
           "shape" = 0.2 * runif(1))
res <- GPtail(x = SD, par.y = par.y, threshold.y = u, lambda = 1)
class(res)
plot(res)
plot(res, which = 3)

## ----SplineDens2a, fig.show = "hide", fig.keep="all", warning=FALSE------
pred <- predict(SD, newdata = seq(from = 200, to = 300, by = 0.1))
names(pred)

## ----SplineDens3, fig.show = "hide", fig.keep="all", warning=FALSE-------
set.seed(1234)
SDrand <- rSplineDensity(order = 4, xmax = 10)
plot(SDrand, main = "random spline density")

## ----label=FAQ1----------------------------------------------------------
dMydist <- function(x, bar) foodens(x, bar = bar)

## ----label=BrestAsympt1, fig.show="hide", fig.keep="all", warning=FALSE----
theta2.y <- c("rate" = 0.10)
conv.asympt <- convSL(dens.x = Brest.tide,
                      threshold.y = 50, 
                      distname.y = "exponential",
                      lambda = lambda, 
                      par.y = theta2.y, 
                      show.asympt = TRUE,
                      Tlim = c(5, 1000000),
                      main = "Asymptotic curve: exponential Y")

## ----label=BrestAsympt2, fig.show="hide", fig.keep="all", warning=FALSE----
theta3.y <- c("scale" = 10, "shape" = 0.03)
conv.asympt <- convSL(dens.x = Brest.tide,
                      threshold.y = 50, 
                      distname.y = "GPD",
                      lambda = lambda, 
                      par.y = theta3.y, 
                      show.asympt = TRUE,
                      Tlim = c(5, 1000000),
                      main = sprintf("Asymptotic curve: GPD Y with shape %4.2f", 
                          theta3.y["shape"]))


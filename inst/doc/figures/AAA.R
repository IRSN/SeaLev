library(SeaLev)
data(Brest.tide)
library(RColorBrewer)
convSL <- convSL2
data(Brest.tide)
fit.gpd1 <- Renouv(x = Brest$OTdata$Surge,
                   effDuration = as.numeric(Brest$OTinfo$effDuration),
                   threshold = 50, distname.y = "GPD",
                   plot = FALSE)
cov1 <- vcov(fit.gpd1)
u <- 50 
theta.y <- coef(fit.gpd1)[-1]
lambda <- coef(fit.gpd1)[1]

theta.y["shape"] <- 0.02
par(mar = c(5, 10, 5, 5), col.axis = "darkgray", las = 0,
    col.lab = "darkgray")

set.seed(123)
w <- 200
nz <- rpois(1, lambda = lambda * w)
F <- cumsum(Brest.tide$y * diff(Brest.tide$x[1:2]))
x <- approx(x = F, y = Brest.tide$x,  xout = runif(nz))$y
## x <- Brest.tide$x[1] + (Brest.tide$x[length(Brest.tide$x)] - Brest.tide$x[1]) * x
y <-  rGPD(nz, loc = u, shape = theta.y["shape"], scale = theta.y["scale"]) 
z <- x  + y

res <- convSL(dens.x = Brest.tide,
               threshold.y = u,
               distname.y = "GPD",
               lambda = lambda, par.y = theta.y, 
               covpar.y = cov1,
               pct.conf = 95,
               ylim = c(320, 550),
               Tlim = c(10, 4000),
               main = "",
               rl.mark = NULL,
               z = z,
               duration = w,
               pch.points = 16,
               col.points = translude("SeaGreen", alpha = 0.6),
               plot = TRUE,
               grid = FALSE,
               show.asympt = TRUE,
               legend = FALSE)

## lines(log(res$Tz), res$z, log = "x", lwd = 3, type = "l", col = "orangered")

T0 <- 1500
z0 <- approx(x = res$Tz, y = res$z, xout = T0)$y
x0 <- approx(x = res$Tz, y = res$condExp.x, xout = T0)$y

abline(v = log(T0), lty = "dashed", lwd = 1)
abline(h = c(x0, z0), lty = "dashed", lwd = 1)


mtext("$\\textrm{E}[X \\vert Z = z(T)]$", las = 1, line = 2,
     side = 2, at = x0, col = "black", cex = 1.5)

mtext("$z(T)$", line = 2,
      las = 1, side = 2, 
      at = z0, col = "black", cex = 1.5)

mtext("$T$", at = log(T0), side = 1, line = 1.5, cex = 1.5)

SD <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y, order = 2, nKnots = 24)
res2 <- GPtail(x = SD, threshold.y = u, par.y = theta.y, lambda = lambda)


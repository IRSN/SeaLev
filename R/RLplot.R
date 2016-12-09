##================================================================
## Author : Y. Deville
##
## Return level plot using provided data (no computations)
##
##================================================================

RSLplot <- function(data,
                    z = NULL,
                    duration = NULL,
                    lambda,
                    xlab = "period",
                    ylab = "level",
                    mono  = TRUE,
                    lty.quant = "solid",
                    col.quant = ifelse(mono, "black", "SteelBlue"),
                    pct.conf = c(95, 70),
                    col.conf = NULL,
                    filled.conf = FALSE,
                    fill.conf = NULL,
                    lty.conf = ifelse(rep(filled.conf, 2), rep(NA, 2), c("dashed", "dotted")),
                    rl.mark = 100,
                    text.mark = rl.mark,
                    col.mark = NULL,
                    Tlim = NULL,
                    problim = NULL,
                    below = NA,
                    alpha.below = 0.5,
                    grid = TRUE,
                    legend = TRUE,
                    pch.points = c(24, 22, 23),
                    col.points = c("black", "darkgray"),
                    ...) {
    
  mc <- match.call(expand.dots = TRUE)
  nmc <- names(substitute(list(...)))

  if (filled.conf && is.unsorted(rev(pct.conf)))
    warning("'pct.conf' not in descending order. Some conf. bands may not be visible")
  
  if (mono) {
    if (is.null(col.conf)) col.conf <- c("black", "black", "black")
    ## if (is.null(col.points)) col.points <- c("black", "dargkgray")
  } else {
    if (is.null(col.conf)) {
        col.conf <- c("SteelBlue4", "orangered", "SpringGreen3", "purple",
                      "firebrick")
    }
    if (is.null(col.points)) {
        col.points <- c("red3", "SpringGreen3", "SteelBlue3")
    }
  }

  ## pch.points <- c(24, 22, 23)

  col.conf <- rep(col.conf, length.out = length(pct.conf))
  lty.conf <- rep(lty.conf, length.out = length(pct.conf))
  
  if (filled.conf) {
    if (is.null(fill.conf)) {
      fill.conf <- rep(NA, length(pct.conf))
      for (i in 1L:length(fill.conf)) {
        w.fill <- 0.5 * i / (length(pct.conf) + 1L)
        ## cat("w.fill = ", w.fill, "\n")
        rgbc <- (col2rgb(col.conf[i]) * w.fill +  col2rgb("white") * (1 - w.fill)) / 256
        rgbc <- rgb(rgbc[1L], rgbc[2L], rgbc[3L])
        fill.conf[i] <- rgbc
      }
    } else {
      fill.conf <- rep(fill.conf, length.out = length(pct.conf))
    }
  }

  cnames <- colnames(data)
  candLnames <- match(paste("L", pct.conf, sep = "."), cnames)
  candUnames <- match(paste("U", pct.conf, sep = "."), cnames)
  
  ## freq.g = inverse return period (= freqency...)
  freq.g <- lambda * (1 - data$prob)
  x.g <- -log(freq.g)

  labs <- c(1, 2,  5, 10, 20,  50, 100, 200, 500, 1000,
            2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)

  ## change of 2010-09-10
  ## colnames checks: confidence limits?
  
  ind <- grep("quant", colnames(data))
  if (!length(ind)) stop("colnames(data) must contain \"quant\"")

  ind1 <- grep("([L,U])[[:punct:]][[:digit:]]{1,2}", colnames(data))
  ind <- union(ind, ind1) 

  ry <- range(data[ , ind], na.rm = TRUE)
  
  ## if ( ("L.95" %in% colnames(data)) &&
  ##     ("U.95" %in% colnames(data)) ) 
  ##   ry <- c(min(data$L.95, na.rm = TRUE),
  ##           max(data$U.95, na.rm = FALSE))
  ## else
  ##   ry <- range(data$quant, na.rm = TRUE)
  
  if (!is.null(problim)) {
      if (!is.numeric(problim) || length(problim) != 2L ||
          any(is.na(problim)) ||
          any(problim <= 0.0) || any(problim >= 1.0) ) stop("invalid limits in 'problim'.")
      if (!is.null(Tlim)) stop("only one of 'problim' and 'Tlim' can be provided")
      rx <-  -log(lambda * c(1.0 - problim))
  } else {
      if (!is.null(Tlim)) {
          if ( !is.numeric(Tlim) || length(Tlim) != 2L ||
              any(is.na(Tlim)) ||
              any(Tlim <= 0.0) ) stop("invalid limits in 'Tlim'.")
          rx <-  log(Tlim)
      } else rx <- range(x.g)
  }
  
  ## NOT USED YET
  if (!is.null(z)) {
      ry2 <- range(z, na.rm = TRUE)
      if (ry2[1L] < ry[1L]) ry[1L] <- ry2[1L]
      if (ry2[2L] > ry[2L]) ry[2L] <- ry2[2L]
  }
  
  ## prepare plot
  plot(x = c(x.g[1], x.g, x.g[length(x.g)]),
       y = c(ry[1L], data$quant, ry[2L]),
       type = "n",
       xlab = xlab, ylab = ylab,
       xlim = rx,
       xaxt = "n",
       ...)
  
  ind <- !is.na(candLnames) & !is.na(candUnames)
  
  for (i in 1:length(ind)) {
      if (all(is.na(data[ , candLnames[i]])) &&
          all(is.na(data[ , candUnames[i]]))) {
          warning("All values missing in ", cnames[candLnames[i]], " and ",
                  cnames[candUnames[i]])
          ind[i] <- FALSE
      }
  }
  
  for (i in 1L:length(pct.conf))  {
      
    if (ind[i]) {
      
      iL <- candLnames[i]
      iU <- candUnames[i]

      if (filled.conf) {
          polygon(x = c(x.g, rev(x.g)),
                  y = c(data[ , iL], rev(data[ ,iU])),
                  col = fill.conf[i], border = NA)
      }
      
      if (!is.na(lty.conf[i])) {
          lines(x = x.g, y = data[ , iL],
                type = "l", lwd = 2, col = col.conf[i], lty = lty.conf[i])      
          lines(x = x.g, y = data[ , iU],
                type = "l", lwd = 2, col = col.conf[i], lty = lty.conf[i])
      }
        
    } else warning("confidence limits for level ", pct.conf[i], "% not found in data")
    
  }
  
  lines(x = x.g,
        y = data$quant,
        type = "l", lwd = 2, col = col.quant[1])
  
  ## abline(h = threshold, col = "SeaGreen3", lwd = 2)
  
  ## Added in version 0.2-6
  ## semi-tranparent rectangle under the level given in 'below'
  
  if (!is.na(below)) {
      
      colMask <- rgb(red = 1, green = 1, blue = 1, alpha = alpha.below)
      
      coords <- par()$usr
      
      rect(xleft = coords[1],
           xright = coords[2],
           ybottom = coords[3],
           ytop = below,
           border = NA,
           col = colMask)
      
  }
  
  ##================================
  ## add empirical z levels, if any.
  ##================================
  
  if (!is.null(z)) {
      
      if (is.list(z)) {
          if (length(duration) != length(z)) {
              stop("When 'z' is list, 'duration' must be a list or vector with the same length")
          }
          rz <- unlist(lapply(z, length))
          lapply(z, sort, decreasing = TRUE)
          col.points <- rep(col.points, length.out = length(z))
          pch.points <- rep(pch.points, length.out = length(z))
      } else {
          rz <- length(z)
          z <- list(sort(z, decreasing = TRUE))
      }
      
      Tz <- list()
      
      for (iz in 1:length(z)) {
          
          Tz[[iz]]  <- (duration[iz] * 705.8  + 1.0) / (1L:rz[iz]) / 705.8
          
          points(x = log(Tz[[iz]]), y = z[[iz]], 
                 pch = pch.points[iz], col = col.points[iz])
          
      }
      
  }
  
  
  ## add x-axis
  axis(side = 1, at = log(labs), labels = labs)
  
  ## add upper x axis ??
  ## axis(side = 3, at = x.g, labels = data$prob)

  if (grid) {
      abline(v = log(labs), col = "gray", lty = "dotted")
      abline(h = pretty(par()$usr[3L:4L]), col = "gray", lty = "dotted")
  }
      
  ## Not very useful...
  ##
  ## if (length(spec.x)) {
  ##   spec.f <- log(spec.rl)
  ##  points(x = spec.f, y = spec.x,
  ##         pch = spec.pch, bg = spec.bg, cex = spec.cex, col = spec.col)
  ## }
  
  Cols <- c("SteelBlue3", "orangered")

  if (legend) {
      if (sum(ind)) {
          if (filled.conf) {
              legend("topleft",
                     legend = paste(pct.conf[ind], "%", sep = ""),
                     col = col.conf[ind],
                     lwd = 2,
                     lty = lty.conf[ind],
                     fill = fill.conf)
          } else {
              legend("topleft",
                     legend = c("theo.", paste(pct.conf[ind], "%", sep = "")),
                     col = c(col.quant[1L], col.conf[ind]),
                     lwd = 2,
                     lty = c(lty.quant[1L], lty.conf[ind]))
          }
      } else {
          legend("topleft",
                 legend = "theo.",
                 col = col.quant[1L],
                 lwd = 2,
                 lty = lty.quant[1L])
      }
  }
      
  if (length(rl.mark)) {
      
      if (mono || (is.null(col.mark))) col.mark <- rep("black", length(rl.mark))
      else col.mark <- rep(col.mark, length.out = length(rl.mark))
      
      for (i in 1L:length(rl.mark)) {
          abline(v = log(rl.mark[i]), col = col.mark[i])
          text(x = log(rl.mark[i]),
               y = ry[1] + (ry[2] - ry[1]) / 6,
               labels = text.mark[i],
               col = col.mark[i],
               pos = 4L)
      }
      
  }
  
}


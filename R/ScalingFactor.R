kurtosis <- function(x, na.rm = FALSE) {
  # kurtosis calculation function from an R package "moments"
  if (is.matrix(x)) {
    apply(x, 2, kurtosis, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    n * sum((x - mean(x))^4) / (sum((x - mean(x))^2)^2)
  }
  else if (is.data.frame(x)) {
    sapply(x, kurtosis, na.rm = na.rm)
  } else {
    kurtosis(as.vector(x), na.rm = na.rm)
  }
}
ndotproduct <- function(x, y, m = 0L, n = 1, na.rm = TRUE, ...) {
  .weightxy <- function(x, y, m = 0, n = 1) {
    x ^ m * y ^ n
  }
  wx <- .weightxy(x[, 1L], x[, 2L], m, n)
  wy <- .weightxy(y[, 1L], y[, 2L], m, n)
  sum(wx * wy, na.rm = na.rm)^2L /
    (sum(wx^2L, na.rm = na.rm) * sum(wy^2L, na.rm = na.rm))
}
opt.weight <- function(rdata, wmz = wmzVector, wint = wintVector, plot = TRUE, fold = 1, verbose = FALSE, nseed = 1, R = 100) {
  # rdata # it should be intensity data with m/z values (row) by compounds (column)
  # wmz = wmzVector # mz weight factors
  # wint = wintVector # int weight factors
  # plot = TRUE # the scatter plot between weight factors and the conditional means of skewness/kurtosis
  # fold = 1 # the size of the sub-data set = the total size of a library / fold
  # verbose = FALSE # information on the evolution of the iterative algorithm is printed
  # nseed = 1 # necessary only when fold > 1
  # R = 100 # necessary only when fold > 1
  # optimal weight factors
  # calculate the conditional means and find the optimal weight factors
  wmz <- sort(unique(wmz))
  wint <- sort(unique(wint))
  lenmz <- length(wmz) # total number of m/z weight factors
  lenint <- length(wint) # total number of intensity weight factors
  
  # random selection
  if (fold > 1) {
    skdata <- skew.kurt.random(rdata = rdata, fold = fold, nseed = nseed, R = R, wmz = wmz, wint = wint, verbose = verbose)
  } else { # using the whole reference library
    skdata <- skew.kurt.norandom(rdata = rdata, wmz = wmz, wint = wint, verbose = verbose)
  }
  sks <- skdata$stat
  
  if (verbose) {
    cat("\n#*** FINDING OPTIMAL WEIGHT FACTORS ***#\n")
  }
  
  # optimal weight factor for intensity
  int.opt <- c()
  for (i in 1:lenint) {
    td <- sks[sks$int == wint[i], ]
    int.opt <- rbind(int.opt, c(wint[i], mean(td$sk), sd(td$sk), mean(td$skew), sd(td$skew), mean(td$kurt), sd(td$kurt)))
  }
  int.opt <- data.frame(int.opt)
  dimnames(int.opt)[[2]] <- c("int", "sk", "sdsk", "skew", "sdskew", "kurt", "sdkurt")
  int.opt <- int.opt[order(int.opt$sk, decreasing = T), ]
  ow.int.v <- int.opt[int.opt$sk == max(int.opt$sk)[1], ]
  
  # optimal weight factor for m/z value
  mz.opt <- c()
  for (i in 1:lenmz) {
    td <- sks[sks$mz == wmz[i], ]
    mz.opt <- rbind(mz.opt, c(wmz[i], mean(td$sk), sd(td$sk), mean(td$skew), sd(td$skew), mean(td$kurt), sd(td$kurt)))
  }
  mz.opt <- data.frame(mz.opt)
  dimnames(mz.opt)[[2]] <- c("mz", "sk", "sdsk", "skew", "sdskew", "kurt", "sdkurt")
  mz.opt <- mz.opt[order(mz.opt$sk, decreasing = T), ]
  ow.mz.v <- mz.opt[mz.opt$sk == max(mz.opt$sk)[1], ]
  
  if (plot) {
    # dev.new()
    svg(paste(Name, "Cutoff", Cutoff, " int", ow.int.v[1], "mz", ow.mz.v[1], ".svg"))
    xrange <- range(int.opt[, 1], mz.opt[, 1])
    yrange <- range(int.opt[, 2], mz.opt[, 2])
    plot(int.opt[order(int.opt$int), 1], int.opt[order(int.opt$int), 2], xlab = "Weight factors", ylab = "E(skewness/kurtosis)", type = "o", col = 2, pch = 19, lwd = 2)
    points(mz.opt[order(mz.opt$mz), 1], mz.opt[order(mz.opt$mz), 2], type = "o", col = 3, pch = 19, lwd = 2)
    legend("topright", c("intensity", "m/z"), pch = 19, col = c(2, 3), lwd = 2, lty = 1, bty = "n")
    points(ow.int.v[1], ow.int.v[2], pch = 19, col = 5)
    points(ow.mz.v[1], ow.mz.v[2], pch = 19, col = 4)
    legend("right", c(paste("Optimal intensity factor = ", ow.int.v[1], sep = ""), paste("Optimal m/z factor = ", ow.mz.v[1], sep = "")),
           bty = "n",
           pch = 19, col = c(5, 4)
    )
    title("The optimal weight factors")
    dev.off()
  }
  
  ow.int <- as.numeric(ow.int.v[1])
  ow.mz <- as.numeric(ow.mz.v[1])
  
  if (verbose) {
    cat("# The optimal intensity weight factor = ", ow.int, "\n")
    cat("# The optimal m/z weight factor = ", ow.mz, "\n")
    cat("# The optimal factors up to the third rank \n")
    cat("# * intensity * = ", int.opt$int[int.opt$sk >= sort(int.opt$sk)[5]], "\n")
    cat("# *    m/z    * = ", mz.opt$mz[mz.opt$sk >= sort(mz.opt$sk)[5]], "\n\n")
  }
  
  list(
    opt.int = ow.int,
    opt.mz = ow.mz,
    int.sk = int.opt,
    mz.sk = mz.opt,
    skdata = skdata,
    opt = list(int = ow.int.v, mz = ow.mz.v)
  )
}

skew.kurt.norandom <- function(rdata , wmz = wmzVector,wint = wintVector,verbose = FALSE) {
  #rdata: it should be intensity data with m/z values (row) by compounds (column)
  # when the size of a reference library is tractable
  # calculate skewness and kurtosis using the whole reference library
  wmz <- sort(unique(wmz))
  wint <- sort(unique(wint))
  lenmz <- length(wmz) # total number of m/z weight factors
  lenint <- length(wint) # total number of intensity weight factors
  lenwt <- lenmz * lenint # total size of results
  
  rlt <- skew.kurt.id(rdata = rdata, id = 1, fold = 1, nseed = -1, wmz = wmz, wint = wint, verbose = verbose)
  if (dim(rlt)[1] != lenwt) {
    cat("*** ERROR: it (", dim(rlt)[1], ") has fewer number of rows than ", lenwt, " ***\n")
    stop()
  }
  if (verbose) {
    cat("# it is done\n")
  }
  
  list(stat = rlt, mz = wmz, int = wint)
}

skew.kurt.id <- function(rdata, id = 1,fold = 10,nseed = 1,wmz = wmzVector,wint = wintVector,verbose = FALSE) {
  # it should be intensity data with m/z values (row) by compounds (column)
  # skewness and kurtosis calculation
  mz <- 1:dim(rdata)[1]
  rdatasize <- dim(rdata)[2]
  newsize <- round(rdatasize / fold)
  
  # select a random subset
  if (!is.na(nseed)) {
    set.seed(nseed)
  }
  
  if (nseed != -1) {
    pos <- sample(1:rdatasize, newsize)
    rdata <- rdata[, pos]
    if (verbose) {
      cat("# id=", id, ": the total number of compounds selected = ", length(pos), "\n")
    }
  }
  
  rlt <- c()
  for (i in 1:length(wmz)) {
    for (j in 1:length(wint)) {
      if (verbose) {
        if (fold > 1) {
          cat("# current weights for ", id, "-th RUN: (int,mz) = (", wint[j], ",", wmz[i], ")\n")
        } else {
          cat("# current weights: (int,mz) = (", wint[j], ",", wmz[i], ")\n")
        }
      }
      # i <- 1
      # j <- 1
      wrdata <- t((rdata^wint[j]) * (mz^wmz[i]))
      wrdata <- wrdata / sqrt(apply(wrdata^2, 1, sum))
      
      n.cor <- wrdata %*% t(wrdata)
      rm(wrdata)
      n.cor <- as.numeric(n.cor[upper.tri(n.cor)])
      
      tskew <- skewness(n.cor)
      tkurt <- kurtosis(n.cor)
      tsk <- tskew / tkurt
      
      trlt <- as.numeric(c(wint[j], wmz[i], id, tskew, tkurt, tsk))
      rlt <- rbind(rlt, trlt)
      rm(n.cor)
    }
  }
  rdim <- dim(rlt)
  rlt <- data.frame(matrix(c(as.numeric(rlt)), rdim[1], rdim[2]))
  dimnames(rlt)[[2]] <- c("int", "mz", "id", "skew", "kurt", "sk")
  
  if (nseed == -1) {
    rlt <- rlt[, -3]
  }
  
  rlt
}


skew.kurt.random <- function(rdata , fold = 10,nseed = 1,R = 100 , wmz = wmzVector,wint = wintVector,verbose = FALSE) {
  # when the size of a reference library is not tractable #
  # calculate skewness and kurtosis with random selection
  # rdata: it should be intensity data with m/z values (row) by compounds (column)
  # R： total number of sub-data generated randomly
  wmz <- sort(unique(wmz))
  wint <- sort(unique(wint))
  lenmz <- length(wmz) # total number of m/z weight factors
  lenint <- length(wint) # total number of intensity weight factors
  lenwt <- lenmz * lenint # total size of results
  
  rlt <- c()
  for (i in 1:R) {
    trlt <- skew.kurt.id(rdata = rdata, id = i, fold = fold, nseed = nseed, wmz = wmz, wint = wint, verbose = verbose)
    if (dim(trlt)[1] != lenwt) {
      cat("*** ERROR in ", i, "-th RUN: it (", dim(trlt)[1], ") has fewer number of rows than ", lenwt, " ***\n")
      cat("*** SAVE the results up to ", (i - 1), "-th RUN into IncompleteRLT.txt ***\n")
      save(rlt, file = "IncompleteRLT.txt")
      stop()
    }
    if (verbose) {
      cat("# it is done for ", i, "-th RUN\n")
    }
    
    rlt <- rbind(rlt, trlt)
  }
  
  srlt <- c()
  for (i in 1:length(wint)) {
    for (j in 1:length(wmz)) {
      td <- rlt[rlt$int == wint[i] & rlt$mz == wmz[j], ]
      if (dim(td)[1] != R) {
        stop("# The smaple size is not equal to ", R, " (", dim(td)[1], ")\n")
      }
      mskew <- mean(td$skew)
      sdskew <- sd(td$skew)
      mkurt <- mean(td$kurt)
      sdkurt <- sd(td$kurt)
      msk <- mean(td$sk)
      sdsk <- sd(td$sk)
      srlt <- rbind(srlt, c(wint[i], wmz[j], mskew, mkurt, msk, sdskew, sdkurt, sdsk))
    }
  }
  srlt <- as.data.frame(srlt)
  dimnames(srlt)[[2]] <- c("int", "mz", "skew", "kurt", "sk", "sdskew", "sdkurt", "sdsk")
  
  list(stat = srlt, mz = wmz, int = wint, R = R, seed = nseed)
}
skewness <- function(x, na.rm = FALSE) {
  # calculate skewness and kurtosis for a reference library
  # skewness calculation function from an R package "moments"
  if (is.matrix(x)) {
    apply(x, 2, skewness, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    (sum((x - mean(x))^3) / n) / (sum((x - mean(x))^2) / n)^(3 / 2)
  }
  else if (is.data.frame(x)) {
    sapply(x, skewness, na.rm = na.rm)
  } else {
    skewness(as.vector(x), na.rm = na.rm)
  }
}
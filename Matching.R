#' @export
F_CalPMMT<-function(mz,MRP=17500,RefMZ=200){
  A <- 1 / (MRP * (RefMZ^0.5))
  B <- A / 2.35482
  return(B * mz^1.5 + mz / 1000000)
  # MRP=70000
  # mz<-369.35383
}
#' @export
F_CSV<-function(Data){write.csv(Data,paste('Test.csv',seq=''),row.names = F)}
#' @export
F_FixedMatching <- function(top = top, bottom = bottom, mztol = 0.005) {## Fixed matching
  for (i in 1:nrow(bottom)) {
    top[, 1][abs(bottom[, 1][i] - top[, 1]) < mztol] <- bottom[, 1][i]
  }
  alignment <- merge(top, bottom, by = 1, all = TRUE) # use the bottom (library) as reference
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0 # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  # print(alignment)
  return(alignment)
}


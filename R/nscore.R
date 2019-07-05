#' Normal Score Transformation
#'
#' Takes a vector of values x and calculates their normal scores. Returns
# a list with the scores and an ordered table of original values and
# scores, which is useful as a back-transform table. See backtr().
#'
#' @param x numeric, vector
#'
#' @return list
#' nscore contains normal score transformed vector
#' trn.table contains an ordered table of the original values for
#' use with backtransformation
#' @export
nscore <- function(x) {
  nscore <- qqnorm(x, plot.it = FALSE)$x # normal score
  trn.table <- data.frame(x = sort(x), nscore = sort(nscore))

  return(list(nscore = nscore, trn.table = trn.table))
}


#' Back transformation function
#'
#' Given a vector of normal scores and a normal score object
# (from nscore), the function returns a vector of back-transformed
# values.
#' @param scores a vector of normal scores and a normal score object (from nscore)
#' @param nscore the nscore object
#' @param tails One major issue is how to extrapolate to the tails. Options
# other than none may result in dramatically incorrect tail estimates!
# tails options:
# 'none' : No extrapolation; more extreme score values will revert
# to the original min and max values.
# 'equal' : Calculate magnitude in std deviations of the scores about
# initial data mean. Extrapolation is linear to these deviations.
# will be based upon deviations from the mean of the original
# hard data - possibly quite dangerous!
# 'separate' :  This calculates a separate sd for values
# above and below the mean.
#' @param draw show qq plot
#'
#' @return numeric
#' @export
#' @importFrom graphics plot
#' @importFrom stats approxfun cor lm predict qqnorm sd
backtr <- function(scores, nscore, tails = "none", draw = TRUE) {
  if (tails == "separate") {
    mean.x <- mean(nscore$trn.table$x)
    small.x <- nscore$trn.table$x < mean.x
    large.x <- nscore$trn.table$x > mean.x
    small.sd <- sqrt(sum((nscore$trn.table$x[small.x] - mean.x)^2) /
      (length(nscore$trn.table$x[small.x]) - 1))
    large.sd <- sqrt(sum((nscore$trn.table$x[large.x] - mean.x)^2) /
      (length(nscore$trn.table$x[large.x]) - 1))
    min.x <- mean(nscore$trn.table$x) + (min(scores) * small.sd)
    max.x <- mean(nscore$trn.table$x) + (max(scores) * large.sd)
    
    # check to see if these values are LESS extreme than the
    # initial data - if so, use the initial data.
    # print(paste('lg.sd is:',large.sd,'max.x is:',max.x,'max nsc.x is:',max(nscore$trn.table$x)))
    if (min.x > min(nscore$trn.table$x)) {
      min.x <- min(nscore$trn.table$x)
    }
    if (max.x < max(nscore$trn.table$x)) {
      max.x <- max(nscore$trn.table$x)
    }
  }
  
  if (tails == "equal") { # assumes symmetric distribution around the mean
    mean.x <- mean(nscore$trn.table$x)
    sd.x <- sd(nscore$trn.table$x)
    min.x <- mean(nscore$trn.table$x) + (min(scores) * sd.x)
    max.x <- mean(nscore$trn.table$x) + (max(scores) * sd.x)
    
    # check to see if these values are LESS extreme than the
    # initial data - if so, use the initial data.
    if (min.x > min(nscore$trn.table$x)) {
      min.x <- min(nscore$trn.table$x)
    }
    if (max.x < max(nscore$trn.table$x)) {
      max.x <- max(nscore$trn.table$x)
    }
  }
  
  if (tails == "none") { # No extrapolation
    min.x <- min(nscore$trn.table$x)
    max.x <- max(nscore$trn.table$x)
  }
  
  min.sc <- min(scores)
  max.sc <- max(scores)
  x <- c(min.x, nscore$trn.table$x, max.x)
  nsc <- c(min.sc, nscore$trn.table$nscore, max.sc)

  if (draw) {
    plot(nsc, x, main = "Transform Function")
  }
  back.xf <- approxfun(nsc, x) # Develop the back transform function
  val <- back.xf(scores)

  return(val)
}


#' Builds stats
#' 
#' Given a hard (primary) dataset and a soft (secondary) vector
# of colocated values, return a list object with max, min,
# correlation, and linear model coefficients for prediction.
#'
#' @param hard numeric vector
#' @param soft numeric vector
#'
#' @return model object
#' @export
modelSoftHard <- function(hard, soft) {
  soft.model <- lm(hard ~ soft)
  min.max <- predict(soft.model, data.frame(soft = c(min(soft), max(soft))))
  model.object <- list(
    min = min.max[1], max = min.max[2],
    coefficients = soft.model$coeff,
    r = cor(hard, soft), summary = summary(soft.model)
  )
  return(model.object)
}

#' Writes a data.frame to a GEO_EAS text file
#' 
#' This function writes a data frame to a GEO-EAS text file.
#' Useful for interacting with GSLIB.
#' @param dat data.frame
#'
#' @param outfile character, file path
#'
#' @export
#' @importFrom utils write.table
write.data.frame.geoeas <- function(dat, outfile) {

  # Write the header
  cat("GSLIB file created in R
", file = outfile)
  cat(length(names(dat)), file = outfile, append = TRUE)
  cat("
", file = outfile, append = TRUE)
  write(cbind(names(dat)), file = outfile, append = TRUE)
  write.table(dat, file = outfile, append = TRUE, sep = "	", col.names = FALSE, row.names = FALSE)
}

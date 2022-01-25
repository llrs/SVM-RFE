#' Wrapper to run svmRFE function while omitting a given test fold
#'
#' @param test.fold Indexes of the samples.
#' @param X Input data to be subseted by test.fold.
#' @inheritParams svmRFE
#' @param \dots other arguments passed to svmRFE.
#' @return A list with contains the feature rankings for that fold (`feature.ids`),
#' as well as the training set row names used to obtain them (`train.data.ids`).
#' The remaining test set row names are included as well (`test.data.ids`)
#' @seealso [svmRFE()]
#' @export
svmRFE.wrap <- function(test.fold, X, response, ...) {
  train.data <- X[-test.fold, ]
  test.data <- X[test.fold, ]

  # Rank the features
  features.ranked <- svmRFE(train.data, ...)

  return(list(feature.ids = features.ranked,
              train.data.ids = row.names(train.data),
              test.data.ids = row.names(test.data)))
}

#' Feature ranking via Feature Elimination (RFE) algorithm
#'
#' @param X input data as 2 dimensional data, samples in rows.
#' @param response Character vector with the response vector of X or a numeric
#' value with the column id where the response variable is.
#' @param k K-fold number
#' @param halve.above Max number of remaining features
#' @return A list of the features selected.
#' @export
#' @importFrom utils txtProgressBar flush.console setTxtProgressBar
#' @importFrom stats sd na.omit
svmRFE <- function(X, response, k = 1, halve.above = 5000) {
    if (is.numeric(response) && length(response) == 1 && missing(response)) {
        stop("Response is not a numeric vector of length 1.")
    } else if (length(response) != nrow(X)) {
        stop("Response is not a character vector of length nrow(X)")
    }


  n <- ncol(X) - 1

  # Scale data up front so it doesn't have to be redone each pass
  message("Scaling data...")

  if (is.character(response) || is.factor(response)) {
      X <- scale(X)
  } else {
      response <- X[, 1]
      X[, -1] <- scale(X[, -1])
  }
  message("Done!\n")
  flush.console()

  pb <- txtProgressBar(1, n, 1, style = 3)

  i.surviving <- seq_len(n)
  i.ranked <- n
  ranked.list <- vector(length = n)

  # Recurse through all the features
  while (length(i.surviving) > 0) {
    if (k > 1) {
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)
      folds <- rep(seq_len(k), length.out = nrow(X))[sample(nrow(X))]
      folds <- split(seq_along(folds), folds)

      # Obtain weights for each training set
      w <- lapply(folds, getWeights, X[, i.surviving, drop = FALSE], response = response)
      w <- do.call(rbind, w)

      # Normalize each weights vector
      w <- t(apply(w, 1, function(x) x / sqrt(sum(x^2))))

      # Compute ranking criteria
      v <- w * w
      vbar <- apply(v, 2, mean)
      vsd <- apply(v, 2, sd)
      c <- vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w <- getWeights(NULL, X[, i.surviving, drop = FALSE], response = response)
      c <- w * w
    }

    # Rank the features
    ranking <- order(c)
    if (length(i.surviving) == 1) {
      ranking <- 1
    }

    if (length(i.surviving) > halve.above) {
      # Cut features in half until less than halve.above
      nfeat <- length(i.surviving)
      ncut <- round(nfeat / 2)
      n <- nfeat - ncut

      message("\nFeatures halved from ", nfeat, " to ", n, "\n")
      flush.console()

      pb <- txtProgressBar(1, n, 1, style = 3)
    } else {
      ncut <- 1
    }

    # Update feature list
    ranked.list[i.ranked:(i.ranked - ncut + 1)] <- i.surviving[ranking[1:ncut]]
    i.ranked <- i.ranked - ncut
    i.surviving <- i.surviving[-ranking[1:ncut]]

    setTxtProgressBar(pb, n - length(i.surviving))
    flush.console()
  }

  close(pb)

  return(ranked.list)
}

#' Fit a linear SVM model and obtain feature weights
#' @param test.fold Index of the input data used for testing
#' @param X Train data.
#' @return A matrix with the weights of the SVM
#' @seealso [svm()]
#' @export
#' @importFrom e1071 svm
getWeights <- function(test.fold, X, response) {
  train.data <- X
  if (missing(response) || is.null(response)) {
      response <- 1
  }
  if (!is.null(test.fold)) {
      if (is.null(nrow(X)) || nrow(X) == 0) {
          browser()
      }
      train.data <- X[-test.fold, , drop = FALSE]
      response <- response[-test.fold]
  }
  if (!(is.character(response) || is.factor(response))) {
      train.data <- train.data[, -response, drop = FALSE]
  }

  svmModel <- svm(train.data, response,
    cost = 10, cachesize = 500,
    scale = FALSE, type = "C-classification", kernel = "linear"
  )

  t(svmModel$coefs) %*% svmModel$SV
}

#' Compile feature rankings across multiple folds
#' @param results A list of `svmRFE.wrap` results.
#' @param input The input data
#' @param save Logical value whether to create a file or not.
#' @param file The file name of the file to be created.
#' @return Either a file (if save is true) or a matrix with the information about the features.
#' @export
#' @importFrom utils write.table
WriteFeatures <- function(results, input, save = TRUE, file = "features_ranked.txt") {
  # Compile feature rankings across multiple folds
  featureID <- order(rowMeans(sapply(results, function(x) order(x$feature))))
  avg.rank <- order(rowMeans(sapply(results, function(x) order(x$feature))))
  feature.name <- colnames(input[, -1])[featureID]
  features.ranked <- data.frame(FeatureName = feature.name, FeatureID = featureID, AvgRank = avg.rank)
  if (save == T) {
    write.table(features.ranked, file = file, quote = F, row.names = F)
  } else {
    features.ranked
  }
}

#' Estimate generalization error
#'
#' Does so across all hold-out folds, for a given number of top features.
#' @param i Top  i features to estimate the error.
#' @param results A list of `svmRFE.wrap` results.
#' @param input The original data
#' @return Contains the generalization error estimates in the external 10-fold CV.
#' These accuracies are averaged as `error`.
#' @export
#' @importFrom e1071 tune tune.control
FeatSweep.wrap <- function(i, results, input) {
  # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
  svm.list <- lapply(results, function(x) {
    tune(svm,
      train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
      train.y = input[x$train.data.ids, 1],
      validation.x = input[x$test.data.ids, 1 + x$feature.ids[1:i]],
      validation.y = input[x$test.data.ids, 1],
      # Optimize SVM hyperparamters
      ranges = tune(svm,
        train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
        train.y = input[x$train.data.ids, 1],
        ranges  = list(gamma = 2^(-12:0), cost = 2^(-6:6))
      )$best.parameters,
      tunecontrol = tune.control(sampling = "fix")
    )$performances
  })

  error <- mean(sapply(svm.list, function(x) x$error))
  return(list(svm.list = svm.list, error = error))
}

#' Plot average generalization error vs. number of top features
#' @param errors The results of `FeatSweep.wrap`.
#' @param errors2 A different result of `FeatSweep.wrap` to be compared with (colored in gray).
#' @param no.info The `min(prop.table(table(input[,1])))`
#' @param ylim Graphical parameter to limit the size of the plot
#' @param xlab Name of the x axis
#' @param ylab Name of the y axis
#' @return Called for its side effects of plotting.
#' @export
#' @importFrom graphics points text plot abline lines
PlotErrors <- function(errors, errors2 = NULL, no.info = 0.5,
                       ylim = range(c(errors, errors2), na.rm = T),
                       xlab = "Number of Features", ylab = "10x CV Error") {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col = "black") {
    lines(which(!is.na(errors)), na.omit(x), col = col)
    points(which.min(x), min(x, na.rm = T), col = "red")
    text(which.min(x), min(x, na.rm = T), paste(which.min(x), "-", format(min(x, na.rm = T), dig = 3)), pos = 4, col = "red", cex = 0.75)
  }

  plot(errors, type = "n", ylim = ylim, xlab = xlab, ylab = ylab)
  AddLine(errors)
  if (!is.null(errors2)) AddLine(errors2, "gray30")
  abline(h = no.info, lty = 3)
}

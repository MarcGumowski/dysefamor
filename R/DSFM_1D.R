# ------------------------------------------------------------------------------ #
# --------------- Yield Curve Estimation using DSFM1D -------------------------- #
# ------------------------------------------------------------------------------ #

#######
# ------------------------------------------------------------------------------ #
# DSFM One-dimension ----------------------------------------------------------- #
# ------------------------------------------------------------------------------ #
#######

#' Estimation of Dynamic Semiparametric Factor Model for One-Dimensional Data
#'
#' \code{DSFM1D} performs a model estimation using Dynamic Semiparametric Factor
#' mechanics with one-dimensional covariates. This function is called by the
#' \code{\link{DSFM}} main function.
#'
#' Dynamic Semiparametric Factor Models (DSFM) are defined as
#' \eqn{ Y_{t,j} = m_0(X_{t,j}) + \sum_{l=1}^L Z_{t,l} m_l(X_{t,j}) +
#' \varepsilon_{t,j}}.
#' DSFM estimation is performed using kernel density for the non-parametric
#' functions \eqn{m_l}. The estimation is performed using  the iterative
#' algorithm of Fengler and al. (2007).
#'
#' The function has predefined arguments that can be changed for better
#' approximation.
#' First, the number of data points on the estimation grid is set
#' to 25. Larger grid can significantly increase the computation time without
#' necesseraly improve the fit.
#' Secondly, the bandwidth \eqn{h} is basically set to 0.05 but optimal
#' bandwidth has to be found externally. The algorithm always normalize the
#' covariates to work on an estimation grid bounded beetween [0,1].
#'
#' For model selection, different criteria are computed.
#'
#' For number of factors selection, the function compute the Explained Variance,
#' for bandwidth selection, two criteria are computed, a weighted Akaike
#' Information Criterion (AIC) and a weighted Schwarz Criterion (SC).
#' The goodness-of-fit is measured by the Root Mean Squared Error (RMSE). The
#' proper definition of each criterion can be found in references.
#'
#' @param data a matrix containing time indicator in first row, value
#' \eqn{Y_{t,j}} in second row, and the coordinates \eqn{X_{t,j}} in the
#' remaining row. Proper formatting has to be done using the
#' \code{\link{dataDSFM1D}}.
#' @param numDataPoints the number of points in the axis of the grid to perform
#' the kernel density estimation.
#' @param h the bandwidth used to perform the kernel density estimation. Can be
#' either a single global parameter, or a vector of the same length of
#' numDataPoints to perform local kernel estimation.
#' @param L the number of underlying factors.
#' @param initialLoad the type of initial loadings to choose between White Noise
#' \code{"WN"}, and AR(1) process \code{"AR"}. Required as starting value of the
#' algorithm. Changing the \code{initialLoad} can sligthly improve the
#' convergence rate.
#' @param tol the tolerance for the algorithm to stop.
#' @param maxIt the maximum number of iterations for the algorithm to break.
#'
#' @return \code{DSFM1D} returns an object of class \code{"DSFM1D"}.
#'
#' The generic functions \code{print}, \code{summary}, \code{plot} and
#' \code{predict} are available.
#'
#' An object of class \code{"DSFM1D"} is a list containing the following
#' components:
#' \item{\code{Data}}{the input data.}
#' \item{\code{Y}}{the input data in a more usual format, i.e. a matrix with a
#' time indicator as first row and the following rows being the value
#' \eqn{Y_{t,j}} for each covariates \eqn{X_{t,j}}.}
#' \item{\code{YHat}}{the estimated \eqn{\hat{Y}_{t,j}} with the same format,
#' i.e. a matrix with a time indicator as first row and the following rows being
#' the value \eqn{\hat{Y}_{t,j}} for each covariates \eqn{X_{t,j}}.}
#' \item{\code{ZHat}}{the estimated factor loadings \eqn{\hat{Z}_{t,j}}.}
#' \item{\code{mHat}}{the estimated factor functions \eqn{\hat{m}_l}.}
#' \item{\code{EV}}{gives the Explained Variance, used to select the approriate
#' number of factors.}
#' \item{\code{RMSE}}{gives the Root Mean Squared Error,
#' used to compare models.}
#' \item{\code{Bw}}{gives the bandwidth \eqn{h} used and two selection criteria
#'  to select the optimal bandwidth.}
#' \item{\code{x1}}{the vector of the covariates.}
#' \item{\code{Density}}{the kernel density estimation performed.}
#' \item{\code{Convergence}}{the value of the algorithm stopping criterion at
#' each loop.}
#' \item{\code{Time}}{an indicator of the time taken by the function to perform
#' the fit.}
#'
#' @author The implementation of model by Marc Gumowski was based on
#' Fengler and al. (2007).
#'
#' @references Borak, Szymon, Matthias R. Fengler, and Wolfgang K. Haerdle
#' (2005). "DSFM Fitting of Implied Volatility Surfaces". In: \emph{5th
#' International Conference on Intelligent Systems Design and Applications
#' (ISDA'05)}, pp. 526-531. IEEE.
#'
#' Fengler, Matthias R, Wolfgang K Haerdle, and Enno Mammen (2007).
#' "A Semiparametric Factor Model for Implied Volatility Surface Dynamics".
#' In: \emph{Journal of Financial Econometrics 5.2}, pp. 189-218.
#'
#' Haerdle, Wolfgang K., and Piotr Majer (2014).
#' "Yield Curve Modeling and Forecasting using Semiparametric Factor Dynamics".
#' In: \emph{The European Journal of Finance}, pp. 1-21.
#'
#' @seealso \code{\link{summary.DSFM1D}} for
#' summaries and \code{\link{plot.DSFM1D}} for plot
#' possibilities.
#'
#' \code{\link{predict.DSFM1D}} provide succint
#' predictions.
#'
#' \code{\link{dataDSFM1D}} has to be used before
#' using the \code{\link{DSFM}} function to ensure that the data are correctly
#' formated.
#'
#' \code{\link{simulateDSFM1D}} is a function to simulate data that can be used
#' as simple example purposes.
#'
#' @examples
#' # Prepare the data --------------------------------------------------------- #
#' # Interest rate of zero-coupon bond yield curves. Data from Bank of Canada.
#' data(canadianYieldCurves)
#' maturity <- c(1/4, 1/2, 3/4, 1:10, 20, 30)
#' dsfmData <- dataDSFM1D(canadianYieldCurves[1:100, ], maturity)
#' dsfmData
#' plot(dsfmData)
#'
#' # Set the parameters ------------------------------------------------------- #
#' h        <- 0.167
#' L        <- 3
#'
#' # Fit the model ------------------------------------------------------------ #
#' dsfmFit  <- DSFM(dsfmData, h = h, L = L)
#' summary(dsfmFit)
#' plot(dsfmFit)
#'
#' # Perform prediction ------------------------------------------------------- #
#' horizon  <- 5
#' predict(dsfmFit, nAhead = horizon)
#'
#' @export
#'
DSFM1D <- function(data, numDataPoints = 25, h = 0.5, L = 3,
                   initialLoad = "WN", tol = 1e-5, maxIt = 301) {

  Time1 <- Sys.time()                                 # Get starting time

  # Initial Settings
  date    <- unique(data$Date)                        # Get the dates
  I       <- length(date)                             # T
  L0      <- L + 1                                    # L with factor L_0
  y       <- matrix(data$y, I, byrow = T)          # Get the y
  namesX1 <- unique(data$x1)
  # Normalization of the covariates
  x1      <- unique(data$x1) / max(unique(data$x1))
  J       <- length(x1)

  # Initial Loadings Z_t,j
  switch(initialLoad,
         AR = {
           # AR(1), beta = 0.4
           ZHat <- matrix(0, I, L)
           for (l in 1:L) {
             ZHat[ ,l] <- arima.sim(list(ar = 0.4, ma = 0), I)
           }
         },
         WN = {
           # White Noise N(0,1)
           ZHat <- matrix(rnorm(I * L, 0, 1), I, L)
         },
         stop("Invalid initialLoad Parameter")
  )
  ZHat   <- cbind(rep(1, I), ZHat)

  # Create a regular grid of points u covering the whole space
  minX   <- min(x1)
  maxX   <- max(x1)
  delta  <- (maxX - minX) / (numDataPoints - 1)
  u      <- seq(minX, maxX, delta)
  U      <- length(u)

  # Loop parameters
  ZHatOld    <- matrix(0, I, L0)
  mHatOld    <- matrix(0, U, L0)
  it         <- 0
  stopCriterion     <- 1
  plotStopCriterion <- c()

  # KERNEL -------------------------------------------------------------------- #

  if (length(h) == 1) {
    h <- rep(h, numDataPoints)
  }
  Kernel <- KernelDensity1D(y, I, J, x1, u, U, h)

  # LOOP ---------------------------------------------------------------------- #

  while (stopCriterion >= tol) {

    it = it + 1

    if (it == maxIt) {
      stop("Too many iterations, no convergence")
    }

    # Matrix B
    Bu <- matrix(0, L0, L0)
    BList <- replicate(U, list(Bu))
    for (n in 1:U){
      BList[[n]] <- J * t(ZHat) %*% (ZHat * Kernel$pHat[ ,n])
    }

    # Vector Q
    Q <- vector(mode="numeric", L0)
    QList <- replicate(U, list(Q))
    for (n in 1:U) {
      QList[[n]] <- J * t(ZHat) %*% Kernel$qHat[ ,n]
    }

    # Vector mHat for each grid points of length L
    mHat <- mapply(function(x, y) x %*% y, lapply(BList, solve), QList)
    mHat <- t(mHat)

    # Matrix M
    M <- matrix(0, L, L)
    MList <- replicate(I, list(M))
    for (t in 1:I) {
      MList[[t]] <- (t(mHat[ ,2:L0]) %*% (mHat[ ,2:L0] * Kernel$pHat[t, ]))  *
        delta
    }

    # Vector S
    S  <- vector(mode="numeric", L)
    SList <- replicate(I, list(S))
    for (t in 1:I) {
      S1 <- colSums(Kernel$qHat[t, ] * data.frame(mHat[ ,2:L0]))
      S2 <- colSums(c(Kernel$pHat[t, ] * data.frame(mHat[ ,1])) *
                      data.frame(mHat[ ,2:L0]))
      SList[[t]] <- (S1 - S2) * delta
    }

    # Vector ZHat for each date t of length L
    ZHat <- mapply(function(x, y) x %*% y, lapply(MList, solve), SList)
    if (L0 > 2) {
      ZHat <- t(ZHat)
    }
    ZHat <- cbind(1, ZHat)


    # Stopping Criterion
    stopCriterionAll <- vector(mode="numeric",I)
    for (i in 1:I) {
      stopCriterionAll[[t]] <- sum((rowSums(ZHat[t, ] * mHat - ZHatOld[t, ] *
                                              mHatOld)^2))  * delta
    }
    stopCriterion <- sum(stopCriterionAll)

    ZHatOld <- ZHat
    mHatOld <- mHat
    plotStopCriterion[it] <- stopCriterion
  }

  # Orthogonalization --------------------------------------------------------- #

  # First Step ------- #

  # pHat
  pHat <- (1 / I) * colSums(Kernel$pHat)

  # Vector gamma
  g <- colSums(mHat[ ,1] * as.matrix(mHat[ ,2:L0]) * pHat) * delta

  # Matrix Gamma
  G <- (t(mHat[ ,2:L0]) %*% (mHat[ ,2:L0] * pHat)) * delta

  # mHat new
  mHat[ ,1] <- mHat[ ,1] - t((t(g) %*% solve(G)) %*% t(mHat[ ,2:L0]))
  mHat[ ,2:L0] <- t((G%^%(-0.5)) %*% t(mHat[ ,2:L0]))

  # ZHat new
  for (t in 1:I){
    ZHat[t,2:L0] <- t(G%^%(0.5) %*% (ZHat[t,2:L0] + solve(G) %*% g))
  }

  # Second Step ------- #

  # Matrix B
  B <- t(ZHat[ ,2:L0]) %*% ZHat[ ,2:L0]

  # Eigenvectors Z
  Z <- eigen(B)$vectors

  # mHat new
  mHat[ ,2:L0] <- t(t(Z) %*% t(mHat[ ,2:L0]))

  # ZHat new
  ZHat[ ,2:L0] <- t(t(Z) %*% t(ZHat[ ,2:L0]))

  # Global fit of the model --------------------------------------------------- #

  YBar <- mean(t(t(y)))
  # Interpolation of function m_l
  m_l <- rep(0, J)
  m_l <- replicate(L0, list(m_l))
  for (l in 1:L0) {
    m_l[[l]] <- approx(u, mHat[ ,l], xout = x1, rule = 2)$y
  }
  # Explained Variance - Model Size Selection ----------- #
  YHat <- matrix(0, I, J)
  for (t in 1:I) {
    YTHat <- list()
    for (l in 1:L0) {
      YTHat[[l]] <- ZHat[t,l] * m_l[[l]]
    }
    YTHat <- Reduce("+",YTHat)
    YHat[t, ] <- YTHat
  }
  numerator <- (y - YHat)^2
  numerator <- sum(numerator)
  denominator <- (y - YBar)^2
  denominator <- sum(denominator)
  EV <- 1 - numerator / denominator

  # Goodness-of-fit - Root Mean Squared Error - RMSE --- #
  RMSE <- data.frame(sqrt((1 / (I * J)) * numerator))
  names(RMSE)     <- "RMSE"

  # Fit of the sliced model --------------------------------------------------- #

  if (length(unique(h)) == 2) {

    Slice       <- length(h[which(h == unique(h)[1])])
    xS          <- list()
    xS[[1]]     <- which(x1 <= u[Slice])
    xS[[2]]     <- which(x1 >  u[Slice])

    RMSESlice <- c()
    for (i in 1:2) {
      numerator <- (y[ ,xS[[i]]] - YHat[ ,xS[[i]]])^2
      numerator <- sum(numerator)

      # Goodness-of-fit - Root Mean Squared Error - RMSE[[i]] -- #
      RMSESlice[i] <- sqrt((1 / (I* length(xS[[i]]))) * numerator)
    }
    RMSE  <- data.frame(t(c(RMSE, RMSESlice)))
    names(RMSE) <- c("RMSE - Overall",
                     paste("RMSE - Slice", 1:(length(RMSE) - 1)))
  }


  # Bandwidth Selection  ------------------------------- #
  w     <- 1 / pHat
  N     <- I * J
  if (length(unique(h)) == 1) {
    Kh  <- (1 / unique(h)) * QuarticKernel1D(0 / unique(h))
    # Weighted AIC_2
    AIC2 <- (1 / N) * numerator * exp(2 *
                                        (L / N) * Kh *
                                        ((sum(w) * delta) /
                                           (sum(w * pHat)) * delta))
    # Weighted SC_1
    SC1  <- (1 / N) * numerator * exp(log(N) *
                                        (L / N) * Kh *
                                        ((sum(w) * delta) /
                                           (sum(w * pHat)) * delta))
  } else if (length(unique(h)) == 2) {
    # Special Case - 2 Different Bandwidths
    Kh    <- (length(h[which(h == unique(h)[1])]) / length(h)) *
      (1 / unique(h)[1]) * QuarticKernel1D(0 / unique(h)[1]) +
      (length(h[which(h == unique(h)[2])]) / length(h)) *
      (1 / unique(h)[2]) * QuarticKernel1D(0 / unique(h)[2])
    # Weighted AIC_2
    AIC2 <- (1 / N) * numerator * exp(2 *
                                        (L / N) * Kh *
                                        ((sum(w) * delta) /
                                           (sum(w * pHat) * delta)))
    # Weighted SC_1
    SC1  <- (1 / N) * numerator * exp(log(N) *
                                        (L / N) * Kh *
                                        ((sum(w) / (sum(w * pHat)) * delta)))
  } else {
    # If Local Bandwidth does not compute selection criteria then
    AIC2 <- NA
    SC1  <- NA
  }



  # Outputs ------------------------------------------------------------------- #

  # Create data frames as outputs
  ZHat            <- data.frame(date, ts(ZHat))
  names(ZHat)     <- c("Date", paste0("Z_t", 0:L, ".hat"))
  mHat            <- data.frame(u * max(unique(data$x1)), mHat)
  names(mHat)     <- c("u1", paste0("m_", 0:L, ".hat"))
  y               <- as.matrix(y, I, J)
  y               <- data.frame(date, y)
  names(y)        <- c("Date", namesX1)
  YHat            <- data.frame(date, YHat)
  names(YHat)     <- names(y)
  EV              <- data.frame(EV)
  rownames(EV)    <- NULL
  names(EV)       <- paste0("EV(L = ", L, ")")
  Bw              <- data.frame(t(c(unique(h), AIC2, SC1)))
  names(Bw)       <- c(paste0("h", 1:length(unique(h))), "wAIC_2", "wSC_1")

  Time2 <- Sys.time() - Time1

  model <- list(Data = data, Y = y, YHat = YHat, ZHat = ZHat, mHat = mHat,
                EV = EV, RMSE = RMSE, Bw = Bw, x1 = namesX1, Density = pHat,
                Convergence = plotStopCriterion, Iterations = it, Time = Time2)
  model$call   <- match.call()
  class(model) <- "DSFM1D"

  return(model)
}





#######
# ----------------------------------------------------------------------------- #
# S3 Methods (print.DSFM1D, summary.DSFM1D, plot.DSFM1D, predict.DSFM1D) ------ #
# ----------------------------------------------------------------------------- #
#######

# Print Function -------------------------------------------------------------- #

#' @export
#'
print.DSFM1D <- function(x, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\nY:\n")
  print(x$Y)
  cat("\nEstimated Y:\n")
  print(x$YHat)
  cat("\nEstimated Z:\n")
  print(x$ZHat)
  cat("\nEstimated m:\n")
  print(x$mHat)
  cat("\nExplained Variance:\n")
  print(x$EV, row.names = F)
  cat("\nRoot Mean Squared Error:\n")
  print(x$RMSE, row.names = F)
  cat("\nBandwidth:\n")
  print(x$Bw, row.names = F)
  cat("\nCovavriates:",x$x1,"\n")
  cat("\nIterations:",x$Iterations,"\n")
  cat("\n")
  print(x$Time)
}

# Summary Function ------------------------------------------------------------ #

#' Summarizing DSFM Fits for One-Dimensional Data
#'
#' \code{summary} method for class \code{"DSFM1D"}.
#'
#' @param object an object of class \code{"DSFM1D"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function \code{summary.DSFM1D} returns a list of summary
#' statistics of the fitted DSFM given in \code{object}, using the components
#' (list elements) \code{"call"} from its arguments, plus:
#' \item{\code{ZHat}}{a summary of the estimated factor loadings.}
#' \item{\code{EV}}{the Explained Variance, used to select the approriate
#' number of factors.}
#' \item{\code{RMSE}}{the Root Mean Squared Error, used to compare the goodness-
#' of-fit between models.}
#' \item{\code{Bw}}{The bandwidth and its selection criteria.}
#'
#' @seealso \code{\link{dataDSFM1D}}, \code{\link{DSFM}}, \code{\link{DSFM1D}},
#' \code{\link{plot.DSFM1D}}, \code{\link{predict.DSFM1D}}.
#'
#' @export
#'
summary.DSFM1D <- function(object, ...) {

  ZHat <- object$ZHat[ ,2:length(object$ZHat)]
  L0   <- length(ZHat)
  # Table of the statistics of the Z_t,j
  loadingsSummary <- data.frame(rbind(apply(ZHat[2:L0], 2, summary),
                                      apply(ZHat[2:L0], 2, sd),
                                      apply(ZHat[2:L0], 2, skewness),
                                      apply(ZHat[2:L0], 2, kurtosis)),
                                check.names = T, row.names = NULL)
  row.names(loadingsSummary) <- c("Min.", "1st Qu.", "Median", "Mean",
                                  "3rd Qu.", "Max.", "Std.", "Skew.", "Kurt.")

  cat("Call:\n")
  print(object$call)

  cat("\nNumber of observations: ")
  cat(dim(object$Y)[1], "\n")

  cat("\nEstimated m:\n")
  print(object$mHat)

  cat("\nSummary of the estimated factor loadings:\n")
  print(loadingsSummary)

  cat("\nExplained Variance:\n")
  print(object$EV, row.names = F)

  cat("\nRoot Mean Squared Error:\n")
  print(object$RMSE, row.names = F)

  if (is.na(object$Bw$wAIC_2) == F) {
    cat("\nBandwidth selection criteria:\n")
    print(object$Bw, row.names = F)
  }
}

# Plot Function  -------------------------------------------------------------- #

#' Plot Diagnostics for a DSFM1D Object
#'
#' Four plots (selectable by \code{which}) are currently available: a plot of
#' the convergence rate of the algorithm, a plot of the time series of the
#' estimated loadings for each factor, a plot of the estimated factor functions,
#' and a plot of the fit. By default, the four plots are provided.
#'
#' @param x an \code{object} of class \code{"DSFM1D"}.
#' @param which to choose between \code{"all"},\code{"loadings"}.
#' \code{"functions"},\code{"convergence"}, and \code{"fit"}.
#' @param ask logical; if TRUE, the user is asked before each plot, see
#' \code{\link{par}}(ask=.).
#' @param pal the color palette for the fit plot. To choose between
#' \code{"pink"},\code{"blue"},\code{"light"},\code{"dark"}.
#' @param col what color should be drawn.
#' @param type what type of plot should be drawn:
#' see \code{\link{plot.default}}.
#' @param theta angles defining the viewing direction of the fit plot.
#' @param border the color of the line drawn around the surface facets of the
#' fit plot. The default, NULL, corresponds to par("fg"). A value of NA will
#' disable the drawing of borders: this is sometimes useful when the surface is
#' shaded.
#' @param box should the bounding box for the surface be displayed.
#' The default is \code{FALSE}.
#' @param shade the shade at a surface facet is computed as ((1+d)/2)^shade,
#' where d is the dot product of a unit vector normal to the facet and a unit
#' vector in the direction of a light source. Values of shade close to one
#' yield shading similar to a point light source model and values close to zero
#' produce no shading. Values in the range 0.5 to 0.75 provide an approximation
#' to daylight illumination.
#' @param expand a expansion factor applied to the z coordinates. Often used
#' with 0 < expand < 1 to shrink the plotting box in the z direction.
#' @param ticktype character: "simple" draws just an arrow parallel to the axis
#' to indicate direction of increase; "detailed" draws normal ticks as per 2D
#' plots.
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @seealso \code{\link{dataDSFM1D}}, \code{\link{DSFM}}, \code{\link{DSFM1D}},
#' \code{\link{summary.DSFM1D}}, \code{\link{predict.DSFM1D}}.
#'
#' @export
#'
plot.DSFM1D <- function(x, which = "all", ask = TRUE, pal = "pink",
                        col = "#016C59", type = "l", theta = 40, border = NA,
                        box = T, shade = .2, expand = .5, ticktype = "detailed",
                        ...) {

  date  <- unique(x$Y[ ,1])
  I     <- dim(x$Y)[1] # T
  L     <- dim(x$ZHat)[2] - 2
  L0    <- L + 1
  numDataPoints <- dim(x$mHat)[1]
  u     <- x$mHat[,1]

  par(ask = ask)

  if ("all" %in% which) {
    which = c("convergence", "loadings", "functions", "fit")
  }
  # Convergence Plot
  if ("convergence" %in% which) {
    N <- length(x$Convergence)
    xAxis <- 1:N
    layout(matrix(1))
    plot(xAxis[(N - N*0.95):N], x$Convergence[(N - N * 0.95):N],
         main = "Convergence of the algorithm",
         ylab = "Stopping Criterion", xlab = "Iteration",
         col = col, type = type, ...)
    layout(matrix(1))
  }
  # Factor Loadings
  if ("loadings" %in% which) {
    labDates <- seq(as.Date(date[1]), as.Date(tail(date, 1)), by = "year")
    layout(matrix(1:L, L, 1))
    for (l in 1:L) {
      plot(as.Date(date), x$ZHat[ ,l + 2],
           main = bquote(hat(Z)[t*","*.(l)]),
           ylab = "", xlab = "Time", xaxt ="n",
           col = col, type = type, ...)
      axis.Date(1, date, at = labDates)
    }
    layout(matrix(1))
  }
  # Factor Functions
  if ("functions" %in% which) {
    layout(matrix(1:(L0 + L0 %% 2), L0 %% 2 + L0 %/% 2, byrow = T))
    for (l in 1:L0) {
      plot(u, x$mHat[[l + 1]],
           main = bquote(hat(m)[.(l - 1)]),
           ylab = "", xlab = "u",
           col = col, type = type, ...)
    }
    layout(matrix(1))
  }
  # Y
  if ("fit" %in% which) {
    I     <- dim(x$Y)[1]
    J     <- dim(x$Y)[2] - 1
    x1    <- seq(1, I) / I
    x2    <- x$x1
    Y     <- t(t(x$Y[ ,2:dim(x$Y)[2]]))
    YHat  <- t(t(x$YHat[ ,2:dim(x$Y)[2]]))

    switch(pal,
           pink = {
             jetColors <- colorRampPalette(c("#1E0000", "#402525",
                                             "#C27E7E", "#E8E8B4", "#FFFFFF"))
           },
           blue = {
             jetColors <- colorRampPalette(c("#F6EFF7", "#BDC9E1",
                                             "#67A9CF", "#1C9099", "#016C59"))
           },
           light = {
             jetColors <- colorRampPalette(c("#3FB8AF", "#7FC7AF",
                                             "#DAD8A7", "#FF9E9D", "#FF3D7F"))
           },
           dark = {
             jetColors <- colorRampPalette(c("#556270", "#4ECDC4",
                                             "#C7F464", "#FF6B6B", "#C44D58"))
           },
           stop('Invalid pal Parameter. Valid Parameters are: "pink", "blue",
                "light, "dark".')
    )
    color <- jetColors(200)

    layout(matrix(1:2, 1, 2))

    nrz <- nrow(Y)
    ncz <- ncol(Y)
    zFacet <- (Y[-1, -1] + Y[-1, -ncz] + Y[-nrz, -1] + Y[-nrz, -ncz]) / 4
    facetCol <- cut(zFacet, 200)
    persp(x1, x2, Y, xlab = "Time", ylab = "x", zlab = "Y",
          col = color[facetCol], main ="Y", theta = theta, border = border,
          box = box, shade = shade, expand = expand, ticktype = ticktype, ...)

    nrz <- nrow(YHat)
    ncz <- ncol(YHat)
    zFacet <- (YHat[-1, -1] + YHat[-1, -ncz] + YHat[-nrz, -1] +
                 YHat[-nrz, -ncz]) / 4
    facetCol <- cut(zFacet, 200)
    persp(x1, x2, YHat, xlab = "Time", ylab = "x", zlab = "Y",
          col = color[facetCol], main = expression(hat(Y)),
          theta= theta, border = border, box = box, shade = shade,
          expand = expand, ticktype = ticktype, ...)
    layout(matrix(1))
  }

  layout(matrix(1))
  par(ask = FALSE)
}

# Predict Function  ----------------------------------------------------------- #

#' Predict Method for One-Dimensional DSFM Fits
#'
#' Predicted values based on DSFM object using Vector Autoregressive Processes
#' (VAR(p)).
#'
#' This function makes uses of package \code{vars} to fit and predict VAR(p)
#' processes, before pluging in the predicted factor loadings in the DSFM. The
#' factors functions are interpolated by \code{\link{approx}}.
#'
#' @param object an object of class \code{"DSFM1D"}.
#' @param nAhead the number of steps ahead for which prediction is required.
#' @param p the order of the Vector Autoregressive Process to be fitted.
#' @param ... other parameters to be passed through the \code{\link{VAR}}
#' and \code{\link{predict.varest}} functions.
#'
#' @return \code{predict.DSFM1D} returns an object of class
#' \code{"predict.DSFM1D"}. This class is a list containing:
#' \item{\code{YHatForecast}}{the forecasted responses.}
#' \item{\code{YHatForecastMatrix}}{the forecasted responses in a more usual
#' format.}
#' \item{\code{ZHatForecast}}{the predicted factors loadings.}
#' \item{\code{nAhead}}{the number of steps ahead.}
#'
#' @seealso \code{\link{VAR}},\code{\link{predict.varest}}, \code{\link{DSFM}},
#' \code{\link{DSFM1D}},\code{\link{dataDSFM1D}}.
#'
#' @references Bernhard Pfaff (2008). VAR, SVAR and SVEC Models: Implementation
#' Within R Package vars. In: \emph{Journal of Statistical Software 27(4)}.
#' URL http://www.jstatsoft.org/v27/i04/.
#'
#' @export
#'
predict.DSFM1D <- function(object, nAhead = 12, p = 1, ...) {

  date  <- object$Y[ ,1]
  if (date[2]-date[1] <= 1) {
    timeDiff <- "day"
  } else if (date[2] - date[1] == 7) {
    timeDiff <- "week"
  } else if (date[2] - date[1] <= 31 & date[2] - date[1] >= 28) {
    timeDiff <- "month"
  } else if (date[2] - date[1] <= 123 & date[2] - date[1] >= 118) {
    timeDiff <- "quarter"
  } else { timeDiff <- "year"
  }
  dateForecast  <- seq.Date(date[length(date)],
                            length.out = nAhead + 1,
                            by = timeDiff)[2:(nAhead + 1)]

  # Initial parameters
  ZHat  <- object$ZHat[ ,3:length(object$ZHat)]
  mHat  <- object$mHat[ ,2:length(object$mHat)]
  L     <- dim(ZHat)[2]
  L0    <- L + 1
  x1    <- object$x1
  J     <- length(x1)
  u1    <- object$mHat[ ,1]

  # VAR(p) estimation of the loadings Z_tl^forecast using package vars
  ZVarEst         <- vars::VAR(ZHat,p = p, ...)
  ZVarEstForecast <- predict(ZVarEst, n.ahead = nAhead, ...)
  ZVarEstForecast <- lapply(ZVarEstForecast$fcst, function(x) x[ ,1])
  ZVarEstForecast <- do.call(cbind, ZVarEstForecast)
  ZHatForecast    <- cbind(1, ZVarEstForecast)

  # Interpolation of function m_l
  m_l <- rep(0, J)
  m_l <- replicate(L0, list(m_l))
  for (l in 1:L0) {
    m_l[[l]] <- approx(u1, mHat[ ,l], xout = x1)$y
  }

  # YHat^forecast
  YHatForecastMatrix <- matrix(0, nAhead, J)
  for (t in 1:nAhead) {
    YTHatForecastList <- list()
    for (l in 1:L0) {
      YTHatForecastList[[l]] <- ZHatForecast[t,l] * m_l[[l]]
    }
    YTHatForecastMatrix <- Reduce("+", YTHatForecastList)
    YHatForecastMatrix[t, ] <- YTHatForecastMatrix
  }

  ZHatForecast         <- data.frame(dateForecast, ts(ZHatForecast))
  names(ZHatForecast)  <- c("Date", paste0("Z_t", 0:L, ".hat"))
  YHatForecastMatrix   <- data.frame(dateForecast, YHatForecastMatrix)
  names(YHatForecastMatrix)  <- names(object$YHat)
  YHatForecast         <- dataDSFM1D(YHatForecastMatrix,object$x1)

  predict <- list(YHatForecast = YHatForecast,
                  YHatForecastMatrix = YHatForecastMatrix,
                  ZHatForecast = ZHatForecast,
                  nAhead = nAhead)
  predict$call   <- match.call()
  class(predict) <- "predict.DSFM1D"

  return(predict)
}

# Print Function -------------------------------------------------------------- #

#' @export
#'
print.predict.DSFM1D <- function(x, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\nForecasted Y:\n")
  print(x$YHatForecastMatrix)
  cat("\nForecasted Z:\n")
  print(x$ZHatForecast)
  cat("\nHorizon:", x$nAhead)

}






#######
# ----------------------------------------------------------------------------- #
# Simulation ------------------------------------------------------------------ #
# ----------------------------------------------------------------------------- #
#####



# 1 Dimension Simulation ------------------------------------------------------ #

#' Simulate Responses for One-Dimensional Factor Models
#'
#' This function simulates responses for one-dimensional factor models.
#'
#' Two different way of generating data are available: using an Extended
#' Nelson-Siegel model following Bliss (1997) \code{"ns"} with a predefined
#' VAR(1) process, or a Dynamic Semiparametric Factor Model \code{"dsfm"}
#' with predefined factors functions.
#' This function is used for example purpose, only few parameters are available
#' to control the simulation.
#'
#' The starting values to simulate the Extended
#' Nelson-Siegel model are taken following Linton and al. (2001).
#'
#' The factors loadings of the DSFM are simulated from three independant AR(1)
#' process and the factors functions are predefined to be orthogonals.
#'
#' @param model the type of model to be generated. To choose between \code{"ns"}
#' and \code{"dsfm"}.
#' @param n the number of observations.
#' @param x1 a vector of covariates.
#' @param L the number of factors for the DSFM.
#' @param var the error \eqn{\varepsilon_{t,j}} of the models. Allows to control
#' the noise.
#' @param beta a vector of dimension (1x3) controlling the starting values of
#' the VAR(1) process.
#' @param tau a vector of dimension (1x2) specifying the parameters \eqn{\Tau_1}
#' and \eqn{\Tau_2} of the Extended Nelson-Siegel model. These parameters are
#' constant.
#'
#' @return \code{simulateDSFM1D} returns a list containing:
#' \item{\code{dataSim}}{an object of class \code{"DSFM1DData"}, output of the
#' \code{\link{dataDSFM1D}} function. This object can be immediatly used by the
#' \code{\link{DSFM}} algorithm.}
#' \item{\code{YSim}}{the simulated data in a more usual format.}
#' \item{\code{Z_tl}}{the simulated factor loadings.}
#' \item{\code{x1}}{the vector of the covariates.}
#'
#' Depending on the model simulated, the functions returns also:
#' \item{\code{m_l}}{the factors functions used to compute the DSFM.}
#' \item{\code{tau}}{the constant \eqn{\Tau_1} and \eqn{\Tau_2} used to compute
#' the Extended Nelson-Siegel model.}
#'
#' @references Linton, Oliver et al. (2001). "Yield Curve Estimation by Kernel
#' Smoothing Methods". In: \emph{Journal of Econometrics 105.1}, pp. 185-223.
#'
#' Bliss, Robert R.(1997). "Testing Term Structure Estimation Methods".
#' In: \emph{Advances in Futures and Options Research 9}, pp. 197-231.
#'
#' @seealso \code{\link{dataDSFM1D}}, \code{\link{DSFM}}, \code{\link{DSFM1D}}.
#'
#' @export
#'
simulateDSFM1D <- function(model = "ns", n = 100, x1 = 1:30, L = 3,
                           var = 0.0005, beta = c(0.065, -0.015, 0.05),
                           tau=c(0.5, 6)) {

  L0   <- L + 1

  # Simulate DSFM
  if (model =="dsfm") {
    # Z_tl
    Z0 <- rep(1,n)
    Z1 <- arima.sim(list(ar = 0.7, ma = 0), n)
    Z2 <- arima.sim(list(ar = 0.2, ma = 0), n)
    Z3 <- arima.sim(list(ar = 0.2, ma = 0), n)
    Zs <- data.frame(cbind(Z0, Z1, Z2, Z3))
    Zs <- list(Z0, Z1, Z2, Z3)

    # Function m_l
    m <- function(u1 = 0,u2 = 0) {
      I  <- length(u1)
      m0 <- c()
      m1 <- m2 <- m3 <- m0

      for(i in 1:I) {
        m0[i] <- 0
        m1[i] <- 1
        m2[i] <- cos(u1[i] / 8)
        m3[i] <- sin(u1[i] / 8)
      }
      return(list(m0, m1, m2, m3))
    }
    # DSFM
    J <- length(x1)
    YSim <- matrix(0, n, J)
    for (t in 1:n) {
      YTSim <- list()
      for (l in 1:L0) {
        YTSim[[l]] <- m(x1)[[l]] * Zs[[l]][t]
      }
      YTSim <- Reduce("+", YTSim) + rnorm(J, 0, var)
      YSim[t, ] <- YTSim
    }
    namesZs   <- paste0("Z_t", 0:L)
    date      <- seq(as.Date("2000-01-01"), by = "month", length.out = n)
    YSim      <- data.frame(date, YSim)
    Zs        <- data.frame(date, Zs)
    names(Zs) <- c("Date", namesZs)
    dataSim   <- dataDSFM1D(YSim, x1)

    listData  <- list(dataSim = dataSim, YSim = YSim, Z_tl = Zs, m_l = m,
                      x1 = x1)
  }

  # Simulate Nelson-Siegel
  if (model == "ns") {

    # Z_tl
    # VAR(1)
    K <- 3 # Number of variables
    p <- 1 # Lag-order

    # Variance of the Error term - Sigma
    Sigma   <- diag(rep(1e-5, K))
    # A - Parameter Matrix
    A       <- matrix(c(0.98, 0.03, -0.03, -0.04, 0.90, 0.02, 0.24, 0.09, 0.47),
                      K, K, byrow = T)

    # Generate VAR
    Zlag    <- beta # Linton
    Zs      <- Zlag
    mu      <- rep(0, K) # Intercept
    for (i in 1:(n - 1)) {
      temp  <- mu + A %*% Zlag + MASS::mvrnorm(1, rep(0, K), Sigma)
      Zs    <- cbind(Zs, temp)
      Zlag  <- temp
    }
    Zs      <- t(Zs)

    # NS
    J <- length(x1)
    YSim <- matrix(0, n, J)
    for (t in 1:n) {
      YSim[t, ] <- (Zs[t,1] +
                      Zs[t,2] * ((1 - exp(-x1 / tau[1])) / (x1 / tau[1])) +
                      Zs[t,3] * (((1 - exp(-x1 / tau[2])) /
                                    (x1 / tau[2])) - exp(-x1 / tau[2])) +
                      rnorm(J, 0, var))
    }
    namesZs   <- paste0("Z_t", 0:2)
    date      <- seq(as.Date("2000-01-01"), by = "month", length.out = n)
    YSim      <- data.frame(date, YSim)
    Zs        <- data.frame(date, Zs, row.names = NULL)
    names(Zs) <- c("Date", namesZs)
    dataSim   <- dataDSFM1D(YSim, x1)
    listData  <- list(dataSim = dataSim, YSim = YSim, Z_tl = Zs, x1 = x1,
                      tau = tau)
  }

  return(listData)
}


#######
# ----------------------------------------------------------------------------- #
# Prepare Data ---------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
#####

# Format data function -------------------------------------------------------- #

#' DSFM Data Set Generation
#'
#' This function generates a data set to be used in the \code{\link{DSFM}}
#' function.
#'
#' @param y a numeric matrix of the data. The first column is a time indicator
#' with a proper \code{date} format. The remainings columns are the responses
#' \code{y} for each covariates.
#' @param x1 a numeric vector of the covariates.
#'
#' @return \code{dataDSFM1D} returns a list, which belongs to the class
#' \code{"DSFM1DData"}. The list contains the dates, the responses, and
#' the covariates in three distinct columns.
#'
#' The generic functions \code{print},\code{summary}, and \code{plot}
#' are available for this class.
#'
#' @note A data set with the correct \code{DSFM1D} class format can be
#' inserted as \code{y} to return an \code{object} of class \code{"DSFM1DData"}.
#'
#' @seealso \code{\link{summary.DSFM1DData}}, \code{\link{plot.DSFM1DData}}.
#'
#' @export
#'
dataDSFM1D <- function(y, x1=NULL) {

  if (dim(y)[2] == 3) {
    data        <- data.frame(y)
    data[,1]    <- as.Date(data[,1])
  } else {
    x1   <- data.frame(x1)
    date <- y[ ,1]
    if (inherits(date, "Date")) {
      date <- data.frame(as.Date(y[ ,1] %x% rep(1, dim(x1)[1]),
                                 origin = "1970-01-01"))

      # If Dates is not as format Date. Tries to deduce the correct dates
    } else {
      format <- gsub("[[:punct:]]", "-", y[ ,1])
      format <- data.frame(read.table(text = format, sep = "-",
                                      colClasses = "factor"))
      if (max(as.numeric(format[ ,2])) <= 12) {
        names(format)[2] <- "months"
        if (nchar(as.character(format[1,1])) == 4) {
          names(format)[1] <- "years"
          names(format)[3] <- "days"
        } else {
          if (max(as.numeric(format[ ,1])) %in% 28:31) {
            names(format)[1] <- "days"
            names(format)[3] <- "years"
          } else {
            names(format)[3] <- "days"
            names(format)[1] <- "years"
          }
        }
      } else {
        names(format)[2] <- "days"
        if (nchar(as.character(format[1,1])) == 4) {
          names(format)[1] <- "years"
          names(format)[3] <- "months"
        } else {
          if (max(as.numeric(format[ ,1])) <= 12) {
            names(format)[1] <- "months"
            names(format)[3] <- "years"
          } else {
            names(format)[1] <- "years"
            names(format)[3] <- "months"
          }
        }
      }

      date      <- paste(format$years, format$months, format$days, sep = "-")
      date      <- c(t(matrix(rep(date, dim(x1)[1]), ncol = dim(x1)[1])))
      date      <- as.Date(date)
    }
    y           <- data.frame(c(t(y[ ,2:length(y)])))
    data        <- data.frame(date, y, x1)
  }
  names(data)   <- c("Date", "y", "x1")
  class(data)   <- "DSFM1DData"
  data$call     <- match.call()
  return(data)
}

# Print Function -------------------------------------------------------------- #

#' @export
#'
print.DSFM1DData <- function(x, ...) {

  cat("Dataset of class DSFM1DData to be used in the DSFM() function.\n\n")
  cat("Number of observations: ", length(unique(x$Date)), ".\n",
      sep = "")
  cat("Covariates range: from ", min(x$x1), " to ",
      max(x$x1), ".", sep = "", "\n")
  cat("Time range: from ", as.character(min(x$Date)),
      " to ", as.character(max(x$Date)), sep = "", ".\n")
}

# Summary Function ------------------------------------------------------------ #

#' Summarizing One-Dimensional DSFM Data Set
#'
#' \code{summary} method for class \code{"DSFM1DData"}.
#'
#' @param object an object of class \code{"DSFM1DData"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function \code{summary.DSFM1DData} returns summary
#' statistics of the data set given in \code{object}.
#'
#' @seealso \code{\link{dataDSFM1D}}, \code{\link{plot.DSFM1DData}}.
#'
#' @export
#'
summary.DSFM1DData <- function(object, ...) {

  x1Summary <- data.frame(rbind(matrix(unlist(by(object$y,
                                                 object$x1,
                                                 summary),
                                              use.names = F), nrow = 6),
                                by(object$y,
                                   object$x1, sd),
                                by(object$y,
                                   object$x1, skewness),
                                by(object$y,
                                   object$x1, kurtosis)),
                          check.names = F, row.names = NULL)
  row.names(x1Summary) <- c("Min.", "1st Qu.", "Median", "Mean",
                            "3rd Qu.", "Max.", "Std.", "Skew.", "Kurt.")

  cat("Call:\n")
  print(object$call)

  cat("\nSummary of the covariates:\n")
  print(x1Summary)
}

# Plot Function  -------------------------------------------------------------- #

#' Plot Method for an \code{DSFM1DData} Object
#'
#' Plots the 3D vizualisation of the data set contained in \code{object} of class
#' \code{DSFM1DData}.
#'
#' @param x an \code{object} of class \code{"DSFM1DData"}.
#' @param pal the color palette for the fit plot. To choose between
#' \code{"pink"},\code{"blue"},\code{"light"},\code{"dark"}.
#' @param theta angles defining the viewing direction of the fit plot.
#' @param border the color of the line drawn around the surface facets of the
#' fit plot. The default, NULL, corresponds to par("fg"). A value of NA will
#' disable the drawing of borders: this is sometimes useful when the surface is
#' shaded.
#' @param box should the bounding box for the surface be displayed.
#' The default is \code{FALSE}.
#' @param shade the shade at a surface facet is computed as ((1+d)/2)^shade,
#' where d is the dot product of a unit vector normal to the facet and a unit
#' vector in the direction of a light source. Values of shade close to one
#' yield shading similar to a point light source model and values close to zero
#' produce no shading. Values in the range 0.5 to 0.75 provide an approximation
#' to daylight illumination.
#' @param expand a expansion factor applied to the z coordinates. Often used
#' with 0 < expand < 1 to shrink the plotting box in the z direction.
#' @param ticktype character: "simple" draws just an arrow parallel to the axis
#' to indicate direction of increase; "detailed" draws normal ticks as per 2D
#' plots.
#' @param ... other parameters to be passed through to plotting
#' \code{\link{persp}} function.
#'
#' @seealso \code{\link{dataDSFM1D}}, \code{\link{summary.DSFM1DData}}.
#'
#' @export
#'
plot.DSFM1DData <- function(x, pal = "pink", theta = 40, border = NA,
                            box = T, shade = .2, expand = .5,
                            ticktype = "simple",  ...) {

  I     <- length(unique(x$Date))
  J     <- length(unique(x$x1))
  x1    <- seq(1, I) / I
  x2    <- unique(x$x1)
  Y     <- matrix(x$y, I, J, byrow = T)

  switch(pal,
         pink = {
           jetColors <- colorRampPalette(c("#1E0000", "#402525",
                                           "#C27E7E", "#E8E8B4", "#FFFFFF"))
         },
         blue = {
           jetColors <- colorRampPalette(c("#F6EFF7", "#BDC9E1",
                                           "#67A9CF", "#1C9099", "#016C59"))
         },
         light = {
           jetColors <- colorRampPalette(c("#3FB8AF", "#7FC7AF",
                                           "#DAD8A7", "#FF9E9D", "#FF3D7F"))
         },
         dark = {
           jetColors <- colorRampPalette(c("#556270", "#4ECDC4",
                                           "#C7F464", "#FF6B6B", "#C44D58"))
         },
         stop('Invalid pal Parameter. Valid Parameters are: "pink", "blue",
                "light, "dark".')
  )
  color <- jetColors(200)
  nrz <- nrow(Y)
  ncz <- ncol(Y)
  zFacet <- (Y[-1, -1] + Y[-1, -ncz] + Y[-nrz, -1] + Y[-nrz, -ncz]) / 4
  facetCol <- cut(zFacet, 200)
  persp(x1, x2, Y, xlab = "Time", ylab = "x", zlab = "Y",
        col = color[facetCol], main = "Data", theta = theta, border = border,
        box = box, shade = shade, expand = expand, ticktype = ticktype, ...)
}

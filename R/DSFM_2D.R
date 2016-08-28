# ----------------------------------------------------------------------------- #
# --------------- Yield Curve Estimation using DSFM2D ------------------------- #
# ----------------------------------------------------------------------------- #

#######
# ----------------------------------------------------------------------------- #
# DSFM Two-dimensions  -------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
#######

#' Estimation of Dynamic Semiparametric Factor Model for Two-Dimensional Data
#'
#' \code{DSFM2D} performs a model estimation using Dynamic Semiparametric Factor
#' mechanics with two-dimensional covariates. This function is called by the
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
#' \eqn{Y_{t,j}} in second row, and the coordinates \eqn{X_{t,j}} in the two
#' remaining rows. Proper formatting has to be done using the
#' \code{\link{dataDSFM2D}}.
#' @param numDataPoints the number of points in the axis of the grid to perform
#' the kernel density estimation.
#' @param h the bandwidth used to perform the kernel density estimation. can be
#' a single global parameter, a vector of two bandwidths (one by dimension) or a
#' matrix of size \eqn{numDataPoints x 2} for local bandwidth.
#' @param L the number of underlying factors.
#' @param initialLoad the type of initial loadings to choose between White Noise
#' \code{"WN"}, and AR(1) process \code{"AR"}. Required as starting value of the
#' algorithm. Changing the \code{initialLoad} can sligthly improve the
#' convergence rate.
#' @param tol the tolerance for the algorithm to stop.
#' @param maxIt the maximum number of iterations for the algorithm to break.
#'
#' @return \code{DSFM2D} returns an object of class \code{"DSFM2D"}.
#'
#' The generic functions \code{print},\code{summary}, \code{plot} and
#' \code{predict} are available.
#'
#' An object of class \code{"DSFM2D"} is a list containing the following
#' components:
#' \item{\code{Data}}{the input data.}
#' \item{\code{YHat}}{the estimated \eqn{\hat{Y}_{t,j}}.}
#' \item{\code{ZHat}}{the estimated factor loadings \eqn{\hat{Z}_{t,j}}.}
#' \item{\code{mHat}}{the estimated factor functions \eqn{\hat{m}_l}.}
#' \item{\code{EV}}{gives the Explained Variance, used to select the approriate
#' number of factors.}
#' \item{\code{RMSE}}{gives the Root Mean Squared Error, used to compare models.}
#' \item{\code{Bw}}{gives the bandwidth \eqn{h} used and two selection criteria
#'  to select the optimal bandwidth.}
#' \item{\code{Density}}{the kernel density estimation performed.}
#' \item{\code{Convergence}}{the value of the algorithm stopping criterion at
#' each loop.}
#' \item{\code{Time}}{an indicator of the time taken by the function to perform
#' the fit.}
#'
#' @author The implementation of model by Marc Gumowski was based on
#' Fengler and al. (2007).
#'
#' @references Fengler, Matthias R, Wolfgang K Haerdle, and Enno Mammen (2007).
#' "A Semiparametric Factor Model for Implied Volatility Surface Dynamics".
#' In: \emph{Journal of Financial Econometrics 5.2}, pp. 189-218.
#'
#' Haerdle, Wolfgang K., and Piotr Majer (2014).
#' "Yield Curve Modeling and Forecasting using Semiparametric Factor Dynamics".
#' In: \emph{The European Journal of Finance}, pp. 1-21.
#'
#' Borak, Szymon, Matthias R. Fengler, and Wolfgang K. Haerdle (2005)."DSFM
#' Fitting of Implied Volatility Surfaces". In: \emph{5th International
#' Conference on Intelligent Systems Design and Applications (ISDA'05)},
#' pp. 526-531. IEEE.
#'
#' @seealso \code{\link{summary.DSFM2D}} for
#' summaries and \code{\link{plot.DSFM2D}} for plot
#' possibilities.
#'
#' \code{\link{predict.DSFM2D}} provide succint
#' predictions.
#'
#' \code{\link{dataDSFM2D}} has to be used before
#' using the \code{\link{DSFM}} function to ensure that the data are correctly
#' formated.
#'
#' \code{\link{simulateDSFM2D}} is a function to simulate data that can be used
#' as simple example purposes.
#'
#' @examples
#' \dontrun{
#' # Prepare the data --------------------------------------------------------- #
#' simulatedData <- simulateDSFM2D()
#' dsfmData      <- simulatedData$dataSim
#'
#' # Set the parameters ------------------------------------------------------- #
#' h        <- c(0.05,0.05)
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
#' }
#'
DSFM2D <- function(data, numDataPoints = 25, h = c(0.5, 0.5), L = 3,
                   initialLoad = "WN", tol = 1e-5, maxIt = 301) {

  Time1 <- Sys.time()                                 # Get starting time

  # Initial Settings
  date  <- unique(data$Date)                          # Get the dates
  I     <- length(date)                               # T
  L0    <- L + 1                                      # L with factor L_0
  y     <- list()
  x1    <- x2 <- J <- y
  for (t in 1:I){
    x1[[t]] <- unique(data$x1[which(data$Date == date[t])])
    x2[[t]] <- unique(data$x2[which(data$Date == date[t])])
    J[[t]]  <- length(x1[[t]]) * length(x2[[t]])
    y[[t]]  <- matrix(data$y[which(data$Date == date[t])], length(x1[[t]]),
                      length(x2[[t]]), byrow = T) # Get the y
  }
  # Normalization of the covariates
  x1 <- lapply(x1, function(x) x / max(x))
  x2 <- lapply(x2, function(x) x / max(x))

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
  ZHat     <- cbind(rep(1, I), ZHat)

  # Create a regular grid of points u covering the whole space
  minx1    <- min(unlist(x1))
  minx2    <- min(unlist(x2))
  maxx1    <- max(unlist(x1))
  maxx2    <- max(unlist(x2))
  delta    <- ((maxx1 - minx1) * (maxx2 - minx2)) / numDataPoints^2
  u1       <- u1Points <-seq(minx1, maxx1, length.out = numDataPoints)
  u2       <- u2Points <-seq(minx2, maxx2, length.out = numDataPoints)
  u1Points <- rep(u1Points, length(u2Points))
  u1Points <- sort(u1Points)
  u        <- data.frame(u1Points, u2Points)
  U        <- length(u1) * length(u2)

  # Loop parameters
  ZHatOld    <- matrix(0, I, L0)
  mHatOld    <- matrix(0, U, L0)
  it         <- 0
  stopCriterion     <- 1
  plotStopCriterion <- c()

  # KERNEL -------------------------------------------------------------------- #

  if (length(h) == 1) {
    h <- matrix(h, numDataPoints * numDataPoints, 2)
  }
  if (length(h) == 2){
    h <- matrix(h, numDataPoints * numDataPoints, 2, byrow = T)
  }

  Kernel <- KernelDensity2D(y, I, J, x1, x2, u, U, h)

  # LOOP ---------------------------------------------------------------------- #

  while(stopCriterion >= tol) {

    it = it + 1

    if (it == maxIt) {
      stop("Too many iterations, no convergence")
    }

    # Matrix B
    Bu <- matrix(0, L0, L0)
    BList <- replicate(U, list(Bu))
    for (n in 1:U) {
      BList[[n]] <- t(ZHat) %*% (ZHat * Kernel$jPHat[ ,n])
    }

    # Vector Q
    Q <- vector(mode="numeric", L0)
    QList <- replicate(U, list(Q))
    for (n in 1:U) {
      QList[[n]] <- t(ZHat) %*% Kernel$jQHat[ ,n]
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
    for (t in 1:I){
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
    stopCriterionAll <- vector(mode="numeric", I)
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
  mHat[ ,2:L0] <- t((G %^% (-0.5)) %*% t(mHat[ ,2:L0]))

  # ZHat new
  for (t in 1:I) {
    ZHat[t,2:L0] <- t(G %^% (0.5) %*% (ZHat[t,2:L0] + solve(G) %*% g))
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

  # Fit the model ------------------------------------------------------------- #

  Y    <- data$y
  YBar <- mean(Y)
  YHat <- list()
  # Interpolation of function m_l
  for (t in 1:I) {
    m_l <- matrix(0, I, J[[t]])
    m_l <- replicate(L0, list(m_l))
    for (l in 1:L) {
      m_l[[l]] <- akima::interp(u[ ,2], u[ ,1], mHat[ ,l], xo = x1[[t]],
                                yo = x2[[t]], linear = T, extrap = F)$z
    }
    # Explained Variance - Model Size Selection --------- #
    YTHat <- c()
    YTHat <- list(L0, YTHat)
    for (l in 1:L) {
      YTHat[[l]] <- ZHat[t,l] * m_l[[l]]
    }
    YHat[[t]] <- Reduce("+", YTHat)
  }
  YHatRMSE  <- YHat
  YHat  <- unlist(YHat)

  numerator <- (Y - YHat)^2
  numerator <- sum(numerator)
  denominator <- (Y - YBar)^2
  denominator <- sum(denominator)
  EV <- 1 - numerator / denominator
  rownames(EV) <- NULL

  # Goodness-of-fit - Root Mean Squared Error - RMSE ---- #
  numeratorRMSE <- list()
  for (t in 1:I) {
    numeratorRMSE[[t]] <- (1 / J[[t]]) * (y[[t]] - YHatRMSE[[t]])
  }
  numeratorRMSE <- Reduce("+", numeratorRMSE)
  numeratorRMSE <- sum(numeratorRMSE^2)
  RMSE <- sqrt((1 / I) * numeratorRMSE)


  # Bandwidth Selection  -------------------------------- #
  if ((length(unique(h[ ,1])) == 1) & (length(unique(h[ ,2])) == 1)) {
    w     <- 1 / pHat
    N     <- length(Y)
    Kh    <- (1 / (h[1,1] * h[1,2])) * NormalKernel1D(0 / h[1,1])  %*%
      t(NormalKernel1D(0 / h[1,2]))
    # Weighted AIC_2
    AIC2 <- (1 / N) * numerator * exp(2 * (L / N) * Kh * sum(w) /  sum(w*pHat))
    # Weighted SC_1
    SC1  <- (1 / N) * numerator * exp(log(N) * (L / N) * Kh * sum(w) /
                                        sum(w*pHat))
  }else{
    AIC2 <- NA
    SC1  <- NA
  }

  # Outputs ------------------------------------------------------------------- #

  # Create data frames as outputs
  ZHat           <- data.frame(date, ZHat)
  names(ZHat)    <- c("Date", paste0("Z_t", 0:L, ".hat"))
  mHat           <- data.frame(u * c(max(data$x1),max(data$x2)), mHat)
  names(mHat)    <- c("u1", "u2", paste0("m_", 0:L, ".hat"))
  YHat           <- data.frame(data$Date, YHat, data$x1, data$x2)
  names(YHat)    <- c("Date", "y.hat", "x1", "x2")
  EV             <- data.frame(EV)
  names(EV)      <- paste0("EV(L = ", L, ")")
  names(RMSE)    <- "RMSE"
  Bw             <- data.frame(t(c(unique(h), AIC2, SC1)))
  names(Bw)      <- c("h1", "h2", "wAIC_2", "wSC_1")

  Time2 <- Sys.time() - Time1 # Calculate time difference

  model <- list(Data = data, YHat = YHat, ZHat = ZHat, mHat = mHat, EV = EV,
                RMSE = RMSE, Bw = Bw, Density = pHat,
                Convergence = plotStopCriterion, Iterations = it, Time = Time2)
  model$call   <- match.call()
  class(model) <- "DSFM2D"

  return(model)
}



#######
# ----------------------------------------------------------------------------- #
# S3 Methods (print.DSFM2D, summary.DSFM2D, plot.DSFM2D, predict.DSFM2D) ------ #
# ----------------------------------------------------------------------------- #
#######

# Print Function -------------------------------------------------------------- #

print.DSFM2D <- function(object, ...){

  cat("Call:\n")
  print(object$call)
  cat("\nY:\n")
  print(object$Y)
  cat("\nEstimated Y:\n")
  print(object$YHat)
  cat("\nEstimated Z:\n")
  print(object$ZHat)
  cat("\nEstimated m:\n")
  print(object$mHat)
  cat("\nExplained Variance:\n")
  print(object$EV, row.names = F)
  cat("\nRoot Mean Squared Error:\n")
  print(object$RMSE, row.names = F)
  cat("\nBandwidth:\n")
  print(object$Bw, row.names = F)
  cat("\nIterations:",object$Iterations,"\n")
  cat("\n")
  print(object$Time)
}

# Summary Function ------------------------------------------------------------ #

#' Summarizing DSFM Fits for Two-Dimensional Data
#'
#' \code{summary} method for class \code{"DSFM2D"}.
#'
#' @param object an object of class \code{"DSFM2D"}.
#'
#' @return The function \code{summary.DSFM2D} returns a list of summary
#' statistics of the fitted DSFM given in \code{object}, using the components
#' (list elements) \code{"call"} from its arguments, plus:
#' \item{\code{ZHat}}{a summary of the estimated factor loadings.}
#' \item{\code{EV}}{the Explained Variance, used to select the approriate
#' number of factors.}
#' \item{\code{RMSE}}{the Root Mean Squared Error, used to compare the goodness-
#' of-fit between models.}
#' \item{\code{Bw}}{The bandwidth and its selection criteria.}
#'
#' @seealso \code{\link{dataDSFM2D}}, \code{\link{DSFM}}, \code{\link{DSFM2D}},
#' \code{\link{plot.DSFM2D}}, \code{\link{predict.DSFM2D}}.
#'
summary.DSFM2D <- function(object) {

  # INPUTS -------------------------------------------------------------------- #
  #
  # object : An object of class DSFM2D, output of function DSFM2D().
  #
  # --------------------------------------------------------------------------- #

  ZHat   <- object$ZHat[ ,2:length(object$ZHat)]
  L0     <- length(ZHat)
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
  cat(dim(object$Y)[1],"\n")

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

#' Plot Diagnostics for a DSFM2D Object
#'
#' Four plots (selectable by \code{which}) are currently available: a plot of
#' the convergence rate of the algorithm, a plot of the time series of the
#' estimated loadings for each factor, a plot of the estimated factor functions,
#' and a plot of the fit for a given time. By default, the four plots are
#' provided.
#'
#' @param object an object of class \code{"DSFM1D"}.
#' @param which to choose between \code{"all"},\code{"loadings"},
#' #' \code{"functions"},\code{"convergence"}, and \code{"fit"}.
#' @param n number of time indicator to be plotted.
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
#' @seealso \code{\link{dataDSFM2D}}, \code{\link{DSFM}}, \code{\link{DSFM2D}},
#' \code{\link{summary.DSFM2D}}, \code{\link{predict.DSFM2D}}.
#'
plot.DSFM2D <- function(object, which = "all", n = 1, ask = TRUE, pal = "pink",
                        col = "#016C59", type = "l", theta = 40, border = NA,
                        box = T, shade = .2, expand = .5, ticktype = "detailed",
                        ...) {

  date  <- unique(object$Data$Date)
  L     <- dim(object$ZHat)[2] - 2
  L0    <- L + 1
  u1    <- unique(object$mHat[ ,1])
  u2    <- unique(object$mHat[ ,2])

  par(ask = ask)

  if ("all" %in% which) {
    which = c("convergence", "loadings", "functions", "fit")
  }

  # Convergence Plot
  if ("convergence" %in% which) {
    layout(matrix(1))
    N <- length(object$Convergence)
    x <- 1:N
    plot(x[(N - N * 0.95):N], object$Convergence[(N - N * 0.95):N],
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
      plot(as.Date(date),object$ZHat[ ,l + 2],
           main = bquote(hat(Z)[t*","*.(l)]),
           ylab = "", xlab = "Time", xaxt ="n",
           col = col, type = type, ...)
      axis.Date(1, date, at = labDates)
    }
    layout(matrix(1))
  }
  # Factor Functions
  if ("functions" %in% which) {
    layout(matrix(1:(L0 + L0 %% 2),L0 %% 2 + L0 %/% 2, byrow = T))
    mHat <- list()
    for (l in 1:L0) {
      mHat[[l]] <- matrix(object$mHat[ ,l + 2], length(u1), length(u2),
                          byrow = T)
      persp(u1, u2, mHat[[l]], main = bquote(hat(m)[.(l-1)]), zlab="",
            theta = theta, border = border, box = box, shade = shade,
            expand = expand, ticktype = ticktype, ...)
    }
    layout(matrix(1))
  }
  if ("fit" %in% which) {

    x1    <- unique(object$Data$x1[which(object$Data$Date == date[n])])
    x2    <- unique(object$Data$x2[which(object$Data$Date == date[n])])
    J     <- length(x1) * length(x2)
    Y     <- matrix(object$Data$y[which(object$Data$Date == date[n])],
                    length(x1), length(x2), byrow = T) # Get the y
    YHat  <- matrix(object$YHat[which(object$YHat[ ,1] == date[n]), 2],
                    length(x1), length(x2), byrow = T)
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
    persp(x1, x2, Y, xlab = names(object$Data)[3], ylab = names(object$Data)[4],
          zlab = names(object$Data)[2], col = color[facetCol],
          main = paste("Y -", date[n]), theta = theta, border = border,
          box = box, shade = shade, expand = expand, ticktype = ticktype, ...)

    nrz <- nrow(YHat)
    ncz <- ncol(YHat)
    zFacet <- (YHat[-1, -1] + YHat[-1, -ncz] + YHat[-nrz, -1] +
                 YHat[-nrz, -ncz]) / 4
    facetCol <- cut(zFacet, 200)
    persp(x1, x2, YHat, xlab = names(object$Data)[3],
          ylab = names(object$Data)[4], zlab = names(object$Data)[2],
          col = color[facetCol],
          main = bquote(hat(Y) ~"-"~.(as.character(date)[n])), theta = theta,
          border = border, box = box, shade = shade, expand = expand,
          ticktype = ticktype, ...)
    layout(matrix(1))
  }

  layout(matrix(1))
  par(ask = FALSE)
}

# Predict Function  ----------------------------------------------------------- #

#' Predict Method for Two-Dimensional DSFM Fits
#'
#' Predicted values based on \code{"DSFM2D"} object using Vector Autoregressive
#' Processes (VAR(p)).
#'
#' This function makes uses of package \code{vars} to fit and predict VAR(p)
#' processes, before pluging in the predicted factor loadings in the DSFM. The
#' factors functions are interpolated by the function \code{\link{interp}} of
#' package \code{akima}.
#'
#' @param object an object of class \code{"DSFM1D"}.
#' @param nAhead the number of steps ahead for which prediction is required.
#' @param x1Forecast the vector of covariates in the first dimension to be
#'  forecasted.
#' @param x2Forecast the vector of covariates in the second dimension to be
#' forecasted.
#' @param p the order of the Vector Autoregressive Process to be fitted.
#' @param ... other parameters to be passed through the \code{\link{VAR}}
#' and \code{\link{predict.varest}} functions.
#'
#' @return \code{predict.DSFM2D} returns an object of class
#' \code{"predict.DSFM2D"}. This class is a list containing:
#' \item{\code{YHatForecast}}{the forecasted responses.}
#' \item{\code{YHatForecastList}}{the forecasted responses in a more usual
#' format.}
#' \item{\code{ZHatForecast}}{the predicted factors loadings.}
#' \item{\code{nAhead}}{the number of steps ahead.}
#'
#' @seealso \code{\link{VAR}},\code{\link{predict.varest}},
#' \code{\link{interp}}, \code{\link{DSFM}}, \code{\link{DSFM2D}},
#' \code{\link{dataDSFM2D}}.
#'
#' @references Bernhard Pfaff (2008). VAR, SVAR and SVEC Models: Implementation
#' Within R Package vars. In: \emph{Journal of Statistical Software 27(4)}.
#' URL http://www.jstatsoft.org/v27/i04/.
#'
#' Hiroshi Akima and Albrecht Gebhardt (2015). akima: Interpolation of
#' Irregularly and Regularly Spaced Data. R package version 0.5-12.
#' URL http://CRAN.R-project.org/package=akima
#'
predict.DSFM2D <- function(object, nAhead = 12,
                           x1Forecast = unique(object$mHat[ ,1]),
                           x2Forecast = unique(object$mHat[ ,2]), p = 1, ...) {

  date  <- unique(object$Data$Date)
  if (date[2] - date[1] <= 1) {
    timeDiff <- "day"
  }else if (date[2] - date[1] == 7) {
    timeDiff <- "week"
  }else if (date[2] - date[1] <= 31 & date[2] - date[1] >= 28){
    timeDiff <- "month"
  }else if (date[2] - date[1] <= 123 & date[2] - date[1] >= 118){
    timeDiff <- "quarter"
  }else{
    timeDiff <- "year"
  }
  dateForecast  <- seq.Date(date[length(date)],
                            length.out = nAhead + 1,
                            by = timeDiff)[2:(nAhead + 1)]

  # Initial Parameters
  ZHat  <- object$ZHat[ ,3:length(object$ZHat)]
  mHat  <- object$mHat[ ,3:length(object$mHat)]
  L     <- dim(ZHat)[2]
  L0    <- L + 1
  u       <- data.frame(object$mHat[ ,1], object$mHat[ ,2])

  # VAR(p) estimation of the loadings Z_tl^forecast using vars
  ZVarEst         <- vars::VAR(ZHat, p = p, ...)
  ZVarEstForecast <- predict(ZVarEst, n.ahead = nAhead, ...)
  ZVarEstForecast <- lapply(ZVarEstForecast$fcst, function(x) x[ ,1])
  ZVarEstForecast <- do.call(cbind, ZVarEstForecast)
  ZHatForecast    <- cbind(1, ZVarEstForecast)

  YHatForecastList <- list()
  # Interpolation of function m_l
  for (t in 1:nAhead) {
    m_l <- matrix(0, length(x1Forecast), length(x2Forecast))
    m_l <- replicate(L0, list(m_l))
    for (l in 1:L0) {
      m_l[[l]] <- akima::interp(u[ ,2], u[ ,1], mHat[ ,l],
                                xo = x1Forecast, yo = x2Forecast,
                                linear = T, extrap = F)$z
    }
    # YHat_hj
    YTHat <- c()
    YTHat <- list(L0,YTHat)
    for (l in 1:L0){
      YTHat[[l]] <- ZHatForecast[t,l] * m_l[[l]]
    }
    YHatForecastList[[t]]  <- Reduce("+", YTHat)
  }

  ZHatForecast            <- data.frame(dateForecast, ts(ZHatForecast))
  names(ZHatForecast)     <- c("Date", paste0("Z_t", 0:L, ".hat"))
  names(YHatForecastList) <- dateForecast
  YHatForecast            <- dataDSFM2D(YHatForecastList,x1Forecast,x2Forecast)

  predict <- list(YHatForecast = YHatForecast,
                  YHatForecastList = YHatForecastList,
                  ZHatForecast = ZHatForecast,
                  nAhead = nAhead)
  predict$call   <- match.call()
  class(predict) <- "predict.DSFM2D"

  return(predict)
}

# Print Function -------------------------------------------------------------- #

print.predict.DSFM2D <- function(object, ...) {

  cat("Call:\n")
  print(object$call)
  cat("\nForecasted Y:\n")
  print(object$YHatForecastList)
  cat("\nForecasted Z:\n")
  print(object$ZHatForecast)
  cat("\nHorizon:",object$nAhead)

}






#######
# ----------------------------------------------------------------------------- #
# Simulation ------------------------------------------------------------------ #
# ----------------------------------------------------------------------------- #
#####

# 2 Dimensions Simulation ----------------------------------------------------- #

#' Simulate Responses for Two-Dimensional DSFM
#'
#' This function simulates responses for two-dimensional DSFM.
#'
#' This function is used for example purpose, only few parameters are available
#' to control the simulation. The factors loadings are generated using three
#' independant AR(1) process and the factors functions are predefined to be
#' orthogonals.
#'
#' @param n the number of observations.
#' @param x1 a vector of the first dimension of covariates. They will be constant
#'  through time
#' @param x2 a vector of the second dimension of covariates. They will be constant
#'  through time
#' @param L the number of factors for the DSFM.
#' @param var the error \eqn{\varepsilon_{t,j}} of the models. Allows to control
#'  the noise.
#'
#' @return \code{simulateDSFM2D} returns a list containing:
#' \item{\code{dataSim}}{an object of class \code{"DSFM2DData"}, output of the
#' \code{\link{dataDSFM2D}} function. This object can be immediatly used by the
#' \code{\link{DSFM}} algorithm.}
#' \item{\code{YSim}}{the simulated data in a more usual format.}
#' \item{\code{Z_tl}}{the simulated factor loadings.}
#' \item{\code{m_l}}{the factors functions used to compute the DSFM.}
#'
#' @seealso \code{\link{dataDSFM2D}}, \code{\link{DSFM}}, \code{\link{DSFM2D}}.
#'
simulateDSFM2D <- function(n = 10 , x1 = 1:25 / 25, x2 = 1:25 / 25, L = 3,
                           var = 0.05) {

  L0 <- L + 1

  # Z_tl
  Z0 <- rep(1, n)
  Z1 <- arima.sim(list(ar = 0.7, ma = 0), n) # 0.7
  Z2 <- arima.sim(list(ar = 0.2, ma = 0), n) # 0.2
  Z3 <- arima.sim(list(ar = 0.2, ma = 0), n)
  Zs <- data.frame(cbind(Z0, Z1, Z2, Z3))
  Zs <- list(Z0, Z1, Z2, Z3)

  # Function m_l. Park 2009, p288.
  m <- function(u1 = 0,u2 = 0) {
    m0 <- matrix(0, length(u1), length(u2))
    m1 <- m2 <- m3 <- m0

    for(i in 1:length(u1)){
      for (j in 1:length(u2)){
        # Park (2009)
        m0[i,j] <- 1
        m1[i,j] <- 3.46 * u1[i] - 0.5
        m2[i,j] <- 9.45 * ((u1[i] - 0.5)^2 + (u2[j] - 0.5)^2) - 1.6
        m3[i,j] <- 1.41 * sin(2 * pi * u2[j])
        #         # Borak (2007)
        #         m0[i,j] <- 0
        #         m1[i,j] <- 1
        #         m2[i,j] <- -5 * u1[i] + 5
        #         m3[i,j] <- -2 * u2[j] + 1
      }
    }
    return(list(m0, m1, m2, m3))
  }

  # Simulate DSFM
  YSim <- list()
  for (t in 1:n) {
    YTSim <- list()
    for (l in 1:L0) {
      YTSim[[l]] <- sweep(m(x1, x2)[[l]], 1, Zs[[l]][t], "*")
    }
    YTSim <- Reduce("+", YTSim) + matrix(rnorm(n * length(x2), 0, var),
                                         length(x1), length(x2))
    YSim[[t]] <- YTSim
  }

  date        <- seq(as.Date("2000/1/1"), by = "month", length.out = n)
  names(YSim) <- date
  Zs          <- data.frame(date, Zs)
  names(Zs)   <- c("Date", paste0("Z_t", 0:L))

  dataSim     <- dataDSFM2D(YSim,x1,x2)
  listData    <- list(dataSim = dataSim, YSim = YSim, Z_tl = Zs, m_l = m)

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
#' @param y a numeric list of matrix containing the response variables.
#' Each matrix represents one date.
#' @param x1 a numeric vector of the covariates x1.
#' @param x2 a numeric vector of the covariates x2.
#'
#' @return \code{dataDSFM2D} returns a list, which belongs to the class
#' \code{"DSFM2DData"}. The list contains the dates, the responses and
#' the covariates in four distinct columns.
#'
#' The generic functions \code{print},\code{summary}, and \code{plot}
#' are available for this class.
#'
#' @note If the data set contains different covariates for each dates, the data
#' has to be formatted manually then inserted in the \code{dataDSFM2D}
#' function to obtain an object of class \code{DSFM2DData}.
#'
#' The correct format is a four columns \code{data.frame} containing a time
#' indicator as first columns, the responses in the second columns and the
#' coordinates of the two dimensions in the two remaining columns.

#' @seealso \code{\link{summary.DSFM2DData}}, \code{\link{plot.DSFM2DData}}.
#'
dataDSFM2D <- function(y,x1=NULL,x2=NULL) {

  # INPUTS -------------------------------------------------------------------- #
  #          y : Insert your y in a list of matrix format. y is a list.
  #              1 matrix by day.
  #
  #         x1 : Vector of the coordinates x1.
  #         x2 : Vector of the coordinates x2.
  # --------------------------------------------------------------------------- #

  if (is.list(y) == FALSE) {
    if (dim(y)[2] == 4) {
      data        <- data.frame(y)
      data[,1]    <- as.Date(data[,1])
      names(data) <- c("Date", "y", "x1", "x2")
    }
  } else {
    x  <- data.frame(x1 %x% rep(1, length(x2)), x2)
    data <- list()
    for (i in 1:length(y)) {
      data[[i]] <- data.frame(as.Date(names(y))[i], c(t(y[[i]])), x)
      names(data[[i]]) <- c("Date", "y", "x1", "x2")
    }
    data <- do.call("rbind", data)
  }
  class(data)    <- "DSFM2DData"
  data$call      <- match.call()
  return(data)
}

# Print Function -------------------------------------------------------------- #

print.DSFM2DData <- function(object, ...) {

  cat("Dataset of class DSFM2DData to be used in the DSFM() function.\n\n")
  cat("Number of observations: ", length(unique(object$Date)), ".\n",
      sep = "")
  cat("Covariates range:\n")
  cat("          - x1: from ", min(object$x1), " to ",
      max(object$x1), ".", sep = "", "\n")
  cat("          - x2: from ", min(object$x2), " to ",
      max(object$x2), ".", sep = "", "\n")
  cat("Time range: from ", as.character(object$Date[1]),
      " to ", as.character(object$Date[length(object$Date)]), sep = "", ".\n")
}

# Summary Function ------------------------------------------------------------ #

#' Summarizing Two-Dimensional DSFM Data Set
#'
#' \code{summary} method for class \code{"DSFM2DData"}.
#'
#' @param object an object of class \code{"DSFM2DData"}.
#'
#' @return The function \code{summary.DSFM2DData} returns a summary
#' statistics of the data set given in \code{object} for each covariates
#' directions.
#'
#' @seealso \code{\link{dataDSFM2D}}, \code{\link{plot.DSFM2DData}}.
#'
summary.DSFM2DData <- function(object) {

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

  x2Summary <- data.frame(rbind(matrix(unlist(by(object$y,
                                                 object$x2,
                                                 summary),
                                              use.names = F), nrow = 6),
                                by(object$y,
                                   object$x2, sd),
                                by(object$y,
                                   object$x2, skewness),
                                by(object$y,
                                   object$x2, kurtosis)),
                          check.names = F, row.names = NULL)
  row.names(x2Summary) <- c("Min.", "1st Qu.", "Median", "Mean",
                            "3rd Qu.", "Max.", "Std.", "Skew.", "Kurt.")

  cat("Call:\n")
  print(object$call)

  cat("\nSummary of the covariates x1:\n")
  print(x1Summary)

  cat("\nSummary of the covariates x2:\n")
  print(x2Summary)
}

# Plot Function  -------------------------------------------------------------- #

#' Plot Method for an \code{DSFM2DData} Object
#'
#' Plots the 3D vizualisation of the data set contained in \code{object} of class
#' \code{DSFM2DData}.
#'
#' @param object an object of class \code{"DSFM2DData"}.
#' @param n number of time indicator to be plotted.
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
#' @seealso \code{\link{dataDSFM2D}}, \code{\link{summary.DSFM2DData}}.
#'
plot.DSFM2DData <- function(object, n = 1, pal = "pink", theta = 40,
                            border = NA, box = T, shade = .2, expand = .5,
                            ticktype = "simple", ...) {

  x1    <- unique(object$x1[which(object$Date == unique(object$Date)[n])])
  x2    <- unique(object$x2[which(object$Date == unique(object$Date)[n])])
  J     <- length(x1) * length(x2)
  Y     <- matrix(object$y[which(object$Date == unique(object$Date)[n])],
                  length(x1), length(x2), byrow = T)

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
  persp(x1, x2, Y, xlab = "x1", ylab = "x2",
        zlab = "Y", col = color[facetCol],
        main = paste("Data -", unique(object$Date)[n]), theta = theta,
        border = border, box = box, shade = shade, expand = expand,
        ticktype = ticktype, ...)
}

# ----------------------------------------------------------------------------- #
# --------------- Yield Curve Estimation using DSFM --------------------------- #
# ----------------------------------------------------------------------------- #




#######
# ----------------------------------------------------------------------------- #
# Main Algorithm  ------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
#######


# DSFM Estimation ------------------------------------------------------------- #

#' Estimation Using Dynamic Semiparametric Factor Model
#'
#' \code{DSFM} performs a model estimation using Dynamic Semiparametric Factor
#' mechanics.
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
#' Secondly, the bandwidth \eqn{h} is basically set to 0.05 but optimal bandwidth
#' has to be found externally. The algorithm always normalize the covariates
#' to work on an estimation grid bounded beetween [0,1].
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
#' remaining rows. Proper formatting has to be done using the
#'  \code{\link{dataDSFM1D}} or \code{\link{dataDSFM2D}} functions.
#' @param numDataPoints the number of points in the axis of the grid to perform
#' the kernel density estimation.
#' @param h the bandwidth used to perform the kernel density estimation. In one
#' dimension, can be either a single global parameter, or a vector of the same
#' length of numDataPoints to perform local kernel estimation. In two dimension,
#' can be a single global parameter, a vector of two bandwidths (one by
#' dimension) or a matrix of size \eqn{numDataPoints x 2} for local bandwidth.
#' @param L the number of underlying factors.
#' @param initialLoad the type of initial loadings to choose between White Noise
#' \code{"WN"}, and AR(1) process \code{"AR"}. Required as starting value of the
#' algorithm. Changing the \code{initialLoad} can sligthly improve the
#' convergence rate.
#' @param tol the tolerance for the algorithm to stop.
#' @param maxIt the maximum number of iterations for the algorithm to break.
#'
#' @return \code{DSFM} returns an object of class \code{"DSFM1D"} or
#' \code{"DSFM2D"} depending on the dimension of the input data.
#'
#' The generic functions \code{print},\code{summary}, \code{plot} and
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
#' \item{\code{RMSE}}{gives the Root Mean Squared Error, used to compare the
#' goodness-of-fit between models.}
#' \item{\code{Bw}}{gives the bandwidth \eqn{h} used and two selection criteria to
#' select the optimal bandwidth.}
#' \item{\code{x1}}{the vector of the covariates.}
#' \item{\code{Density}}{the kernel density estimation performed.}
#' \item{\code{Convergence}}{the value of the algorithm stopping criterion at
#' each loop.}
#' \item{\code{Time}}{an indicator of the time taken by the function to perform
#' the fit.}
#'
#' Object of class \code{"DSFM2D"} provides the same outputs except that the
#' \code{Y} and \code{YHat} are kept in the specific format used by the function.
#'
#' @author The implementation of model by Marc Gumowski was based on
#' Fengler and al. (2007).
#'
#' @references Borak, Szymon, Matthias R. Fengler, and Wolfgang K. Haerdle (2005)."DSFM
#' Fitting of Implied Volatility Surfaces". In: \emph{5th International
#' Conference on Intelligent Systems Design and Applications (ISDA'05)},
#' pp. 526-531. IEEE.
#'
#' Fengler, Matthias R, Wolfgang K Haerdle, and Enno Mammen (2007).
#' "A Semiparametric Factor Model for Implied Volatility Surface Dynamics".
#' In: \emph{Journal of Financial Econometrics 5.2}, pp. 189-218.
#'
#' Haerdle, Wolfgang K., and Piotr Majer (2014).
#' "Yield Curve Modeling and Forecasting using Semiparametric Factor Dynamics".
#' In: \emph{The European Journal of Finance}, pp. 1-21.
#'
#' @seealso \code{\link{summary.DSFM1D}} / \code{\link{summary.DSFM2D}} for
#' summaries and \code{\link{plot.DSFM1D}} / \code{\link{plot.DSFM2D}} for plot
#' possibilities.
#'
#' \code{\link{predict.DSFM1D}} / \code{\link{predict.DSFM2D}} provide succint
#' predictions.
#'
#' \code{\link{dataDSFM1D}} / \code{\link{dataDSFM2D}} have to be used before
#' using the \code{\link{DSFM}} function to ensure that the data are correctly
#' formated.
#'
#' \code{\link{simulateDSFM1D}} / \code{\link{simulateDSFM2D}} are functions
#' to simulate data that can be used as simple example purposes.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis.Date layout par persp plot
#' @importFrom stats approx arima.sim rnorm sd ts
#' @importFrom utils read.table tail
#'
#' @examples
#'
#' ## One-Dimensional Data #################################################### #
#'
#' # Prepare the data --------------------------------------------------------- #
#' # Interest rate of zero-coupon bond yield curves. Data from Bank of Canada.
#' data(canadianYieldCurves)
#' maturity <- c(1/4,1/2,3/4,1:10,20,30)
#' dsfmData <- dataDSFM1D(canadianYieldCurves,maturity)
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
#'
#'
#' ## Two-Dimensional Data #################################################### #
#'
#'#' # Prepare the data --------------------------------------------------------- #
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
#'
#' @export
#'
DSFM <- function(data, numDataPoints = 25, h = 0.05, L = 3,
                         initialLoad = "WN", tol = 1e-5, maxIt = 301) {

  if (class(data) == "DSFM1DData") {
    DSFM1D(data, numDataPoints, h, L, initialLoad, tol, maxIt)
  } else if (class(data) == "DSFM2DData") {
    DSFM2D(data, numDataPoints, h, L, initialLoad, tol, maxIt)
  } else {
    stop(
"Invalid Class Object. Object must be of class 'DSFM1DData' or 'DSFM2DData'")
  }
}




#######
# ----------------------------------------------------------------------------- #
# Inner Functions ------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
#######

# Multivariate Random Normal Distribution ------------------------------------- #
# Dependency, need to import MANDATORY
# require(MASS)

# Bilinear Interpolation ------------------------------------------------------ #
# Dependency, need to import MANDATORY
# require(akima)

# Estimate VAR(p) ------------------------------------------------------------- #
# Dependency, need to import MANDATORY
# require(vars)

# Power of matrix operator ---------------------------------------------------- #

#' Power of matrix operator
#'
#' \code{\%^\%} computes the matrix x at the power n. Inner function of the
#' DSFM iterative algorithm.
#'
#' @param x a matrix.
#' @param n the power.
#'
#' @export
#'
"%^%" <- function(x,n) {
  with(eigen(x), vectors %*% (values^n * t(vectors)))
}

# Summary Statistics (from package "moments" - 05/07/16) ---------------------- #

# Skewness -------------------------------------------------------------------- #

#' Skewness
#'
#' This function computes the skewness of given data.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param na.rm logical. Should missing values be removed?
#'
#' @author Lukasz Komsta
#'
#' @references Lukasz Komsta and Frederick Novomestky (2015). moments: Moments,
#' cumulants, skewness, kurtosis and related tests. R package
#' version 0.14. http://CRAN.R-project.org/package=moments
#'
#' @export
#'
skewness <- function(x, na.rm = FALSE) {
  if (is.matrix(x))
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x))
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

# Kurtosis -------------------------------------------------------------------- #

#' Kurtosis
#'
#' This function computes the estimator of Pearson's measure of kurtosis.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param na.rm logical. Should missing values be removed?
#'
#' @author Lukasz Komsta
#'
#' @references Lukasz Komsta and Frederick Novomestky (2015). moments: Moments,
#' cumulants, skewness, kurtosis and related tests. R package
#' version 0.14. http://CRAN.R-project.org/package=moments
#'
#' @export
#'

kurtosis <- function(x, na.rm = FALSE) {
  if (is.matrix(x))
    apply(x, 2, kurtosis, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    n <- length(x)
    n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
  }
  else if (is.data.frame(x))
    sapply(x, kurtosis, na.rm = na.rm)
  else kurtosis(as.vector(x), na.rm = na.rm)
}

# Non-Parametric Part - Kernel ------------------------------------------------ #

# Gaussian Kernel 1 Dimension ------------------------------------------------- #

#' Gaussian Kernel
#'
#' This function computes a one-dimensional gaussian kernel.
#'
#' @param x a numeric vector.
#'
#' @family kernel functions
#'
#' @export
#'
NormalKernel1D <- function(x) {
  kernel1D <- (1/sqrt(2*pi)) * exp(-0.5 * x^2)
  kernel1D
}

# Quartic Kernel 1 Dimension -------------------------------------------------- #

#' Quartic Kernel
#'
#' This function computes a one-dimensional quartic kernel.
#'
#' @param x a numeric vector.
#'
#' @family kernel functions
#'
#' @export
#'
QuarticKernel1D <- function(x) {
  quartic1D <- ifelse(abs(x) <= 1, (15 / 16) * (1 - x^2)^2, 0)
  quartic1D
}

# One-Dimensional Kernel ------------------------------------------------------ #

#' One-Dimensional Kernel Estimation
#'
#' This function computes the one-dimensional kernel estimation needed for the
#' \code{\link{DSFM1D}} funtion to work. It is using the
#' \code{\link{QuarticKernel1D}} function.
#'
#' @param y a numeric matrix of the data.
#' @param I the total number of rows in \code{y}.
#' @param J the total number of columns in \code{y}.
#' @param x1 a numeric vector of the covariates.
#' @param u a numeric vector of the estimation grid.
#' @param U the length of \code{u}.
#' @param h a numeric vector of the bandwidth.
#'
#' @return It returns \eqn{\hat{p}_t(u)} and \eqn{\hat{q}_t(u)}.
#'
#' @family kernel functions
#'
#' @references Fengler, Matthias R, Wolfgang K Haerdle, and Enno Mammen (2007).
#' "A Semiparametric Factor Model for Implied Volatility Surface Dynamics".
#' In: \emph{Journal of Financial Econometrics 5.2}, pp. 189-218.
#'
#' @export
#'
KernelDensity1D <- function(y, I, J, x1, u, U, h) {

  pTHat <- matrix(0, I, U)
  qTHat <- pTHat
  for (t in 1:I) {
    pTjHat <- rep(0, J)
    pTjHatList <- replicate(U, list(pTjHat))
    qTjHat <- rep(0, J)
    qTjHatList <- replicate(U, list(qTjHat))
    for (n in 1:U) {
      pTjHatList[[n]] <- (1 / h[n]) * QuarticKernel1D((u[n] - x1) / h[n])
      qTjHatList[[n]] <- pTjHatList[[n]] * y[t, ]
    }
    pTHat[t, ] <- unlist(lapply(pTjHatList, function(x) (1 / J) * sum(x)))
    qTHat[t, ] <- unlist(lapply(qTjHatList, function(x) (1 / J) * sum(x)))
  }
  list(pHat = pTHat, qHat = qTHat)
}

# Two-Dimensional Kernel Product ---------------------------------------------- #

#' Two-Dimensional Kernel Estimation
#'
#' This function computes the two-dimensional kernel estimation needed for the
#' \code{\link{DSFM2D}} funtion to work. Following Fengler and al. (2007), it
#' is using the product of two \code{\link{QuarticKernel1D}} functions.
#'
#' @param y a list of numeric matrix of the data at each time \eqn{t}.
#' @param I the total number of matrix in \code{y}.
#' @param J a list of the total number of covariates in each matrix of \code{y}.
#' @param x1 a list of numeric vector for the first dimension at each time
#' \eqn{t}.
#' @param x2 a list of numeric vector for the second dimension at each time
#' \eqn{t}.
#' @param u a numeric matrix of the estimation grid.
#' @param U the length of \code{u}.
#' @param h a numeric matrix of the bandwidth.
#'
#' @return It returns \eqn{\hat{p}_t(u)}, \eqn{\hat{q}_t(u)} and the couples
#' \eqn{J_t\hat{p}_t(u)} and \eqn{J_t\hat{q}_t(u)}.
#'
#' @family kernel functions
#'
#' @references Fengler, Matthias R, Wolfgang K Haerdle, and Enno Mammen (2007).
#' "A Semiparametric Factor Model for Implied Volatility Surface Dynamics".
#' In: \emph{Journal of Financial Econometrics 5.2}, pp. 189-218.
#'
#' Borak, Szymon, Matthias R. Fengler, and Wolfgang K. Haerdle (2005)."DSFM
#' Fitting of Implied Volatility Surfaces". In: \emph{5th International
#' Conference on Intelligent Systems Design and Applications (ISDA'05)},
#' pp. 526-531. IEEE.
#'
#' @export
#'
KernelDensity2D <- function(y, I, J, x1, x2, u, U, h) {

  pTHat <- matrix(0, I, U)
  qTHat <- jQTHat <- jPTHat <- pTHat
  for (t in 1:I) {
    pTjHat <- matrix()
    pTjHatList <- replicate(U, list(pTjHat))
    qTjHat <- matrix()
    qTjHatList <- replicate(U, list(qTjHat))
    for (n in 1:U) {
      pTjHatList[[n]] <- ((1 / h[n,1]) * QuarticKernel1D((u[n,1] - x1[[t]]) /
                                                           h[n,1])) %*%
        ((1 / h[n,2]) * t(QuarticKernel1D((u[n,2] - x2[[t]]) / h[n,2])))

       qTjHatList[[n]] <- pTjHatList[[n]] * as.matrix(y[[t]])
    }
    # p.hat and q.hat function
    pTHat[t, ] <- unlist(lapply(pTjHatList, function(x) (1 / J[[t]]) * sum(x)))
    qTHat[t, ] <- unlist(lapply(qTjHatList, function(x) (1 / J[[t]]) * sum(x)))
    # Ji*p.hat and Ji*q.hat
    jPTHat[t, ] <- unlist(lapply(pTjHatList, function(x) sum(x)))
    jQTHat[t, ] <- unlist(lapply(qTjHatList, function(x) sum(x)))
  }
  list(pHat = pTHat, qHat = qTHat, jPHat = jPTHat, jQHat = jQTHat)
}

#' Fit a fused lasso over a lattice
#'
#' Calculates the fused lasso over a 2 dimensional square,
#' hexagonal lattice, or 3 dimensional cube. Allows for an
#' optional linear constraint on the result, seamlessly
#' handles missing data, and allows for an optional set of
#' observation weights.
#'
#' @param y
#'   the observed data points, as a matrix or 3 dimensional array.
#'   See Details for how this is interpreted for a hexagonal lattice.
#' @param weights
#'   optional sample weights. If provided, it must be the same length
#'   as \code{y}; it can be a vector, matrix, or array. Data will be
#'   assumed to be missing if either \code{y} or \code{weights} is
#'   \code{NA}. If missing, the weights will be assumed to be constant
#'   (unity) across all samples.
#' @param lambda
#'   a sequence of lambda values at which to produce a fit. Can be
#'   left blank (highly recommended for general use), at which point
#'   the algorithm will determine appropriate lambda values.
#' @param nlambda
#'   if \code{lambda} is missing, this determines the number of lambda
#'   values dynamically constructed by the algorithm.
#' @param lambdaMinRatio
#'   if \code{lambda} is missing, this determines the ratio between the
#'   largest and smallest \code{lambda} values. The values are evenly
#'   spaced on a log scale, so this ratio should typically be set fairly
#'   small.
#' @param E
#'   an optional matrix for supplying a set of linear constraints on
#'   the output. Columns corrispond to the (vectorised version of) the
#'   data samples, and each row to a given constraint. See Details for
#'   more information.
#' @param c
#'   an optional vector with one element per row of \code{E}. Gives
#'   the right hand side of the constraint. If missing but \code{E}
#'   is supplied, this will be assumed to be a vector of zeros.
#' @param latticeType
#'   the type of lattice to use. Will default to \code{square} if
#'   \code{y} is a matrix and \code{cube} if it is a 3 dimensional
#'   array.
#' @param rho
#'   a positive number used as a tuning parameter for the ADMM
#"   procedure.
#' @param eps
#'   stopping parameter for the iterative algorithm.
#' @param maxIter
#'   maximal number of iterations for the algorithm.
#' @param beta0
#'   initial starting point for the algorithm; if \code{NULL}, the
#'   default, the mean value of all inputs is used.
#' @param verbose
#'   logical. Should the function print out intermediate results
#'   as it is running.
#'
#' @details
#'
#'   When \code{latticeType} is set to \code{hex}, the algorithm
#'   will assume the data are distributed over a 2 dimensional
#'   hexagonal lattice. The coordinate system used assumes that
#'   even rows are offset by a factor of negative one half from
#'   the odd rows. Other patterns can be achieved by inserting
#'   \code{NA} values in over this grid and shifting data
#'   appropriately.
#'
#'   The \code{E} and \code{c} can be supplied to enforce the
#'   constraint Eb=c, where b is a vector of coefficients
#'   corrisponding to the vectorised version of \code{y}
#'   (i.e., \code{as.numeric(y)}). If \code{c} is missing, it will
#'   be assumed to be a vector of zeros.
#'
#' @references
#'
#'   Johnson, Nicholas A. "A Dynamic Programming Algorithm for the
#'    Fused Lasso and L0-Segmentation." Journal of Computational and
#'    Graphical Statistics 22.2 (2013): 246-260.
#'
#' @return
#'   an object of class 'glmgen', with methods for extracting
#'   coefficents and predictions.
#'
#' @examples
#'
#'  set.seed(1)
#'
#'  # Over a 20x10 grid
#'  y = matrix(rnorm(200),20,10)
#'  y[1:5,1:5] = y[1:5,1:5] + 2
#'  z = fusedLattice(y, lambda = 3, eps=0.0001, maxIter=100)
#'  print(matrix(z$beta, 20, 10))
#'  z = fusedLattice(y, lambda = 0.0001, eps=0.0001, maxIter=100)
#'  print(round(matrix(z$beta, 20, 10) - y),3)
#'
#'  # Linear constraint and missing values
#'  y = matrix(rnorm(8*12),8,12)
#'  y[1:5,1:5] = y[1:5,1:5] + 2
#'  y[2,10] = NA
#'  E = Matrix::Matrix(0,2,length(y))
#'  E[1,1] = 1
#'  E[1,2] = -1
#'  E[2,9] = 1
#'  E[2,11] = -1
#'  c = c(0,1)
#'
#'  z = fusedLattice(y, lambda = 0.00001, eps=0.0001, maxIter=100,
#'                   latticeType="square", E=E, c=c)
#'
#'  # 3D Lattice
#'  y = array(rnorm(3 * 4 * 5),dim=c(5,4,3))
#'  z = fusedLattice(y, lambda = 30, latticeType="cube")
#'  array(z$beta, dim=dim(y))
#'  z = fusedLattice(y, lambda = 0.0001, latticeType="cube")
#'  array(z$beta, dim=dim(y)) - y
#'
#' @author Taylor Arnold
#' @useDynLib glmgen lattice_R
#'
#' @importFrom  Matrix bandSparse t solve
#' @export
fusedLattice = function(y, weights, lambda=NULL, nlambda=20L,
                        lambdaMinRatio=1e-05, E=NULL, c=NULL,
                        latticeType=c("square","hex","cube"),
                        rho=1, eps=0.01, maxIter=20, beta0=NULL,
                        verbose=FALSE) {

  # Extract lattice type and problem dimensions
  if (is.matrix(y) && storage.mode(y) == "double") {
    latticeType = match.arg(latticeType)
  } else if (is.array(y) && length(dim(y)) == 3L &&
              storage.mode(y) == "double") {
    latticeType = "cube"
  } else stop("y must be a numeric matrix or 3D array.")
  latticeType = which(latticeType[[1]] == c("square","hex","cube")) - 1L

  # Create dummy weights if needed
  if (missing(weights))
    weights = rep(1, length(y))

  # Deal with missing values; if there are none we can use
  # faster methods in the C code, so flag the non-NA case
  naflag = TRUE
  index = which(is.na(y) | is.na(weights))
  if (length(index) == 0) {
    naflag = FALSE
  } else {
    y[index] = NA
    weights[index] = 0
  }

  # Deal with initialized value of beta0
  if (is.null(beta0)) {
    beta0 = rep(mean(y, na.rm=TRUE), length(y))
  } else {
    beta0 = as.double(beta0)
    if (length(beta0) != length(y))
      stop("beta0 must have the same length as y, if supplied")
  }

  # Extract the constraint Eb=c, if it exists
  if (missing(E)) {
    E = new("dgTMatrix", i = integer(0), j = integer(0),
            Dim = c(length(y), 0L), Dimnames = list(NULL, NULL),
            x = numeric(0), factors = list())
    E = new("dgTMatrix", i = integer(0), j = integer(0),
            Dim = c(length(y), length(y)), Dimnames = list(NULL, NULL),
            x = numeric(0), factors = list())
  } else {
    if (ncol(E) != length(y))
      stop("E must have a column for each element in y.")
    E = as(E,"dgTMatrix")
  }
  if (missing(c)) c = rep(0, nrow(E))
  if (nrow(E) != length(c))
    stop("Number of rows of E must match the length of c!")

  # Calculate an estimate of the lambda values:
  if (is.null(lambda)) {
    if (latticeType == 0L) {
      D1 = Matrix::bandSparse(nrow(y), m = nrow(y), c(0, 1),
             diagonals = list(rep(-1, nrow(y)), rep(1, nrow(y) - 1)))
      D2 = Matrix::bandSparse(ncol(y), m = ncol(y), c(0, 1),
             diagonals = list(rep(-1, ncol(y)), rep(1, ncol(y) - 1)))
      max_lam = max(max(abs(Matrix::solve(Matrix::t(D1)) %*% y),na.rm=TRUE),
                    max(abs(Matrix::solve(Matrix::t(D1)) %*% y),na.rm=TRUE))
    }

    min_lam = max_lam * lambdaMinRatio
    lambda = exp(seq(log(max_lam), log(min_lam), length.out=nlambda))
  }

  output = list(y=y, w=weights, lambda=lambda, dims=dim(y),
                latticeType=latticeType,
                beta=matrix(NA, ncol=length(lambda), nrow=length(y)))

  for (i in 1:length(lambda)) {
    output$beta[,i] = .Call("lattice_R",
                            sY = y,
                            sW = weights,
                            sLambda = as.double(lambda[i]),
                            sRho = as.double(rho),
                            sEps = as.double(eps),
                            sMaxiter = as.integer(maxIter),
                            sVerbose = as.integer(verbose),
                            sNaflag = as.integer(naflag),
                            sLatticeType = as.integer(latticeType),
                            sE = E,
                            sC = as.double(c),
                            sBeta0 = as.double(beta0),
                            PACKAGE = "glmgen")
  }

  index = which(is.na(y) | is.na(weights))
  output$beta[index,] = NA
  class(output) = c("fusedLattice", "glmgen")
  output
}


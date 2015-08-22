#' Fit a fused lasso over a Decomposed Graph
#'
#' Calculates the fused lasso over a graph, which has been
#' decomposed into a series of chains.
#'
#' @param y
#'   the observed data points as a numeric vector.
#' @param weights
#'   optional sample weights. If provided, it must be the same length
#'   as \code{y}.
#' @param edges
#'    a list object describing the edges in the chains. Each element in
#'    the list should be a numeric vector describing the chains. See
#'    Details for a description of how this should be specified.
#' @param edgeWeights
#'   an optional list of edge weights. Should be a list of the same
#'   length as \code{edges}, with each component having the same length
#'   as well.
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
#' @param rho
#'   a positive number used as a tuning parameter for the ADMM
#'   procedure.
#' @param eps
#'   stopping parameter for the iterative algorithm.
#' @param maxIter
#'   maximal number of iterations for the algorithm.
#' @param beta0
#'   initial starting point for the algorithm; if \code{NULL}, the
#'   default, the mean value of all inputs is used.
#' @param method
#'   method used to solve the 1D chains in the ADMM algorithms.
#'   Currently either \code{dp} for Johnson's dynamic program or
#'   \code{prox} for the proximal optimization method for Barbero and
#'   Suvrit.
#' @param verbose
#'   logical. Should the function print out intermediate results
#'   as it is running.
#'
#' @details
#'
#'   Each element in the list \code{edges} should idealy contain
#'   some permutation of the integers 1 through the length of \code{y}
#'   describing a path over the graph. This path may have holes, which
#'   can be specified by inserting negative, zero, or \code{NA} values
#'   into the list element. It is not strictly required that each chain
#'   visit every node in the graph (these orphaned nodes are simply
#'   considered to be have no constraints), however it is strictly not
#'   allowed for a cycle to visit a node multiple times.
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
#'   Barbero, Alvaro, and Suvrit Sra. "Modular proximal optimization
#'    for multidimensional total-variation regularization." arXiv
#'    preprint arXiv:1411.0589 (2014).
#'
#'   Barbero, Alvaro, and Suvrit Sra. "Fast Newton-type methods for
#'    total variation regularization." Proceedings of the 28th
#'    International Conference on Machine Learning (ICML-11). 2011.
#'
#' @return
#'   an object of class 'glmgen', with methods for extracting
#'   coefficents and predictions.
#'
#' @examples
#'
#'  set.seed(1)
#'
#'
#' @author Taylor Arnold
#' @useDynLib glmgen graph_fused_R
#'
#' @importFrom  Matrix bandSparse t solve
#' @export
fusedGraph = function(y, weights, edges, edgeWeights=NULL,
                        lambda=NULL, nlambda=20L, lambdaMinRatio=1e-05,
                        E=NULL, c=NULL,
                        rho=1, eps=0.01, maxIter=20, beta0=NULL,
                        method=c("prox","dp"),
                        verbose=FALSE) {

  # Create dummy observation weights if needed
  if (missing(weights))
    weights = rep(1, length(y))

  # Test edges and convert to C format (0 indexed; -1 for holes)
  if (!is.list(edges)) stop("edges must be a list object.")
  edges = lapply(edges, function(v) as.integer(v - 1))
  edges = lapply(edges, function(v) {v[is.na(v) | v < 0] = -1; v})
  okay = sapply(edges, function(v) !any(duplicated(v[v >= 0])) &
              !any(is.na(match(v, (-1):(length(y)-1)))))
  if (any(!okay))
    stop("bad element in edges list: ",
         paste(which(!okay),collapse=","))

  # Select the method for solving the sub-problems
  method = match.arg(method)
  if (!is.null(edgeWeights)) {
    if (method != "prox")
    warning("Edge weights currently are not suppored in the",
            "DP algorithm;\nusing proximal method ('prox') instead")
    method = "proxW"
  }
  methodType = which(method == c("dp", "prox", "proxW")) - 1L

  # Test/Create edge weights
  if (is.null(edgeWeights)) {
    edgeWeights = lapply(edges, function(v) rep(1,length(v)))
  } else {
    if (length(edgeWeights) != length(edges) ||
        any(sapply(edges,length) != sapply(edgeWeights,length)))
      stop("length of edgeWeights and its elements must",
           "match that of edges.")
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
    # if (latticeType == 0L) {
    #   D1 = Matrix::bandSparse(nrow(y), m = nrow(y), c(0, 1),
    #          diagonals = list(rep(-1, nrow(y)), rep(1, nrow(y) - 1)))
    #   D2 = Matrix::bandSparse(ncol(y), m = ncol(y), c(0, 1),
    #          diagonals = list(rep(-1, ncol(y)), rep(1, ncol(y) - 1)))
    #   max_lam = max(max(abs(Matrix::solve(Matrix::t(D1)) %*% y),na.rm=TRUE),
    #                 max(abs(Matrix::solve(Matrix::t(D1)) %*% y),na.rm=TRUE))
    # }

    # min_lam = max_lam * lambdaMinRatio
    # lambda = exp(seq(log(max_lam), log(min_lam), length.out=nlambda))
  }

  output = list(y=y, w=weights, lambda=lambda, dims=dim(y),
                beta=matrix(NA, ncol=length(lambda), nrow=length(y)))

  for (i in 1:length(lambda)) {
    output$beta[,i] = .Call("graph_fused_R",
                            sY = as.numeric(y),
                            sW = as.numeric(weights),
                            sEdge = as.integer(unlist(edges)),
                            sWedge = as.numeric(unlist(edgeWeights)),
                            sEdgeLen = as.integer(lapply(edges,length)),
                            sLambda = as.double(lambda[i]),
                            sRho = as.double(rho),
                            sEps = as.double(eps),
                            sMaxiter = as.integer(maxIter),
                            sVerbose = as.integer(verbose),
                            sMethodType = as.integer(methodType),
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


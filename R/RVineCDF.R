#' CDF of an R-Vine Copula Model
#'
#' This function calculates the cumulative distribution function of a
#' d-dimensional R-vine copula using a naïve Monte-Carlo simulation
#'
#' @param data An N x d data matrix that specifies where the cdf shall
#' be evaluated.
#' @param RVM An \code{\link{RVineMatrix}} object including the structure and
#' the pair-copula families and parameters.
#' @param n integer; default \code{n = 1e4}, the number of Monte-Carlo simulations.
#' @param check.pars logical; default is \code{TRUE}; if \code{FALSE}, checks
#' for family/parameter-consistency are ommited (should only be used with
#' care).
#' @author Thibault Vatter
#'
#' @seealso \code{\link{BiCopHfunc}}, \code{\link{RVineMatrix}},
#' \code{\link{RVineMLE}}, \code{\link{RVineAIC}}, \code{\link{RVineBIC}},
#' \code{\link{RVinePDF}}
#'
#' @examples
#'
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#'
#' # define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#'
#' # define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#'
#' # define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#'
#' # define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#'
#' # compute the cdf at (0.1, 0.2, 0.3, 0.4, 0.5)
#' RVineCDF(c(0.1, 0.2, 0.3, 0.4, 0.5), RVM)
#'
RVineCDF <- function(data, RVM, n = 1e4, check.pars = TRUE) {
    args <- preproc(c(as.list(environment()), call = match.call()),
                    check_data,
                    fix_nas,
                    check_if_01,
                    check_RVMs,
                    prep_RVMs)
    if (is.matrix(data)) {
        res <- apply(data, 1, function(x) RVineCDF(x, RVM, n, FALSE))
    } else {
        u <- t(RVineSim(n, RVM)) # it's obviously inefficient to re-simulate
        # for each observation, furthermore, it could be improved with ghalton (qrng)
        return(mean(apply(u <= data, 2, all)))
    }
    return(res)
}

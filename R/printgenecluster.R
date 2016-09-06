#' Print groups of clustered genes.
#'
#' @param x An object of class \code{genecluster}.
#' @param ... Additional arguments to be passed (not in use).
#'
#' @examples
#' unbundleCluster(clustering(exanova, clust.metric = 'manhattan', clust.method = 'average'))
#'
#' @export

print.genecluster <- function (x, ...)
{
  y <- x$groups
  y
}

#' Cluster all interesting (differentially expressed)
#'  genes using hierarchical clustering (based on all groups together).
#' @noRd
#' @param resultAOV A list of data frames of gene expressions
#'  (genes in rows, samples in columns).
#' @param clust.metric The method to calculate a distance measure
#'  used in hierarchical clustering, possible names: \code{"euclidean"},
#' \code{"maximum"}, \code{"manhattan"}, \code{"canberra"},
#'  \code{"binary"} or \code{"minkowski"}.
#' @param clust.method The agglomeration method of clustering genes.
#'  This should be one of: \code{"ward.D"},
#'  \code{"ward.D2"}, \code{"single"},
#'  \code{"complete"}, \code{"average" (= UPGMA)},
#'  \code{"mcquitty" (= WPGMA)}, \code{"median" (= WPGMC)}
#'  or \code{"centroid" (= UPGMC)}.
#' @param dist.matrix A matrix with calculated distances
#' to be used as a metric by
#' \code{hclust} function.
#' The default distance is \code{euclidean}.
#'
#' @return An object of class \code{hclust}.
#'
#' @examples
#' clustering(exanova, clust.metric = 'manhattan', clust.method = 'average')

clustering <- function(resultAOV, clust.metric = NULL,
                       clust.method = NULL, dist.matrix = NULL) {
    clust <- do.call(cbind, resultAOV)
    row.names(clust) <- row.names(resultAOV[[1]])

    if (is.null(dist.matrix)) {
        if (is.character(clust.metric) & (length(clust.metric) == 1)) {
            distance <- dist(clust, method = clust.metric)
        } else {
            distance <- dist(clust, method = "euclidean")
        }
    } else {
        distance <- dist.matrix
    }
    if (is.null(clust.method)) {
        clust.method <- "average"
    }
    cluster1 <- hclust(distance, method = clust.method)
    cluster1
}

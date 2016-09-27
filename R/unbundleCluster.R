#' Find all groups - rearrange the dendrogram:
#'  find which elements were merged in each step.
#' @noRd
#' @param hc An object of class hclust.
#'
#' @return An object of class \code{genecluster},
#'  a list of labels of elements that were merged in each step.
#'
#' @examples
#' clusteredGenes <- clustering(exanova,
#'  clust.metric = 'manhattan', clust.method = 'average')
#' unbundleCluster(clusteredGenes)


unbundleCluster <- function(hc) {
    roots <- hc$labels[hc$order]
    merged <- hc$merge
    leaves <- lapply(1:(length(hc$labels) - 1), function(i) {
        hc$labels[-ordering(i, merged)]
    })
    out <- list(groups = c(roots, leaves), hclust = hc)
    names(out[[1]]) <- paste0("G", seq_along(out[[1]]))
    class(out) <- append("genecluster", class(out))
    out
}

ordering <- function(x, merged) {
    if (merged[x, 1] < 0) {
        a <- merged[x, 1]
    } else {
        a <- ordering(merged[x, 1], merged)
    }
    if (merged[x, 2] < 0) {
        b <- merged[x, 2]
    } else {
        b <- ordering(merged[x, 2], merged)
    }
    c(a, b)
}


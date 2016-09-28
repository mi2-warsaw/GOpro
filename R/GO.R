#' Constructor function for the class \code{GO}.
#' @noRd
#' @param allGO A list of lists of vectors each with a p-value
#'  obtained from Fisher test, GO id and all members of each group.
#' @param topGO A list of lists with gene ontology terms
#'  identifiers characteristic for each group.
#'
#' @return An object of class \code{GO}.
#'
#' @examples
#' clustered <- clustering(exanova, clust.metric = 'manhattan',
#'  clust.method = 'average')
#' genecl <- unbundleCluster(clustered)
#' # find all and then top GO terms
#' # geneUniverse is a universe of all genes of interest
#' geneUniverse <- rownames(exrtcga[[1]])
#' all.gos <- findAllGOs(geneUniverse, clustered,
#'  min = 10, max = 200)
#' top.gos <- findTopGOs(all.gos, sig.level = 0.05, top = 2)
#' # create an object of class \code{GO}
#' GO(all.gos, top.gos)
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom S4Vectors DataFrame List


GO <- function(allGO, topGO) {
    out <- DataFrame(allGO = List(lapply(allGO, CharacterList)), topGO = NumericList(topGO))
    invisible(out)
}

#' Print grouped genes with the information of
#'  the labels and GO term's IDs.
#' @noRd
#' @param x An object of class \code{GO}.
#' @param ... Additional arguments to be passed (not in use).
#'
#' @return A data frame with labels, gene's aliases
#'  and GO terms for each group.
#'
#' @examples
#' clustered <- clustering(exanova)
#' genecl <- unbundleCluster(clustered)
#' # find all and then top GO terms
#' # geneUniverse is a universe of all genes of interest
#' geneUniverse <- rownames(exrtcga[[1]])
#' all.gos <- findAllGOs(geneUniverse, clustered,
#'  max = 300, onto = c('MF', 'BP'))
#' top.gos <- findTopGOs(all.gos, sig.level = 0.05, top = 3)
#' # create an object of class \code{GO}
#' printout <- GO(all.gos, top.gos)
#'

printGO <- function(x, ...) {
    allGO <- x[[1]]
    topGO <- x[[2]]
    y <- cbind(profile = names(topGO), GOs = lapply(topGO, names),
               p.values = lapply(topGO, function(x) round(x, 3)),
               GENES = sapply(allGO,
                              function(x) {
                                  if (length(x) == 0) {
                                      NA
                                  } else {
                                      paste0(x[[1]][3:length(x[[1]])],
                                             collapse = " ")
                                  }
                              }))
    y
}

#' Create output.
#' @noRd
#' @param GO An object returned by \code{GO} function.
#'
#' @return An S4 object.
#'

output <- function(GO){
    printedGO <- printGO(GO)
    DataFrame(printedGO)
}

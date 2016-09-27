#' Calculate frequencies of each ontology in each group (CC, BP, MF).
#' @noRd
#' @param extendedGO A data frame with a one column of group
#'  and another column with GO term (returned by \code{extendGO}).
#'
#' @return A data frame with the most frequently occurring
#'  ontologies in each group.
#'
#' @examples
#' tukey.results <- groupByTukey(exanova, sig.level = 0.1)
#' # find all and then top GO terms
#' # geneUniverse is a universe of all genes of interest,
#' # i.e. all genes measured by microarray
#' geneUniverse <- rownames(exrtcga[[1]])
#' all.gos <- findAllGOs(geneUniverse, tukey.results)
#' top.gos <- findTopGOs(all.gos, sig.level = 0.1, top = 3)
#' extended <- extendGO(top.gos)
#' ontologyOver(extended)
#'

ontologyOver <- function(extendedGO) {
    group <- extendedGO$GROUP
    ontology <- extendedGO$ONTOLOGY
    BP <- vector("numeric", length(unique(group)))
    CC <- vector("numeric", length(unique(group)))
    MF <- vector("numeric", length(unique(group)))
    no <- length(unique(group))
    BP <- sapply(seq_len(no), function(x) {
        sum(extendedGO[extendedGO$GROUP == x, "ONTOLOGY"] == "BP")
    })
    CC <- sapply(seq_len(no), function(x) {
        sum(extendedGO[extendedGO$GROUP == x, "ONTOLOGY"] == "CC")
    })
    MF <- sapply(seq_len(no), function(x) {
        sum(extendedGO[extendedGO$GROUP == x, "ONTOLOGY"] == "MF")
    })
    counted.onto <- cbind.data.frame(group = seq_len(no),
                                     onto = rep(c("BP", "CC", "MF"),
                                                each = no),
                                     counts = c(BP, CC, MF))
    counted.onto$counts[is.na(counted.onto$counts)] <- 0
    sapply(1:no, function(x) {
        three <- counted.onto[c(x, x + no, x + 2 * no), 3]
        if (all(three == 0)) {
            return("No function")
        }
        names(three) <- counted.onto[c(x, x + no, x + 2 * no), 2]
        names(which(three == max(three)))
    })
}

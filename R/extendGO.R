#' Find the description of characteristic GO functions.
#' @noRd
#' @param TopGOs A list of vectors with significant GO IDs.
#'
#' @return A data frame with an extended description
#'  (term, definition, and ontology) of gene ontology terms.
#'
#' @examples
#' tukey.results <- groupByTukey(exanova, sig.level = 0.1)
#' # find all and then top GO terms
#' # geneUniverse is a universe of all genes
#' # of interest, i.e. all genes measured by microarray
#' geneUniverse <- rownames(exrtcga[[1]])
#' all.gos <- findAllGOs(geneUniverse, tukey.results, min = 8, max = 400)
#' top.gos <- findTopGOs(all.gos, sig.level = 0.1, top = 2)
#' extendGO(top.gos)
#'
#'
#' @import GO.db
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges IntegerList CharacterList

extendGO <- function(TopGOs) {
    group.size <- sapply(TopGOs, length)
    nulls <- which(unlist(lapply(TopGOs, is.na)) == TRUE)
    TopG <- vector("character", sum(group.size))
    TopG[nulls] <- NA
    TopG[-nulls] <- na.omit(unlist(lapply(TopGOs, names)))
    GO.descr <- lapply(TopG, function(GOID) {
        if (!is.na(GOID)) {
                                select(GO.db, keys = GOID,
                                  columns = c("TERM", "DEFINITION",
                                              "ONTOLOGY"),
                                  keytype = "GOID")
        } else {
            data.frame(cbind(GOID = NA, TERM = NA,
                             DEFINITION = NA, ONTOLOGY = NA))
        }
    })
    gnames <- seq_along(group.size)
    assignedFunctions <- data.frame(rep(gnames, group.size),
                                    matrix(unlist(GO.descr),
                                           ncol = 4, byrow = TRUE))
    names(assignedFunctions) <- c("GROUP", "GO ID",
                                  "TERM", "DEFINITION", "ONTOLOGY")
    assignedFunctions
}

extendoutput <- function(top.gos){
    ext <- extendGO(top.gos)
    DataFrame(GROUP = IntegerList(ext[, 1]),
                  GOID = CharacterList(ext[, 2]),
                  TERM = CharacterList(ext[, 3]),
                  DEFINITION = CharacterList(ext[, 4]),
                  ONTOLOGY = CharacterList(ext[, 5]))
}

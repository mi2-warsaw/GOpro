#' Find top gene ontology terms for groups of genes
#'  based on the result obtained from \code{findAllGOs}.
#' @noRd
#' @param AllGOs A list of lists of vectors each with a p-value
#'  obtained from Fisher test, GO id and all members of each group.
#' @param sig.level A numeric value, a significance
#'  level used in BH correction for multiple testing.
#' @param top A number of the most characteristic
#'  functions of groups of genes to be returned.
#'
#' @return A list of lists of corrected p-values of Fisher test of
#'  characteristic for each group GO terms named with
#'  gene ontology terms identifiers.
#'
#' @examples
#' tukey.results <- groupByTukey(exanova, parallel = FALSE, 0.1)
#' #geneUniverse is a universe of all genes of interest
#' geneUniverse <- rownames(exrtcga[[1]])
#' all.gos <- findAllGOs(geneUniverse, tukey.results,
#'  min = 5, max = 500)
#' findTopGOs(all.gos, sig.level = 0.05, top = 10)
#'

findTopGOs <- function(AllGOs, sig.level, top) {
    lapply(AllGOs, function(group.no) {
        if (length(group.no) == 0) {
            return(NA)
        }
        p.val <- sapply(group.no, function(p.no) {
            record <- p.no
            if (length(record) == 0) {
                return(NA)
            }
            if (is.na(record[1])) {
                return(NA)
            }
            p.val <- record[1]
            if (length(p.val) == length(record[1])) {
                names(p.val) <- record[2]
            } else {
                if (length(p.val) == 0)
                    return(NA) else
                        names(p.val) <- record[2, which(!is.na(record[1]))]
            }
            p.val
        })

        p.val <- na.omit(p.val)
        p.adj <- p.adjust(p.val, "fdr")
        if (length(p.adj) == 0) {
            return(NA)
        }
        if (min(p.adj) > sig.level) {
            return(NA)
        }

        noOfGOs <- top
        if (sum(p.adj < sig.level) < top) {
            noOfGOs <- sum(p.adj < sig.level)
        }
        sorted <- sort(p.adj[p.adj < sig.level])

        is.equal <- function(noOfGOs) {
            if (noOfGOs == sum(p.adj < sig.level)) {
                return(noOfGOs)
            }
            if (sorted[noOfGOs] == sorted[noOfGOs + 1]) {
                noOfGOs <- noOfGOs + 1
                noOfGOs <- is.equal(noOfGOs)
            }
            noOfGOs
        }
        most.char <- sorted[1:is.equal(noOfGOs)]
        most.char
    })
}

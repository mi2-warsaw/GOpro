#' Find all ontology terms annotated with genes.
#' @noRd
#' @param geneUni A vector of names of all genes that
#'  were primarily included in the analysis,
#'  i.e. all genes measured by the microarray.
#' @param geneSel A list of groups of selected genes
#'  or an hclust object representing clustered genes.
#' @param min A minimum number of genes that a function needs
#'  to be represented by to be considered in the analysis.
#' @param max A maximum number of genes that a function needs
#'  to be represented by to be considered in the analysis.
#' @param parallel A logical value indicating if a task
#'  should be run on more than one core.
#'  The default is \code{FALSE}.
#' @param onto An ontology or ontologies to be searched for
#'  significant GO terms, at least one of 'MF' (molecular function),
#'  'BP' (biological process), 'CC' (cellular component).
#'
#' @return A list of lists with all gene ontology terms
#'  identifiers relevant to each group of genes and p-values for Fisher's test.
#'
#' @examples
#' tukey.results <- groupByTukey(exanova, 0.1)
#' #geneUniverse is a universe of all genes of interest
#' geneUniverse <- rownames(exrtcga[[1]])
#' findAllGOs(geneUniverse, tukey.results)
#'
#'
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import doParallel
#' @import foreach
#' @import parallel

findAllGOs <- function(geneUni, geneSel, onto = c("MF", "BP", "CC"),
                       min = 5, max = 500, parallel = FALSE) {
    stopifnot(is.list(geneSel) &
                  is.vector(geneUni) &
                  all(is.character(geneUni)))

    if (!all(onto %in% c("MF", "BP", "CC")))
        stop(paste(onto, "is not a valid ontology domain."))

    if (any(is.null(geneSel$merge) &
            is.null(geneSel$dist.method) &
            is.null(geneSel$order))) {
        cgroups <- geneSel
    } else {
        cgroups <- unbundleCluster(geneSel)[[1]]
    }
    iter <- length(cgroups)
    all.annot <- select(org.Hs.eg.db, keys = geneUni,
                        keytype = "SYMBOL", columns = "GO")
    all.annot <- all.annot[all.annot$ONTOLOGY %in% onto, ]

    if (nrow(all.annot) == 0)
        stop("No Gene Ontology terms in the given ontology domain.")
    go2gene.lst.p <- aggregate(SYMBOL ~ GO, all.annot, invisible)
    counted.genes <- unname(sapply(go2gene.lst.p$SYMBOL, length))
    names(counted.genes) <- go2gene.lst.p$GO

    go2gene.lst <- go2gene.lst.p[which(
        counted.genes >= min & counted.genes <= max), ]
    uniSize <- length(geneUni)

    allFisherTests <- function(i) {
        in.group <- cgroups[[i]]
        group.gos <- unique(na.omit(
            all.annot$GO[(all.annot$SYMBOL %in% in.group)]))
        len.in <- length(in.group)
        if (length(group.gos) == 0) {
            return(NA)
        }
        lapply(group.gos, function(group.go) {
            if (sum(go2gene.lst$GO == group.go) == 0) {
                return(NA)
            } else {
                this.genes <- unname(unlist(
                    # names of genes with this f
                    go2gene.lst[go2gene.lst[, 1] == group.go, 2]))
                # the number of all genes with this f
                withGO <- unname(counted.genes[group.go])

                pp <- sum(sapply(seq_along(in.group),
                                 function(x) in.group[x] == this.genes))
                #liczba genow w grupie o danej f
                pn <- withGO - pp
                #liczba pozostalych genow o danej funkcji
                np <- len.in - pp
                #liczba genow w grupie bez danej funkcji
                nn <- uniSize - pp - pn - np
                #liczba genow poza grupa bez danej f

                if (!is.na(sum(c(pp, pn, np, nn))) &
                    all(c(pp, pn, np, nn) >= 0)) {
                    result <- fisher_test(pp, pn, np, nn)
                } else {
                    return(NA)
                }
            }
            out_rec <- c(result, group.go, in.group)
            out_rec
        })
    }
    if (missing(parallel)) {
        parallel <- FALSE
    }
    z <- 0
    if (parallel == TRUE) {
        cores <- detectCores() - 1
        cl <- makeCluster(cores)
        registerDoParallel(cl)
        out <- foreach(z = seq_len(iter), .packages = "GOpro") %dopar%
            allFisherTests(z)
        names(out) <- names(cgroups)
        stopCluster(cl)
        out
    } else {
        out <- foreach(z = seq_len(iter)) %do%
            allFisherTests(z)
        names(out) <- names(cgroups)
        out
    }
    remove.NA <- lapply(out, function(x) which(is.na(x)))
    out2 <- lapply(names(remove.NA), function(x) {
        out[[x]][-remove.NA[[x]]]
    })
    names(out2) <- names(cgroups)
    out2
}

#' Find top Gene Ontology terms for given genes.
#'
#' @param groups a MultiAssayExperiment object
#' containing an ExperimentList class object representing
#' gene expressions for at least 3 cohorts.
#'  Rows must be named with genes' aliases.
#'  The order of samples and genes has to be the same for each
#'  ExperimentList class object.
#' @param topAOV A numeric value, a number of most
#'  significantly differentiated genes to be returned.
#' @param sig.levelAOV a numeric value, a significance level
#'  used in BH correction for multiple testing (\code{aovTopTest}).
#' @param sig.levelGO A numeric value, a significance level used
#'  in BH correction for multiple testing (\code{findTopGOs}).
#' @param minGO A minimum number of functions that a gene
#'  needs to represent to be considered as frequent.
#' @param maxGO A maximum number of functions that a gene
#'  needs to represent to be considered as frequent.
#' @param parallel A logical value indicating if a task
#'  should be run on more than one core.
#' @param clust.metric The method to calculate a distance measure
#'  used in hierarchical clustering, possible names: "euclidean",
#'  "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param clust.method The agglomeration method used to cluster genes.
#'  This should be #'one of "ward.D",
#'  "ward.D2", "single", "complete", "average" (= UPGMA),
#'  "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param topGO A number of the most characteristic functions
#'  of groups of genes to be returned.
#' @param grouped A method of grouping genes, one of 'tukey' and 'clustering'.
#' @param dist.matrix A matrix with calculated distances to be
#'  used as a metric by \code{hclust} function.
#' @param sig.levelTUK A numeric value, a significance level used in
#'  Tukey's all pairwise comparison (\code{groupByTukey}).
#' @param onto An ontology or ontologies to be searched
#'  for significant GO terms,
#'  at least one of 'MF' (molecular function),
#'  'BP' (biological process), and 'CC' (cellular component).
#' @param extend A logical value indicating if an extended
#'  version of the output should be presented.
#' @param over.rep A logical value indicating if an over
#'  represented GO terms should be presented in the plot.
#'
#' @return A data frame containing the top gene ontology terms for
#'  each group of genes and the gene aliases.
#'
#' @examples
#' findGO(exrtcga, grouped = 'clustering', topGO = 10, onto = 'MF')
#' findGO(exrtcga, grouped = 'tukey', topGO = 2, extend = TRUE)
#' @export
#'
#' @useDynLib GOpro, .registration=TRUE
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom MultiAssayExperiment assay experiments

findGO <- function(groups, topAOV = 50, sig.levelAOV = 0.05,
                   parallel = FALSE, grouped = "tukey",
                   sig.levelGO = 0.05, minGO = 5, maxGO = 500,
                   clust.metric = NULL, clust.method = NULL,
                   dist.matrix = NULL, topGO = 3, sig.levelTUK = 0.05,
                   onto = c("MF", "BP", "CC"), extend = FALSE, over.rep = FALSE) {
    groups <- assay(experiments(groups))
    aov.results <- aovTopTest(groups, topAOV, sig.levelAOV, parallel)
    geneUniverse <- rownames(groups[[1]])
    if (!(grouped == "tukey" | grouped == "clustering")) {
        stop(paste("There is no grouping method called", grouped))
    }
    if (grouped == "tukey") {
        tukey.results <- groupByTukey(aov.results, parallel, sig.levelTUK)
        AllGOs <- findAllGOs(geneUniverse, tukey.results,
                             onto, minGO, maxGO, parallel)
        TopGOs <- findTopGOs(AllGOs, sig.levelGO, topGO)
        if (extend == TRUE) {
            out <- extendoutput(TopGOs)
        } else {
            out <- output(GO(AllGOs, TopGOs))
        }
    }
    if (grouped == "clustering") {
        cluster.results <- clustering(aov.results,
                                      clust.metric, clust.method, dist.matrix)
        all.gos <- findAllGOs(cluster.results,
                              geneUni = rownames(groups[[1]]))
        top.gos <- findTopGOs(all.gos, sig.levelGO, topGO)
        unbundled <- unbundleCluster(cluster.results)
        printout <- printGO(GO(all.gos, top.gos))
        plotg(printout, unbundled, top.gos,
              over.represented = over.rep)
        if (extend == TRUE) {
            out <- extendoutput(top.gos)
        } else {
            out <- output(GO(all.gos, top.gos))
        }
    }
    out
}

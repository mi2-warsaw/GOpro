#' Expressions of human genes.
#'
#' A dataset containing gene expressions of 300 human genes
#' randomly chosen from a list returned by function
#' \code{prepareData(RTCGA = TRUE, cohorts = c('leukemia', 'colon', 'bladder'))}.
#'
#' @format A MultiAssayExperiment object of 3 listed
#' experiments with user-defined names and respective classes.
#' Containing an ExperimentList class object of length 3:
#' [1] leukemia: matrix with 300 rows and 173 columns
#' [2] colon: matrix with 300 rows and 190 columns
#' [3] bladder: matrix with 300 rows and 122 columns.
#'
#' @return data

"exrtcga"

# rtcga <- prepareData(RTCGA = TRUE, cohorts = c('leukemia','colon','bladder'))
# set.seed(5)
# n <- sort(unique(ceiling(runif(306, 0, 8546))))
# exrtcga <- lapply(rtcga, function(x) x[n, ])
# exanova <- aovTopTest(minirtcga, top = 50, sig.level = 0.1,
#  parallel = FALSE)
# extukey <- groupByTukey(exanova, parallel = FALSE, 0.1)
# geneUniverse <- rownames(minirtcga[[1]])
# exallgoanova <- findAllGOs(geneUniverse, extukey, min = 10,
# max = 200, parallel = FALSE, onto = c('MF', 'BP', 'CC')),
# cluster.results <- clustering(exanova,
#   clust.metric = "euclidean", clust.method = "centroid")
# exallgoclust <- findAllGOs(geneUniverse, cluster.results, min = 10,
#   max = 200, parallel = FALSE, onto = c('MF', 'BP', 'CC'))

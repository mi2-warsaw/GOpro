#' Expressions of human genes.
#'
#' A dataset containing gene expressions of 300 human genes
#' randomly chosen from a list returned by function
#' \code{prepareData(RTCGA = TRUE, cohorts = c('leukemia', 'colon', 'bladder'))}.
#'
#' @format A list of 3 data frames with 300 rows each and 173, 190, and 122 variables respectively. Each variable corresponds to one probe.

"exrtcga"

#' Expressions of differently expressed human genes.
#'
#' A dataset containing gene expressions of 50 human genes returned by function
#' \code{aovTopTest(exrtcga, top = 50, sig.level = 0.1)}.
#'
#' @format A list of 3 data frames with 50 rows each and 173, 190, and 122 variables respectively. Each variable corresponds to one probe.

"exanova"

# rtcga <- prepareData(RTCGA = TRUE, cohorts = c('leukemia', 'colon', 'bladder'))
# set.seed(5)
# n <- sort(unique(ceiling(runif(306, 0, 8546))))
# exrtcga <- lapply(rtcga, function(x) x[n, ])
# exanova <- aovTopTest(minirtcga, top = 50, sig.level = 0.1, parallel = FALSE)
# extukey <- groupByTukey(exanova, parallel = FALSE, 0.1)
# geneUniverse <- rownames(minirtcga[[1]])
# exallgoanova <- findAllGOs(geneUniverse, extukey, min = 10, max = 200, parallel = FALSE, onto = c('MF', 'BP', 'CC')),
# cluster.results <- clustering(exanova, clust.metric = "euclidean", clust.method = "centroid"),
# exallgoclust <- findAllGOs(geneUniverse, cluster.results, min = 10, max = 200, parallel = FALSE, onto = c('MF', 'BP', 'CC'))

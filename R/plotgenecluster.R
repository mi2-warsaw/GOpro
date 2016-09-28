#' Plot gene dendrogram with labels.
#' @noRd
#' @param x an object of class \code{geneclusterplot}.
#' @param over.represented a logical value, if \code{TRUE} then
#'  over represented GO categories are coded with colors for each branch.
#'
#' @return n
#'
#' @examples
#' clustered <- clustering(exanova,
#'  clust.metric = 'manhattan', clust.method = 'average')
#' genecl <- unbundleCluster(clustered)
#' # find all and then top GO terms
#' # geneUniverse is a universe of all genes of interest
#' geneUniverse <- rownames(exrtcga[[1]])
#' all.gos <- findAllGOs(geneUniverse, clustered)
#' top.gos <- findTopGOs(all.gos, sig.level = 0.05, top = 3)
#' # create an object of class \code{GO}
#' printout <- GO(all.gos, top.gos)
#' plotg(geneclusterplot(printout, genecl, top.gos),
#'  over.represented = TRUE)
#'
#'
#' @import dendextend

plotg <- function(printout, genecl, TopGOs, over.represented) {
    hc <- genecl$hclust
    dendro <- as.dendrogram(hc)
    nodes <- get_nodes_xy(dendro)
    roots <- which(nodes[, 2] == 0)
    new.order <- printout[, 1]
    new.order[roots] <- printout[1:length(roots),
                                 1]
    sorted <- sort(nodes[which(nodes[, 2] > 0),
                         2])
    leaves <- sapply(seq_len(length(roots) - 1),
                     function(x) which(nodes[, 2] == sorted[x]))
    new.order[leaves] <- printout[(length(roots) + 1):nrow(nodes), 1]

    # find dominant GO term (MF, CC, BP)
    if (over.represented == TRUE) {
        extendedGO <- extendGO(TopGOs)
        overrep <- ontologyOver(extendedGO)
        names(overrep) <- paste0("G", 1:nrow(nodes))

        bp <- names(which(overrep == "BP"))
        cc <- names(which(overrep == "CC"))
        mf <- names(which(overrep == "MF"))
        both <- overrep[sapply(overrep, length) == 2]
        bpmf <- names(both[which(apply
                                 ((sapply(seq_along(both),
                                          function(x) both[[x]] %in%
                                              c("BP", "MF"))),
                                 2, sum) == 2)])
        bpcc <- names(both[which(apply
                                 ((sapply(seq_along(both),
                                          function(x) both[[x]] %in%
                                              c("BP", "CC"))),
                                 2, sum) == 2)])
        ccmf <- names(both[which(apply
                                 ((sapply(seq_along(both),
                                          function(x) both[[x]] %in%
                                              c("CC", "MF"))),
                                 2, sum) == 2)])
        bpccmf <- names(overrep[sapply(overrep, length) == 3])
        funColors <- rep(1, (1 + nrow(printout)))
        funColors[which(as.character(new.order) %in%
                            bp)] <- 4
        funColors[which(as.character(new.order) %in%
                            cc)] <- 3
        funColors[which(as.character(new.order) %in%
                            mf)] <- 2
        funColors[which(as.character(new.order) %in%
                            bpmf)] <- "cyan4"
        funColors[which(as.character(new.order) %in%
                            ccmf)] <- 6
        funColors[which(as.character(new.order) %in%
                            bpcc)] <- "orchid4"
        funColors[which(as.character(new.order) %in%
                            bpccmf)] <- "goldenrod1"

        dendro <- assign_values_to_branches_edgePar(dend = dendro,
                                                                value = funColors, edgePar = "col")
        dendro %>% set("nodes_pch", 19) %>% set("nodes_col", "cadetblue") %>%
            set("leaves_pch", 19) %>% set("leaves_col", "coral2") %>%
            plot(main = "Dendrogram of differently expressed genes",
                 cex.main = 0.9, ylim = c(0, 1.2 * max(nodes[, 2])),
                 xlim = c(0, max(nodes[, 1])), ylab = "Height")
        text(nodes[, 1], nodes[, 2] + 0.02 * max(nodes[, 2]),
             label = new.order, cex = 0.7,
             pos = 3, srt = 90, font = 2)
        legend(title = "Dominant GO term", x = max(nodes[, 1]),
               y = max(nodes[, 2]), legend = c(" BP", " CC", " MF", " BP, MF",
                                               " BP, CC", " CC, MF",
                                               " BP, CC, MF", " No function"),
               xjust = 1, yjust = 0.5, x.intersp = 0.01,
               y.intersp = 0.8, cex = 0.7, lty = 1,
               col = c(4, 3, 2, "cyan4", 6, "orchid4",
                       "goldenrod1", 1))
    } else {
        dendro <- assign_values_to_branches_edgePar(dend = dendro,
                                                                value = 1,
                                                                edgePar = "col")
        dendro %>% set("nodes_pch", 19) %>% set("nodes_col", "cadetblue") %>%
            set("leaves_pch", 19) %>%
            set("leaves_col", "coral2") %>%
            plot(main = "Dendrogram of differently expressed genes",
                 cex.main = 0.9,
                 ylim = c(0, 1.2 * max(nodes[, 2])),
                 xlim = c(0, max(nodes[, 1])),
                 ylab = "Height")
        text(nodes[, 1], nodes[, 2] + 0.02 * max(nodes[, 2]),
             label = new.order, cex = 0.7,
             pos = 3, srt = 90, font = 2)
    }
    invisible()
}

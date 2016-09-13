#' GOpro: find the most characteristic gene ontology terms for groups of genes
#'
#' Finds the most characteristic gene ontology terms for groups of genes
#'
#' @docType package
#' @name GOpro
#' @useDynLib GOpro
#' @import AnnotationDbi
#' @import dendextend
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import org.Hs.eg.db
#' @import GO.db
#' @importFrom  Rcpp sourceCpp
#' @importFrom graphics legend plot text
#' @importFrom stats TukeyHSD aggregate aov as.dendrogram dist hclust na.omit p.adjust
NULL

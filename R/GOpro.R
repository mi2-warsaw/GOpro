#' GOpro: find the most characteristic gene ontology terms for groups of genes
#'
#' Based on the gene expressions find the structure somewhat
#' comparable to gene signature. From all given genes,
#' determine which are significantly different between sets.
#' These sets may relate to different health conditions of patients,
#' i.e. different types of cancer.
#' Then divide interesting genes into subsets.
#' Genes belong to a particular subset if they share the same feature.
#' There are two implemented methods that can be used to create
#' genes' subsets.
#' The first method is so-called all pairwise comparisons by
#' Tukey's procedure.
#' Genes that have the same profile (a result of all comparisons)
#' are assigned to one subset.
#' The second way of determining subsets is a method of
#' hierarchical clustering.
#' When all genes are divided into subsets, then
#' for each subset all relevant GO terms are searched for
#' in org.Hs.eg.db database.
#' Each found GO terms is tested using Fisher's test to find out
#' which of them are the most characteristic for the given subset of genes.
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
#' @importFrom stats TukeyHSD aggregate aov
#'  as.dendrogram dist hclust na.omit p.adjust
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges CharacterList NumericList
#' @seealso {\link{findGO}}
NULL

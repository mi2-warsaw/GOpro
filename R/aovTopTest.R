#' Find all significantly different genes based on their
#'  normalized expressions' values.
#' @noRd
#' @param expr A list of data frames of gene expressions
#'  with genes in rows and samples in columns.
#'  Rows must be named with genes' aliases.
#'  The order of samples and genes has to be the same for each data frame.
#' @param top A maximum number of genes to be returned for each group.
#' @param sig.level A numeric value, a significance level used in BH
#'  correction for multiple testing. The default value is 0.05.
#' @param parallel A logical value indicating if a task should
#'  be run on more than one core.
#'  The default value is \code{FALSE}.
#'
#' @return A list of significantly different genes and values
#'  of their expressions.
#'
#' @examples
#' aovTopTest(exrtcga, top = 50, sig.level = 0.1, parallel = FALSE)
#'
#'
#' @import doParallel
#' @import foreach
#' @import parallel

aovTopTest <- function(expr, top, sig.level = 0.05, parallel = FALSE) {

    aovCore <- function(z, expr, group) {
        ex <- unlist(sapply(expr, function(x) x[z, ]))
        gene <- cbind.data.frame(ex, group)
        if (length(ex) == length(unique(group))) {
            warning(paste0("The number of observations within
                           groups is too small in the row number ", z, "."))
            return(NA)
        }
        # at least 2 groups with at least 10% of
        # observations in the smallest group
        if (sum((tapply(gene$ex, gene$group, length) -
                 tapply(gene$ex, gene$group, function(x)
                     sum(is.na(x))) < 0.1 *
                 min(tapply(gene$ex, gene$group, length)))) >
            (length(expr) - 2)) {
            warning("AOV test: too many missing values.")
            NA
        } else {
            (summary(aov(ex ~ group, gene, na.action = na.omit))[[1]][[5]])[1]
        }
    }

    aovTest <- function(expr, parallel) {
        if (is.null(names(expr)) | any(names(expr) == "")) {
            group.names <- paste0("group", seq_along(expr))
        } else {
            group.names <- names(expr)
        }
        nc <- sapply(expr, ncol)
        if (is.null(nc)) {
            stop("The data are insufficient.
                 Provide observations for more than one gene.")
        }
        group <- rep(group.names, nc)
        z <- 0
        if (parallel == TRUE) {
            cores <- detectCores() - 1
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            out_a <- foreach(z = seq_len(nrow(expr[[1]])),
                             .export = "aovCore") %dopar%
                aovCore(z, expr, group)
            stopCluster(cl)
            out_a

        } else {
            foreach(z = seq_len(nrow(expr[[1]]))) %do%
                aovCore(z, expr, group)
        }
    }

    findGeneNames <- function(expr, parallel) {
        if (any(sapply(expr, function(x) is.null(rownames(x)))) ||
            !all(sapply(expr, function(x)
                sapply(expr, function(y) row.names(x) == row.names(y))))) {
            stop("Not all genes are named or not all genes
                 are in the same order for all data frames,
                 or some genes are missing.")
        }
        if (missing(expr) | (is.list(expr) & !is.null(dim(expr))) |
            (length(expr) < 2)) {
            stop("At least two groups are needed for this comparison.")
        }
        if (!is.numeric(unlist(expr))) {
            stop("expr can only contain numbers.")
        }
        if (!is.logical(parallel)) {
            parallel <- FALSE
        }
        anova.result.temp <- aovTest(expr, parallel)
        if (any(is.na(unlist(anova.result.temp)))) {
            anova.result <- anova.result.temp[-which(
                is.na(unlist(anova.result.temp)))]
            names(anova.result) <- row.names(expr[[1]][-which(
                is.na(unlist(anova.result.temp))), ])
        } else {
            anova.result <- anova.result.temp
            names(anova.result) <- row.names(expr[[1]])
        }
        anova.result
    }

    # sig.level and top parameters checking
    if (!is.numeric(sig.level) | (sig.level < 0)) {
        sig.level <- 0.05
        warning("No correct significance level specified.
                Taking 0.05 as the sig.level value.")
    }
    if (length(top) > 1) {
        top <- top[1]
        warning(paste("Too many values provided with the top parameter.
                      Taking ", top, " as the top value."))
    }
    if (top%%1 == 0) {
        top <- round(top)
    }
    if (missing(top) | is.null(top) | (top < 1)) {
        top <- nrow(expr[[1]])
        warning(paste("top parameter is missing or is not positive, taking ",
                      top, " as the top value."))
        top <- nrow(expr[[1]])
    }

    is.equal <- function(top) {
        if (top == sum(p.adj < sig.level)) {
            return(top)
        }
        if (top > sum(p.adj < sig.level)) {
            return(sum(p.adj < sig.level))
        }
        if (sorted[top] == sorted[top + 1]) {
            top <- top + 1
            top <- is.equal(top)
        }
        top
    }

    anova.result <- na.omit(findGeneNames(expr, parallel))
    p.adj <- p.adjust(anova.result, method = "fdr")
    sorted <- sort(p.adj[p.adj < sig.level])
    diff.expressed <- which(row.names(expr[[1]]) %in%
                                names(sorted[seq_len(is.equal(top))]))
    if (length(sorted) == 0) {
        message("Found no significantly different observations.")
    } else {
        lapply(expr, function(x) x[diff.expressed, ])
    }
}

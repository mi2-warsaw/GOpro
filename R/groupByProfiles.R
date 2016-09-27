#' Create profiles names of Tukey test results.
#' @noRd
#' @param part.tukey.res A vector of a list of
#' a result returned by \code{TukeyHSDTest}.
#' @param sig.level A numeric significance level of a test.
#' @param n A number of compared groups.
#'
#' @return A character vector of a name denoting
#' results of all pairwise Tukey comparisons.

makeNames <- function(part.tukey.res, sig.level, n) {
    means <- part.tukey.res[1:n]
    pvals <- part.tukey.res[(n + 1):length(part.tukey.res)]
    names <- names(pvals)

    if (IsConflict(part.tukey.res, sig.level, n)) {
        signedNames <- sapply(names, function(x) {
            temp <- unlist(strsplit(x, "-"))
            means_temp <- sapply(temp, function(z) {
                means[which(names(means) == z)]
            })
            if (pvals[which(names == x)] > sig.level) {
                sign <- "="
            } else {
                sign <- ifelse(means_temp[1] > means_temp[2],
                               ">", "<")
            }
            paste(temp, collapse = sign)
        })
        paste(c(paste(c(signedNames, "ascending means: "), collapse = ";"),
                paste(names(sort(means)), collapse = ",")), collapse = "")
    } else {
        name <- names(means)[1]
        whichpairs <- pvals[sapply(names, function(x) {
            temp <- unlist(strsplit(x, "-"))
            if (temp[1] == name || temp[2] == name) {
                TRUE
            } else {
                FALSE
            }
        })]

        relations <- sapply(names(whichpairs), function(x) {
            if (whichpairs[which(names(whichpairs) == x)] > sig.level) {
                return("equal")
            }
            temp <- unlist(strsplit(x, "-"))
            nameInd <- which(temp == name)
            secondInd <- nameInd%%2 + 1
            tempName1 <- temp[nameInd]
            tempName2 <- temp[secondInd]
            mean1 <- means[which(names(means) == tempName1)]
            mean2 <- means[which(names(means) == tempName2)]
            if (mean1 > mean2) {
                "lower"
            } else {
                "greater"
            }
        })
        lowers <- whichpairs[which(relations == "lower")]
        greaters <- whichpairs[which(relations == "greater")]
        equals <- whichpairs[which(relations == "equal")]

        equalnames <- unique(unlist(strsplit(names(equals), "-")))
        lowernames <- unique(unlist(strsplit(names(lowers), "-")))
        lowernames <- lowernames[-which(lowernames == name)]
        greaternames <- unique(unlist(strsplit(names(greaters),
                                               "-")))
        greaternames <- greaternames[-which(greaternames == name)]
        if (length(equals) > 0) {
            equalmeans <- sort(means[sapply(names(means), function(x) {
                if (length(which(equalnames == x)) == 0) {
                    FALSE
                } else {
                    TRUE
                }
            })])
            output <- paste(names(equalmeans), collapse = "=")
        } else {
            output <- name
        }
        if (length(greaternames) == 1) {
            output <- paste(c(output, greaternames), collapse = "<")
        } else if (length(greaternames) > 1) {
            greatermeans <- means[sapply(names(means), function(x) {
                if (length(which(greaternames == x)) == 0) {
                    FALSE
                } else {
                    TRUE
                }
            })]

            greaterpvals <- sapply(names(pvals), function(x) {
                temp <- unlist(strsplit(x, "-"))
                if (length(which(greaternames == temp[1])) ==
                    0 || length(which(greaternames == temp[2])) ==
                    0) {
                    FALSE
                } else {
                    TRUE
                }
            })
            output <- paste(c(output, makeNames(c(greatermeans,
                                                  pvals[greaterpvals]),
                                                sig.level,
                                                length(greatermeans))),
                            collapse = "<")
        }
        if (length(lowernames) == 1) {
            output <- paste(c(lowernames, output), collapse = "<")
        } else if (length(lowernames) > 1) {
            lowermeans <- means[sapply(names(means), function(x) {
                if (length(which(lowernames == x)) == 0) {
                    FALSE
                } else {
                    TRUE
                }
            })]

            lowerpvals <- sapply(names(pvals), function(x) {
                temp <- unlist(strsplit(x, "-"))
                if (length(which(lowernames == temp[1])) == 0 ||
                    length(which(lowernames == temp[2])) == 0) {
                    FALSE
                } else {
                    TRUE
                }
            })
            output <- paste(c(makeNames(c(lowermeans,
                                          pvals[lowerpvals]),
                                        sig.level,
                                        length(lowermeans)),
                              output),
                            collapse = "<")
        }
        output
    }
}

#' Check if there is a conflict between results of the Tukey test,
#' ex. a=b, b=c, c!=a.
#' @noRd
#' @param part.tukey.res A vector of a list of a result
#' returned by \code{TukeyHSDTest}.
#' @param sig.level A numeric significance level of a test.
#' @param n A number of compared groups.
#'
#' @return A logical value.

IsConflict <- function(part.tukey.res, sig.level,
                       n) {
    means <- part.tukey.res[1:n]
    pvals <- part.tukey.res[(n + 1):length(part.tukey.res)]

    for (i in 1:(n - 1)) {
        group1 <- names(means)[i]
        for (j in (i + 1):n) {
            group2 <- names(means)[j]
            pair_str1 <- paste(c(group1,
                                 group2), collapse = "-")
            pair_str2 <- paste(c(group2,
                                 group1), collapse = "-")
            pair_position <- max(which(names(pvals) == pair_str1),
                                 which(names(pvals) == pair_str2))

            mean1 <- means[which(names(means) == group1)]
            mean2 <- means[which(names(means) == group2)]
            for (k in 1:n) {
                if (i != k && j != k) {
                    group3 <- names(means)[k]
                    mean3 <- means[which(names(means) == group3)]

                    pair_str1 <- paste(c(group1, group3), collapse = "-")
                    pair_str2 <- paste(c(group3, group1), collapse = "-")

                    pair_position2 <- max(which(names(pvals) == pair_str1),
                                          which(names(pvals) == pair_str2))

                    pair_str1 <- paste(c(group2, group3), collapse = "-")
                    pair_str2 <- paste(c(group3, group2), collapse = "-")

                    pair_position3 <- max(which(names(pvals) == pair_str1),
                                          which(names(pvals) == pair_str2))

                    if ((pvals[pair_position] > sig.level)) {
                        if ((pvals[pair_position2] < sig.level &&
                             pvals[pair_position3] > sig.level) ||
                            (pvals[pair_position2] > sig.level &&
                             pvals[pair_position3] < sig.level) ||
                            ((pvals[pair_position2] < sig.level &&
                              pvals[pair_position3] < sig.level) &&
                             ((mean1 > mean3 && mean2 < mean3) ||
                             (mean1 < mean3 && mean2 > mean3)))) {
                            return(TRUE)
                        }
                    }
                }
            }
        }
    }
    FALSE
}

#' Genes grouped by profiles determined by the result
#' returned from \code{tukeyHSDTest} function.
#' @noRd
#' @param expr A list of data frames obtained from \code{aovTopTest}.
#' @param tukey.results A data frame returned by \code{TukeyHSDTest}.
#' @param sig.level A numeric value, a significance level used in
#' Tukey's all pairwise comparisons (\code{groupByTukey}).
#'
#' @return A list of groups of genes named by profiles.

groupByProfiles <- function(expr, tukey.results, sig.level) {
    n <- length(expr)
    named.profiles <- sapply(tukey.results, function(x)
        makeNames(x, sig.level, n))
    groups.names <- sort(unique(named.profiles))
    result.groups <- lapply(groups.names, function(x)
        row.names(expr[[1]][which(named.profiles == x), ]))
    names(result.groups) <- groups.names
    result.groups
}

#' @export
#' @importFrom methods is
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom S4Vectors split mcols
#' @importFrom IRanges IRanges
homopolymerFinder <- function(seq, specific, min.homo=2, report0=FALSE)
# Finds homopolymers in a given DNAStringSet object.
{
    if (!is(seq, "DNAStringSet")) {
        stop("'seq' should be a DNAStringSet object")
    }
    all.homo <- .Call(sarlacc:::cxx_find_homopolymers, seq)
    homo.range <- IRanges(all.homo[[2]], all.homo[[2]] + all.homo[[3]] - 1L)
    mcols(homo.range)$base <- all.homo[[4]]
    groupings <- factor(all.homo[[1]] + 1L, levels = seq_along(seq))
    if (!missing(specific)) {
        group.count <- table(groupings)
        match.base <- elementMetadata(homo.range)$base==rep(specific,group.count)
        homo.range <- homo.range[match.base]
        groupings <- groupings[match.base]
    }
    long.enough <- width(homo.range)>=min.homo
    homo.range <- homo.range[long.enough]
    groupings <- groupings[long.enough]
    output <- split(homo.range, groupings, drop = FALSE)
    if (report0 == TRUE) {
        output[lengths(output) == 0] <- IRanges(start = 0, width = 0)
    }
    names(output) <- names(seq)
    return(output)
}


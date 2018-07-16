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

    all.homo <- .Call(cxx_find_homopolymers, seq)
    homo.range <- IRanges(all.homo[[2]], all.homo[[2]] + all.homo[[3]] - 1L)
    mcols(homo.range)$base <- all.homo[[4]]

    groupings <- factor(all.homo[[1]]+1L, levels=seq_along(seq))
    output <- split(homo.range, groupings, drop=FALSE)
    
    if(!missing(specific)){
        output <- lapply(output, FUN=function(x)if(length(x)!=0)x[elementMetadata(x)$base%in%specific] else x)
    }
    
    if(report0 == TRUE){
        output[lengths(output)==0] <- IRanges(start=0, width=0) 
    }
    output <- lapply(output, FUN=function(x)if(width(x)[1]>0)x[width(x)>=min.homo] else x)
    
    names(output) <- names(seq)
    return(output)
}


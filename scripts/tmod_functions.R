#' message + sprintf
smessage <- function(...) message(sprintf(...))



#' tmod_extract_values
#'
#' from a set of intermediate results from DESeq2, extract p-values, log2
#' fold changes and feature_ids. log2 and fold changes preserve the order
#' from the files, while feature_ids are ordered by the respective key.
#' MSD, minimal significant difference, is a measure "between" p-value and
#' log2fold change, very useful to order the genes by the magnitude of the
#' observed effect.
tmod_extract_values <- function(files, path="intermediate_results", extension="_results.rds", sort_by="msd") {

  sort_by <- match.arg(sort_by, c("msd", "pvalue", "log2FoldChange"))

  res <- map(files, ~ {
    fp <- file.path(path, paste0(.x, extension))
    readRDS(fp) %>% as.data.frame %>% rownames_to_column("feature_id") %>% as_tibble %>%
      mutate(msd=abs(log2FoldChange) - 1.96 * lfcSE)
  })

  names(res) <- files

  pvals <- map_dfr(res, ~ .x %>% pull(pvalue) %>% replace_na(1))
  lfcs  <- map_dfr(res, ~ .x %>% pull(log2FoldChange) %>% replace_na(0))
  
  genes <- switch(sort_by, 
    msd=    map(res, ~ .x %>% arrange(desc(msd)) %>% pull(feature_id)),
    pvalue= map(res, ~ .x %>% arrange(pvalue) %>% pull(feature_id)),
    log2FoldChange= map(res, ~ .x %>% mutate(al=abs(log2FoldChange)) %>% arrange(desc(al)) %>% pull(feature_id))
    )

  list( 
    feature_ids=res[[1]]$feature_id,
    pvals=pvals,
    lfcs=lfcs,
    genes=genes)
}


#' Static cache mechanism
#'
#' Checks for the existence of the specified file. If it exists, it is
#' loaded with readRDS. If not, expr is evaluated and the result is stored in
#' the specified path using saveRDS. In any case, the object stored is
#' returned.
#' @param path path to an RDS file
#' @param expr expression which is the "recipe" for the object
#' @param description (optional) name of the variable to report in messages
#' @param quiet if TRUE, no messages are shown
tmod_static_cache <- function(path, expr, description="unnamed", quiet=FALSE) {

  if(file.exists(path)) {
    ret <- readRDS(path)
    if(!quiet) { 
      message(sprintf("%s: Loading existing version from file %s, created on %s", description, path, attr(ret, "creation time")))
    }
  } else {
    if(!quiet) {
      message(sprintf("%s: Generating object and saving to file %s", description, path))
    }

    ret <- expr
    attr(ret, "creation time") <- Sys.time() 
    saveRDS(ret, file=path)
  }

  return(ret)
}

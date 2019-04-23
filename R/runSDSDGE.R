#'run SDSDGE
#' @param object Seurat object
#' @import dplyr tidyr Seurat
#' @export
#'
runSDSDGE <- function(object){
  DGE <- list()
  for(c in names(object@misc$sds$data@curves)){
    object@misc$sds$dge[[c]] <- FetchData(object,append(object@var.genes, c,0)) %>% tidyr::gather(gene,signal, -one_of(c)) %>% dplyr::rename(curve = 1) %>%
      tidyr::nest(-gene) %>%
      mutate(
        fit = purrr::map(data, ~ gam(signal ~ lo(curve), data = .x)),
        tidied = purrr::map(fit, tidy)
      ) %>%
      tidyr::unnest(tidied) %>%
      filter(term !='Residuals')
  }
  return(object)

}

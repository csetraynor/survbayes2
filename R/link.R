#' GAMMA poisson link
#'
#' Link model to linear predictor to individual
#'  survival predictions.
#'
#' @param post posterior draws \cr
#' @param strata strata
#' @param obs original dataframe defaults to ic2surv
#' @return d a plot dataset
#' @export
#' @importFrom magrittr %>%
link.surv <- function(post, longdata){
  plot.list <- lapply(1:nrow(post), function(i){
    longdata$loghaz <- post[i, ]
    longdata %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
      dplyr::filter(row_number()==n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(time) %>%
      dplyr::select(surv)
  })
  plot.matrix <- do.call(cbind, plot.list)
  as.matrix(plot.matrix)
}

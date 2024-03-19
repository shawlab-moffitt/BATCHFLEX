#' BatchFLEX_shiny
#'
#' Used to run a shiny version of `Batch_FLEX` locally.
#'
#' @return Message indicating BatchFLEX is starting
#' @export
#'
#' @examples
#' set.seed(333)
#' BatchFLEX_shiny()
#'
BatchFLEX_shiny <- function(){
  shiny::runGitHub(
    repo = "BATCH-FLEX-shinyApp",
    username = "shawlab-moffitt"
  )
  return(message("Starting Batchflex"))
}

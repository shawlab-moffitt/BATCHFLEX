#' BatchFLEX_shiny
#'
#' @return Message indicating BatchFLEX is starting
#' @export
#'
#' @examples
#' set.seed(333)
BatchFLEX_shiny <- function(){
  shiny::runGitHub(
    repo = "BATCH-FLEX-shinyApp",
    username = "shawlab-moffitt"
  )
  return(message("Starting Batchflex"))
}

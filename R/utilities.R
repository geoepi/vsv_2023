# load all .R functions/scripts in the /R directory
source_dir <- function(directory_path) {
  
  if (!dir.exists(directory_path)) {
    stop("The specified directory does not exist.")
  }
  
  file_list <- list.files(path = directory_path, pattern = "\\.R$", full.names = TRUE)
  
  for (file in file_list) {
    source(file)
  }
}

library(knitr)
knit_print.gt <- function(x, ...) {
  stringr::str_c(
    "<div style='all:initial';>\n", 
    gt::as_raw_html(x), 
    "\n</div>"
  ) |> 
    knitr::asis_output()
}
registerS3method(
  "knit_print", 'gt_tbl', knit_print.gt, 
  envir = asNamespace("gt") 
)
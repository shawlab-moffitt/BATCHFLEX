# library(tidyverse)
# #example meta data
# meta = readRDS("data-raw/meta.rds")
# example_meta = meta %>%
#   group_by(Study, CellType) %>%
#   slice(1)
#
# #example matrix data
# mat = readRDS("data-raw/example_matrix.rds")
#
# #there are duplicated columns in the expression matrix for some reason?
# duped_cols = colnames(mat) %>% .[duplicated(.)]
# example_meta = example_meta %>%
#   filter(!(Name %in% duped_cols))
# example_mat = mat[,colnames(mat) %in% example_meta$Name]


usethis::use_data(example_meta, overwrite = TRUE)
usethis::use_data(example_mat, overwrite = TRUE)

# BATCH-FLEX

An R package for correct batch effects in numeric data (typically RNA sequencing data). Companion to the BatchFlex Shiny Application at:

[shiny link]()

To install package:

install from CRAN:

install.packages("devtools")
install.packages("BiocManager")

```
install.packages("devtools")
install.packages("BiocManager")
#stringi
install.package("https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.3/stringi_1.8.3.tgz",
repos = NULL, type = 'source')

#https://bioconductor.org/packages/release/bioc/html/sva.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/sva_3.50.0.tgz",
repos = NULL, type="source")
#https://bioconductor.org/packages/release/bioc/html/limma.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/limma_3.58.1.tgz",
repos = NULL, type = 'source')
#https://bioconductor.org/packages/release/bioc/html/RUVSeq.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/RUVSeq_1.36.0.tgz",
repos = NULL, type = 'source')
#https://github.com/omnideconv/immunedeconv
devtools::install_github("omnideconv/immunedeconv")

devtools::install_github('shawlab-moffitt/BATCHFLEX')
```



#' Gene expression matrix
#'
#' Data used in the vignette of cit_gsa. The dataset contain normalized single-
#' cell RNA-seq gene expression from Dendritic Cells (DCs).
#'
#'
#' @format
#' A data frame of numeric values with 7458 rows (genes) and 539 columns (cells).
#' 
#' Row names are the genes symbol and column names are cells identifier. 
#' 
#' 
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820>
#' @references Villani, A. C., Satija, R., Reynolds, G., Sarkizova, S., Shekhar,
#'  K., Fletcher, J., ... & Hacohen, N. (2017). Single-cell RNA-seq reveals new 
#'  types of human blood dendritic cells, monocytes, and progenitors. Science, 
#'  356(6335), eaah4573.
#'  
"data_DCMono_discovery"




#' DC subpopulation
#'
#' DC subpopulation for each cells of the gene expression matrix.
#' 
#'
#' @format
#' A vector of character of length 539 that is the number of cells in the 
#' expression matrix. The possible values are: DC1, DC2 & DC3 and pDC.
#' 
#' 
#' 
#' @source <https://www.science.org/doi/10.1126/science.aah4573#supplementary-materials>
#' @references Villani, A. C., Satija, R., Reynolds, G., Sarkizova, S., Shekhar,
#'  K., Fletcher, J., ... & Hacohen, N. (2017). Single-cell RNA-seq reveals new 
#'  types of human blood dendritic cells, monocytes, and progenitors. Science, 
#'  356(6335), eaah4573.
#'  
"DCpop"



#' Gene set data base 
#'
#' Blood Transriptional Modules of Stanford. 
#' 
#'
#' @format
#' A list that contain 2 elements:
#' \describe{
#'  \item{genesets}{A list where each element contain the genes of the gene set}
#'  \item{geneset.names}{A character vector that contain the gene set names}
#' }  
#' 
#' 
#' 
#' @source <https://www.nature.com/articles/ni.2789#Sec25>
#' @references Li, S., Rouphael, N., Duraisingham, S., Romero-Steiner, S., 
#' Presnell, S., Davis, C., ... & Pulendran, B. (2014). Molecular signatures of 
#' antibody responses derived from a systems biology study of five human 
#' vaccines. Nature immunology, 15(2), 195-204.
#'  
"BTM"































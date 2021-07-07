#'@param
#'E: Expression
#'O: Openness
#'Symbol: Symbol Name
#'
"RegNMF" <- function(E,
                   O,
                   Symbol,
                   PeakName,
                   Symbol_location,
                   Peak_location,
                   feature_cut_perc=0.01,
                   corr_cut_k=100000,
                   cell_num_smooth=sqrt(dim(E)[2]),
                   core=8){

  UseMethod("RegNMF");
}

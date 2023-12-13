#------------------------------------------------
#' @title Convert vcf to long format
#'
#' @description Convert a vcf into a long format data frame with sample ID, locus, alleles and read counts for each allele. 
#'
#' @param vcf object of class vcfR
#'
#' @importFrom vcfR extract.gt is.biallelic
#' @importFrom tibble rownames_to_column 
#' @importFrom tidyr pivot_longer unnest
#' @importFrom dplyr rowwise mutate group_by relocate n
#' @importFrom stringr str_split
#' @importFrom rlang .data
#' @export
#' @examples
#' 

vcf2long <- function(vcf) {
  
  # check inputs
  assert_class(vcf, "vcfR")
  
  # print message to console
  message("Converting from vcf to long format...")
  
  # extract allele counts
  ad <- t(extract.gt(vcf, element = 'AD'))
  
  # make df and into long format
  counts_df <- ad |> 
    as.data.frame() |> 
    rownames_to_column("sample_id") |> 
    pivot_longer(cols = -.data$sample_id, names_to = "locus", values_to = "read_count")
  
  # unnest read_count
  long_df <- counts_df |>
    rowwise() |>
    # split all read count values
    mutate(read_count = list(str_split(.data$read_count, ",")[[1]])) |>
    unnest(cols = .data$read_count) |> 
    group_by(.data$sample_id, .data$locus) |> 
    # create new variable 'allele'
    mutate(allele = paste0("allele-", rep(1:n()))) |> 
    relocate(.data$allele, .before = .data$read_count)

  message("Reformatting complete.")
  
  # Check if any loci are not biallelic and record how many
  n_not_biallelic <- length(which(!is.biallelic(vcf)))
  
  # If the vcf is not biallelic, display a warning message
  if(n_not_biallelic > 0){
  warning("Your vcf is not all bi-allelic. Make sure to double check if this is not expected.")
  }
  
  return(long_df)
}

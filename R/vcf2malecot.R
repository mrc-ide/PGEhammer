#####LOAD FUNCTIONS#####

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


#' Convert VCF file to haplotype table
#'
#' This function takes a VCF file as input and returns a haplotype table.
#' The function determines if the data is biallelic or multiallelic and
#' returns the haplotype table in either wide or long format accordingly.
#'
#' @param vcf_file The path to the VCF file.
#'
#' @return A data frame representing the haplotype table.
#' @export
#'
#' @examples
#' vcf_file <- "/path/to/your/vcf/file.vcf.gz"
#' hap_redux <- vcf2malecot(vcf_file)
#'
vcf2malecot <- function(vcf_file) {
  if (!require(geneHapR)) {
    install.packages("geneHapR")
    library(geneHapR)
  }
  
  if (!require(tidyr)) {
    install.packages("tidyr")
    library(tidyr)
  }
  
  vcf = import_vcf(vcf_file)
  
  ref = vcf@fix[,'REF']
  alt = vcf@fix[,'ALT']
  
  hap = vcf2hap(vcf,
                hapPrefix = "H",
                hetero_remove = FALSE,
                na_drop = FALSE)
  
  locus = paste0(hap[1,3:ncol(hap)-1], "_", hap[2,3:ncol(hap)-1])
  accession = hap$Accession[-c(1:4)]
  
  hap_redux = hap[-c(1,2,3,4), -c(1,ncol(hap))]
  colnames(hap_redux) = locus
  rownames(hap_redux) = accession
  
  if(any(grepl(",", alt))) {
    print("Multi-allelic data - will return haplotype table in long format")

    if (!require(checkmate)) {
      install.packages("checkmate")
      library(checkmate)
    }
    
    if (!require(vcfR)) {
      install.packages("vcfR")
      library(vcfR)
    }
    
    if (!require(tidyverse)) {
      install.packages("tidyverse")
      library(tidyverse)
    }
    
    vcf = read.vcfR(vcf_file)
    vcf_final = vcf2long(vcf)
    
    vcf_final$read_count = as.integer(vcf_final$read_count)
    
    vcf_final = vcf_final %>%
      group_by(sample_id, locus) %>%
      mutate(read_count_total = sum(read_count),
             allele_freq = read_count / read_count_total,
             filter = allele_freq > 0.05) %>%
      ungroup()
    
    filtered_vcf_final <- vcf_final %>%
      filter(filter)  %>%
      select(sample_id, locus, allele) %>%
      mutate(allele = as.integer(factor(allele)),
             locus = as.integer(factor(locus)))
    
    vcf_final_integer <- filtered_vcf_final %>%
      mutate(locus = as.integer(factor(locus)))
    
    hap_redux = vcf_final_integer
    colnames(hap_redux) = c('sample_ID', 'locus', 'haplotype')
    return(hap_redux)
    
  } else {
    print("Bi-allelic data - will return haplotype table in wide format")
    for(i in 1:ncol(hap_redux)) {
      locus = hap_redux[,i]
      locus_ref = ref[i]
      locus_alt = alt[i]
      het =  c(paste0(locus_ref, "|", locus_alt), paste0(locus_alt, "|", locus_ref))
      hap_redux[,i] = convert_het_to_numeric_bi(locus_ref, locus_alt, het, locus)
    }
    hap_redux$sample_ID <- rownames(hap_redux)
    hap_redux <- hap_redux[, c(ncol(hap_redux), 1:(ncol(hap_redux)-1))]
    rownames(hap_redux) <- NULL
  }
  return(hap_redux)
}

#------------------------------------------------
#' @title vcf to pedmap format
#'
#' @description Function to convert processed vcf data into pedmap format ready to use for isoRelate
#'
#' @param vcf_file file path to VCF file containing all variants (SNPs and indels)
#' @param gds_file file path to the output of the GDS file for SeqArray
#' @param ped_file file path to the output of ped file
#' @param map_file file path to the output of map file
#'
#' @importFrom SeqArray
#' @export
#' @examples
#' #vcf_to_gds("../SpotMalariapfPanel_simData_snponly_sanger100.vcf.gz", "../SpotMalariapfPanel_simData_snponly_sanger100.gds")  
#' #gds_to_pedmap("../SpotMalariapfPanel_simData_snponly_sanger100.gds", "my_map.map", "my_ped.ped")

globalVariables(c("%>%", "group_by", "chr", "pos", "summarise", "variant_id", "arrange", "mutate", "cM"))

vcf_to_gds <- function(vcf_file, gds_file) {
  
  # read in VCF file as VCF object and perform quality filtration
  vcf.file <- (vcf_file)
  gds.file <- (gds_file)
  seqVCF2GDS(vcf.file, gds.file)
}

gds_to_pedmap <- function(gds_file, ped_file, map_file) {
  showfile.gds(closeall = TRUE)
  snp_data <- seqOpen(gds_file, readonly = TRUE)
  variants <- data.frame(variant_id = seqGetData(snp_data, "variant.id"),
                         chr = seqGetData(snp_data, "chromosome"),
                         pos = seqGetData(snp_data, "position"),
                         alleles = seqGetData(snp_data, "allele")) %>%
    group_by(chr, pos) %>% summarise(variant_id=min(variant_id)) %>% 
    arrange(chr, pos) %>% as.data.frame
  
  variants$allele <- seqGetData(snp_data, "allele")
  samp <- data.frame(sample=seqGetData(snp_data, "sample.id"))
  samples <- as.character(samp$sample)
  
  n_var <- nrow(variants)
  n_samp <- length(samples)

  genotypes <- seqGetData(snp_data, "genotype") %>% matrix(ncol=n_var, nrow=n_samp)
  rownames(genotypes) <- samples
  
  ## recode genotypes so 1 = ref allele, 2 = alt allele, 0 = missing
  genotypes <- genotypes + 1
  genotypes[is.na(genotypes)] <- 0
  
  ## extract genotypes at biallelic loci
  allele_count <- apply(genotypes, 2, max)
  biallelic <- genotypes[, allele_count<=2]
  
  ## label biallelic loci
  colnames(biallelic) <- (variants[allele_count<=2, ] %>% mutate(name=paste(chr, ":", pos, sep="")))$name
  
  ## construct a map file
  my_map <- variants[allele_count<=2,] %>% mutate(cM=pos/17000) %>%
    dplyr::select(chr, variant_id, cM, pos)
  my_map <- na.omit(my_map)
  
  ## write a map file corresponding to biallelic snps
  write.table(my_map, file=map_file, row.names=FALSE, col.names = FALSE, sep="\t", quote=FALSE)
  
  ## write a ped file encoding genotypes at biallelic snps
  for (sample in samples) {
    cat(paste0("Country\t", sample, "\t0\t0\t1\t0\t", paste(rep(biallelic[sample,], each=2), collapse="\t"), "\n"), file=ped_file, append=TRUE)
  }
  
}

#'
#' Test single-cell allele-specific expression data for scLinaX.
#' @description
#'  Test single-cell allele-specific expression data for scLinaX. This data was generated from 10X multiome data using cellsnp-lite and annovar.
#'
#' @format ## `multiome_ASE_df`
#' A dataframe (tibble) containing single-cell allele-specific expression (scASE) data for 10X multiome data
#' \describe{
#'   \item{SNP_ID}{SNP identifier}
#'   \item{POS}{Genomic position of the SNP (GRCh38)}
#'   \item{REF,ALT,OTH}{Allelic expression of the reference/alterntive/other allele}
#'   \item{Gene}{Gene annotated to the SNP by Annovar}
#'   \item{Sample_ID}{Sample_ID}
#'   ...
#' }
#' @source <https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0>
"multiome_ASE_df"





#' A dataframe (tibble) containing X chromosome inactivation (XCI) status.
#' @description
#' This dataframe contains XCI status data, originally sourced from Tukiainen et al., Nature, 2017, with additional manual curation.
#'
#' @format ## `XCI_ref`
#' \describe{
#' A dataframe (tibble) with the following columns:
#'   \item{Gene}{Gene name}
#'   \item{XCI_status}{XCI status: escape, variable, inactive}
#'   ...
#'}
#' @source Tukiainen et al., Nature, 2017; <https://doi.org/10.1038/nature24265>
"XCI_ref"





#' Test data for QC of reference genes.
#' @description
#' This dataframe (tibble) was generated for the quality control (QC) of reference genes using the AIDA dataset. It was processed using the `run_RefGeneQC` function.
#'
#' @format ## `AIDA_QCREF`
#' \describe{
#' A dataframe (tibble) with the following columns:
#'   \item{Gene}{Gene name}
#'   \item{Mean_AR_target}{Average ratio of expression from Xi across other candidate reference genes when the SNPs on the gene were used as references}
#'   \item{Mean_AR_reference}{Average ratio of expression from Xi for the gene when SNPs on other candidate reference genes were used as references}
#'   ...
#'}
#' @source Asian Immune Diversity Network: <https://chanzuckerberg.com/science/programs-resources/single-cell-biology/ancestry-networks/immune-cell-atlas-of-asian-populations/>
"AIDA_QCREF"





#' Test data for QC of reference genes.
#' @description
#' 	A dataframe (tibble) for the annotation of the cells in the 10X multiome dataset. Annotation was performed with reference-based mapping with Azimuth <https://azimuth.hubmapconsortium.org/>.
#'
#' @format ## `multiome_Annotation`
#' \describe{
#' A dataframe (tibble) with the following columns:
#'   \item{cell_barcode}{Cell barcode}
#'   \item{Annotation}{Annotation of the cell}
#'   ...
#'}
#' @source <https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0>
"multiome_Annotation"





#' Test single-cell allele-specific chromatin accessibility data for scLinaX-multi.
#' @description
#' Test single-cell allele-specific chromatin accessibility data for scLinaX-multi. This data was generated from 10X multiome data using cellsnp-lite and annovar.
#'
#' @format ## `multiome_ATAC_ASE_data`
#' A dataframe (tibble) containing single-cell allele-specific chromatin accessibility data for 10X multiome data
#' \describe{
#'   \item{SNP_ID}{SNP identifier}
#'   \item{POS}{Genomic position of the SNP (GRCh38)}
#'   \item{REF,ALT,OTH}{Allelic expression of the reference/alterntive/other allele}
#'   \item{Sample_ID}{Sample_ID}
#'   ...
#' }
#' @source <https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0>
"multiome_ATAC_ASE_data"





#' Test peak data for scLinaX-multi.
#' @description
#' A dataframe (tibble) containing peak information obtained from single-cell chromatin accessibility data. This data was generated from 10X multiome data using Signac, following the tutorial at <https://stuartlab.org/signac/>.
#' \describe{
#' @format ## `multiome_peak_data`
#' A dataframe (tibble) with the following columns:
#'   \item{start}{Genomic position of the start of the peak}
#'   \item{end}{Genomic position of the end of the peak}
#'   \item{Peak_name}{Name of the peak}
#'   \item{Gene}{Nearest gene to the peak}
#'   ...
#'}
#' @source 10x Genomics: <https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0>
"multiome_peak_data"


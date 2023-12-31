% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RefGeneQC.R
\name{make_QC_gene_table}
\alias{make_QC_gene_table}
\title{Create a dataframe for QC of candidate reference genes}
\usage{
make_QC_gene_table(
  QCed_df,
  per_SNP_sample,
  HE_allele_cell_number_THR = 50,
  QC_total_allele_THR = 10
)
}
\arguments{
\item{QCed_df}{\verb{Output from make_QCed_df} function}

\item{per_SNP_sample}{Output from \code{make_per_SNP_sample_df} function}

\item{HE_allele_cell_number_THR}{Threshold for the number of cells expressing reference SNPs. Candidate reference SNPs expressed in at least "HE_allele_cell_number_THR" cells are used for the analysis. Default: 50.}

\item{QC_total_allele_THR}{Threshold for the total allele count (depth) of the SNP used for calculating the ratio of expression from Xi.
Note that this count is calculated with cells successfully assigned to the group based on the inactivated X chromosome.
This filter is applied in the final step of scLinax and differs from "SNP_DETECTION_DP". Default: 10.}
}
\description{
Create a dataframe for QC of candidate reference genes
}

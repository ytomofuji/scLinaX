% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scLinaX.R
\name{run_scLinaX}
\alias{run_scLinaX}
\title{Run scLinaX analysis for all samples}
\usage{
run_scLinaX(
  ASE_df,
  XCI_ref,
  QCREF,
  Inactive_Gene_ratio_THR = 0.05,
  SNP_DETECTION_DP = 30,
  SNP_DETECTION_MAF = 0.1,
  QC_total_allele_THR = 10,
  HE_allele_cell_number_THR = 50,
  REMOVE_ESCAPE = TRUE,
  PVAL_THR = 0.01,
  RHO_THR = 0.5
)
}
\arguments{
\item{ASE_df}{A dataframe (tibble) containing single-cell allele-specific expression (scASE) data for all samples. This dataframe should have the following columns:
\itemize{
\item SNP_ID: SNP identifier
\item POS: Genomic position of the SNP (GRCh38)
\item REF: Reference allele of the SNP {A,T,G,C}
\item ALT: Alternative allele of the SNP {A,T,G,C}
\item cell_barcode: Cell barcode
\item REFcount: Allelic expression of the reference allele
\item ALTcount: Allelic expression of the alternative allele
\item OTHcount: Allelic expression of the other allele
\item Sample_ID: Sample ID
\item Gene: Gene annotated to the SNP
}}

\item{XCI_ref}{A dataframe (tibble) containing X chromosome inactivation status. This dataframe should have the following two columns:
\itemize{
\item Gene: Gene name
\item XCI_status: XCI status {escape, variable, inactive}
}}

\item{QCREF}{A dataframe generated for reference gene QC (please see the \code{run_RefGeneQC} function). QCREF should have the following three columns:
\itemize{
\item Gene
\item Mean_AR_target
\item Mean_AR_reference
}}

\item{Inactive_Gene_ratio_THR}{Threshold for the ratio of inactive genes}

\item{SNP_DETECTION_DP}{Threshold for the total allele count (depth) of the SNP in the scASE data. SNP–Sample pairs with a total allele count of at least "SNP_DETECTION_DP" are used for the analysis. Default: 30.}

\item{SNP_DETECTION_MAF}{Threshold for the minor allele count of the scASE data. SNP–Sample pairs with a minor allele ratio between "SNP_DETECTION_MAF" and "1 - SNP_DETECTION_MAF" are used for the analysis. Default: 0.1.}

\item{QC_total_allele_THR}{Threshold for the total allele count (depth) of the SNP used for calculating the ratio of expression from Xi.
Note that this count is calculated with cells successfully assigned to the group based on the inactivated X chromosome.
This filter is applied in the final step of scLinax and differs from "SNP_DETECTION_DP". Default: 10.}

\item{HE_allele_cell_number_THR}{Threshold for the number of cells expressing reference SNPs. Candidate reference SNPs expressed in at least "HE_allele_cell_number_THR" cells are used for the analysis. Default: 50.}

\item{REMOVE_ESCAPE}{Boolean:
\itemize{
\item TRUE (default): Remove ASE profiles of SNPs on escapee genes when calculating correlations between pseudobulk ASE profiles.
\item FALSE: Include ASE profiles of SNPs on escapee genes when calculating correlations between pseudobulk ASE profiles.
}}

\item{PVAL_THR}{Threshold for P-values and absolute correlation coefficients in the correlation analysis of pseudobulk profiles generated for reference SNPs.}

\item{RHO_THR}{Threshold for absolute correlation coefficients in the correlation analysis of pseudobulk profiles generated for reference SNPs.}
}
\value{
A list of objects ($result, $raw_exp_result, $Fail_list) representing the results of scLinaX.
Additional objects representing intermediate results of scLinaX are also included in the list ($clustering_result, $Max_Num_Table_result, $df_snp_summary, $phasing_result).
While users may not typically need to use these intermediate objects directly, they are provided for reference.
\itemize{
\item $result: A per-sample raw result of scLinaX. This data has redundancy and should be summarized using the \code{summarize_scLinaX} function before analysis.

The dataframe ($result) has the following columns:
\itemize{
\item Sample_ID: Sample ID
\item SNP_ID: SNP identifier
\item CHR: PAR or nonPAR
\item POS: Genomic position of the SNP (GRCh38)
\item REF: Reference allele of the SNP {A,T,G,C}
\item ALT: Alternative allele of the SNP {A,T,G,C}
\item Gene: Gene name
\item XCI_status: XCI status {escape, variable, inactive, unknown}
\item Gene_class: XCI status combined with CHR information {PAR1, nonPAR_escape, nonPAR_variable, nonPAR_inactive, nonPAR_unknown, PAR2}
\item Used_as_refGene: Whether the Gene was used as a reference gene {Yes, No}
\item Used_as_refSNP: Whether the SNP_ID was used as a reference SNP {Yes, No}
\item Total_A_allele, Total_B_allele: Total allele count of allele A and B
\item Total_allele: Total allele count of the SNP
\item Expressing_cells: Number of cells expressing the SNP
\item minor_allele_ratio: Ratio of the expression from the allele (A, B) with lower expression
\item Reference_Gene, Reference_SNP: A list of genes used as reference genes and SNPs
\item Num_Reference_Gene, Num_Reference_SNP: Number of reference genes and SNPs
\item Reference_Cell_Count: Number of cells showing mono-allelic expression of reference SNPs
\item Num_A_cells, Num_B_cells: Number of cells showing mono-allelic expression of A and B alleles of reference SNPs
\item Num_Fail_cells: Number of cells showing bi-allelic expression of reference SNPs (should be removed from the analysis)
}
\item $raw_exp_result: A per-sample raw result of scLinaX.

The dataframe ($raw_exp_result) has the following columns:
\itemize{
\item cell_barcode: Cell barcode
\item Sample_ID: Sample ID
\item SNP_ID: SNP identifier
\item CHR: PAR or nonPAR
\item POS: Genomic position of the SNP (GRCh38)
\item REF: Reference allele of the SNP {A,T,G,C}
\item ALT: Alternative allele of the SNP {A,T,G,C}
\item REFcount: Allelic expression of the reference allele
\item ALTcount: Allelic expression of the alternative allele
\item OTHcount: Allelic expression of the other allele
\item Gene: Gene name
\item XCI_status: XCI status {escape, variable, inactive, unknown}
\item Gene_class: XCI status combined with CHR information {PAR1, nonPAR_escape, nonPAR_variable, nonPAR_inactive, nonPAR_unknown, PAR2}
\item Used_as_refGene: Whether the Gene was used as a reference gene {Yes, No}
\item Used_as_refSNP: Whether the SNP_ID was used as a reference SNP {Yes, No}
\item Xa: Information of the activated X chromosome {Allele_A, Allele_B}
\item Reference_Gene, Reference_SNP: A list of genes used as reference genes and SNPs
}
\item $Fail_list: List of samples for which scLinaX analysis failed.
\item $clustering_result: Result of the grouping of reference SNPs. Cluster names are in the format {Sample_ID}\emph{cluster}{cluster_ID}.
\item $Max_Num_Table_result: A dataframe describing the number of cells.
\item $df_snp_summary: An original dataframe from which Spearman correlation between pseudobulk ASE profiles is calculated.
\item $phasing_result: Result of the Spearman correlation analysis for pseudobulk ASE profiles.
}
}
\description{
Run scLinaX analysis for all samples
}

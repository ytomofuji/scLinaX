#' @importFrom dplyr filter select group_by summarize ungroup mutate left_join
#' @importFrom dplyr arrange distinct slice pull rename inner_join bind_rows
#' @importFrom tibble tibble
#' @importFrom stringr str_c
#' @importFrom magrittr %>%
#' @importFrom stats sd
NULL

utils::globalVariables(c("ALT", "ALT1", "ALT2", "ALTcount", "ALTratio", "A_allele", "A_allele_count",
                         "B_allele", "B_allele_count", "CHR", "Count", "Gene", "Gene_class",
                         "Major_allele", "Major_allele_count", "Mean_AR", "Mean_AR_reference",
                         "Mean_AR_target", "Minor_allele_count", "N_TotalCell", "OTHcount", "P",
                         "POS", "Peak_name", "REF", "REF1", "REF2", "REFcount", "Reference_Cell_Count",
                         "Reference_Gene", "Reference_SNP", "SNP", "SNP1", "SNP2", "SNP_ID", "Sample_ID",
                         "Sample_N", "TOTALcount", "Total_A_allele", "Total_B_allele", "Total_allele",
                         "Total_allele_count", "Used_as_refGene", "Used_as_refSNP", "XCI_annotation",
                         "XCI_ref", "XCI_status", "Xa", "case_when", "cell_barcode", "data", "flip_or_not",
                         "median", "minor_allele_ratio", "n", "n_distinct", "rho", "tag"))

#' Create a dataframe for QC of candidate reference genes
#' @inheritParams run_scLinaX
#' @inheritParams run_clustering

make_QC_gene_table<-function(QCed_df,per_SNP_sample,HE_allele_cell_number_THR=50,QC_total_allele_THR=10){
  HE_inactive<-dplyr::filter(per_SNP_sample,Count>=HE_allele_cell_number_THR)%>%
    dplyr::filter(Gene_class=="nonPAR_inactive")

  result<-tibble::tibble()

  for(i in seq(1,nrow(HE_inactive))){
    tGene<-HE_inactive$Gene[i]
    tSample_ID<-HE_inactive$Sample_ID[i]
    tSNP_ID<-HE_inactive$SNP_ID[i]
    tCell_Count<-HE_inactive$Count[i]
    tXCI_status<-HE_inactive$Gene_class[i]

    tmp_df<-dplyr::filter(QCed_df,Sample_ID==tSample_ID)
    REF_cells<-dplyr::filter(tmp_df,SNP_ID==tSNP_ID & REFcount>0 & ALTcount==0 & OTHcount==0)%>%dplyr::pull(cell_barcode)%>%unique()
    ALT_cells<-dplyr::filter(tmp_df,SNP_ID==tSNP_ID & REFcount==0 & ALTcount>0 & OTHcount==0)%>%dplyr::pull(cell_barcode)%>%unique()
    Other_cells<-dplyr::filter(tmp_df,SNP_ID==tSNP_ID & !(cell_barcode %in% c(REF_cells,ALT_cells)))%>%dplyr::pull(cell_barcode)%>%unique()

    REFexpdf<-dplyr::filter(tmp_df,cell_barcode %in% REF_cells)%>%
      dplyr::mutate(A_allele=REFcount,B_allele=ALTcount)
    ALTexpdf<-dplyr::filter(tmp_df,cell_barcode %in% ALT_cells)%>%
      dplyr::mutate(B_allele=REFcount,A_allele=ALTcount)
    expdf<-dplyr::bind_rows(REFexpdf,ALTexpdf)%>%dplyr::filter(SNP_ID != tSNP_ID)

    tmp_res<-expdf%>%
      dplyr::group_by(SNP_ID,CHR,POS,REF,ALT,Gene,XCI_status,Gene_class)%>%
      dplyr::summarize(Total_A_allele=sum(A_allele),Total_B_allele=sum(B_allele),Total_allele=(sum(A_allele)+sum(B_allele)))%>%
      dplyr::ungroup()%>%
      dplyr::mutate(minor_allele_ratio=ifelse(Total_A_allele>Total_B_allele,Total_B_allele/Total_allele,Total_A_allele/Total_allele))%>%
      dplyr::mutate(Sample_ID=tSample_ID,Reference_Gene=tGene,Reference_SNP=tSNP_ID,Reference_Cell_Count=tCell_Count,
             Num_REF_cells=length(REF_cells),Num_ALT_cells=length(ALT_cells),Num_Other_cells=length(Other_cells))

    result<-dplyr::bind_rows(result,tmp_res)
  }

  QCed_result<-dplyr::filter(result,Total_allele>=QC_total_allele_THR)
  ref_gene_list<-unique(HE_inactive$Gene)
  Ref_Gene_QC<-QCed_result%>%dplyr::filter(Gene %in% ref_gene_list)%>%
    dplyr::filter(Reference_Gene!=Gene) #remove cases where Reference/target gene was same
  return(Ref_Gene_QC)
}



#' Summarize the dataframe generated for QC
#' @param Ref_Gene_QC Output from the make_QC_gene_table function.
#' @param SAMPLE_NUM_THR Threshold for the sample size used in the calculation of the ratio of expression from Xi.
#' Genes evaluated in at least `SAMPLE_NUM_THR` samples are used for the calculation of the ratio of expression from Xi. Default: 3.

summarize_QC_gene_table<-function(Ref_Gene_QC,SAMPLE_NUM_THR=3){

  #QC-based on target gene
  plt_target_gene<-dplyr::group_by(Ref_Gene_QC,Gene,Reference_Gene,Sample_ID)%>% #collapse per gen
    dplyr::arrange(-Total_allele)%>%dplyr::slice(1)%>%dplyr::ungroup()%>%
    dplyr::group_by(Gene)%>%
    dplyr::summarize(Mean_AR=mean(minor_allele_ratio),SD_AR=sd(minor_allele_ratio),
              Mean_Total_allele=mean(Total_allele),SD_Total_allele=sd(Total_allele),Sample_N=n_distinct(Sample_ID),Count=n())%>%
    dplyr::ungroup()%>%
    dplyr::filter(Sample_N>=SAMPLE_NUM_THR)%>%
    dplyr::arrange(-Mean_AR)

  #QC-based on target gene
  plt_reference_gene<-dplyr::group_by(Ref_Gene_QC,Gene,Reference_Gene,Sample_ID)%>% #collapse per gene
    dplyr::arrange(-Total_allele)%>%slice(1)%>%ungroup()%>%
    dplyr::group_by(Reference_Gene)%>%
    dplyr::summarize(Mean_AR=mean(minor_allele_ratio),SD_AR=sd(minor_allele_ratio),
              Mean_Total_allele=mean(Total_allele),SD_Total_allele=sd(Total_allele),Sample_N=n_distinct(Sample_ID),Count=n())%>%
    dplyr::ungroup()%>%
    dplyr::filter(Sample_N>=SAMPLE_NUM_THR)%>%
    dplyr::arrange(-Mean_AR)

  #Comparison between reference/target-based QC
  colnames(plt_target_gene)[2:ncol(plt_target_gene)]<-str_c(colnames(plt_target_gene)[2:ncol(plt_target_gene)],"_target")
  colnames(plt_reference_gene)[2:ncol(plt_reference_gene)]<-str_c(colnames(plt_reference_gene)[2:ncol(plt_reference_gene)],"_reference")
  reference_gene_QC_summary<-dplyr::full_join(plt_target_gene,plt_reference_gene,by=c("Gene"="Reference_Gene"))

  return(reference_gene_QC_summary)
}





#' Run Gene Quality Control (QC) function
#' @inheritParams run_scLinaX
#' @inheritParams summarize_QC_gene_table
#' @export


run_RefGeneQC<-function(ASE_df,XCI_ref,SNP_DETECTION_DP=30,SNP_DETECTION_MAF=0.1,SAMPLE_NUM_THR=3,HE_allele_cell_number_THR=50,QC_total_allele_THR=10){
  QCed_df<-make_QCed_df(ASE_df,XCI_ref,SNP_DETECTION_DP,SNP_DETECTION_MAF)
  per_SNP_sample<-make_per_SNP_sample_df(QCed_df)
  QC_gene_table<-make_QC_gene_table(QCed_df,per_SNP_sample,HE_allele_cell_number_THR,QC_total_allele_THR)
  QC_result<-summarize_QC_gene_table(QC_gene_table,SAMPLE_NUM_THR)
  return(QC_result)
}






#' @importFrom dplyr filter bind_rows mutate left_join inner_join select group_by ungroup distinct slice pull rename arrange
#' @importFrom purrr map
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps IRanges
#' @importFrom stringr str_c
#' @importFrom magrittr %>%
#' @importFrom stats cor.test
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



#' Run scLinaX-multi using the scATAC data and output from scLinaX
#' @param scLinaX_summary Output from the summarize_scLinaX function run with Annotaion=NULL option
#' @param ATAC_ASE_data A dataframe (tibble) containing single-cell chromatin accessibility data. This dataframe should include the following columns:
#'   * SNP_ID: SNP identifier
#'   * POS: Genomic position of the SNP (GRCh38)
#'   * REF: Reference allele of the SNP {A,T,G,C}
#'   * ALT: Alternative allele of the SNP {A,T,G,C}
#'   * cell_barcode: Cell barcode
#'   * REFcount: Allelic read count of the reference allele
#'   * ALTcount: Allelic read count of the alternative allele
#'   * OTHcount: Allelic read count of the other allele
#'   * Sample_ID: Sample ID
#' @param peak_bed A dataframe (tibble) containing single-cell chromatin accessibility data. This dataframe should include the following columns:
#'   * start: Genomic position of the start of the peak (1-base)
#'   * end: Genomic position of the end of the peak (1-base)
#'   * Peak_name: Name of the peak
#'   * Gene: Nearest gene to the peak
#' @param SNP_DETECTION_DP Threshold for the total allele count (depth) of the SNP in the scATAC data. SNP–Sample pairs with a total allele count of at least "SNP_DETECTION_DP" are used for the analysis. Default: 30.
#' @inheritParams summarize_scLinaX
#' @inheritParams run_scLinaX
#' @export

run_scLinaX_multi<-function(scLinaX_summary,scLinaX_obj,ATAC_ASE_data,peak_bed,SNP_DETECTION_DP=30,SNP_DETECTION_MAF=0.1){

#prepare scATAC ASCA data
ref_snp<-unique(scLinaX_summary$Reference_SNP)
phase<-scLinaX_obj$raw_exp_result%>%
  filter(Reference_SNP==ref_snp)%>%
  dplyr::select(cell_barcode,Xa)%>%
  dplyr::distinct()

ase_data<-dplyr::inner_join(ATAC_ASE_data,phase,by="cell_barcode")

per_sample_SNP<-ase_data%>%
  dplyr::group_by(SNP_ID,Sample_ID)%>%
  dplyr::summarize(REFcount=sum(REFcount),ALTcount=sum(ALTcount))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(TOTALcount=REFcount+ALTcount)%>%
  dplyr::mutate(ALTratio=ALTcount/TOTALcount)

QCpassed_var_cell<-filter(per_sample_SNP,TOTALcount>=SNP_DETECTION_DP)%>%
  filter(ALTratio>SNP_DETECTION_MAF & ALTratio<(1-SNP_DETECTION_MAF))%>%
  dplyr::select(SNP_ID,Sample_ID)

# prepare peak data
annotated_peak_range<-GenomicRanges::GRanges(seqnames=rep('chrX',nrow(peak_bed)),
                              ranges=IRanges::IRanges(
                                start=peak_bed$start,end=peak_bed$end),
                              strand=rep('+',nrow(peak_bed)))
names(annotated_peak_range)<-peak_bed$Peak_name

# extract only SNPs overlapped with Peak
pos_list<-dplyr::select(ase_data,SNP_ID,POS)%>%dplyr::distinct()
ASE_GR<-GenomicRanges::GRanges(seqnames=rep('chrX',nrow(pos_list)),
                ranges=IRanges(
                  start=pos_list$POS,end=pos_list$POS),
                strand=rep('+',nrow(pos_list)))
names(ASE_GR)<-pos_list$SNP_ID

overlaps<-IRanges::findOverlaps(ASE_GR, annotated_peak_range)
query<-ASE_GR[overlaps@from,]%>%names()
subjects<-annotated_peak_range[overlaps@to,]%>%names()

peak_gene_info<-peak_bed%>%dplyr::select(Peak_name,Gene)
SNP_Peak_overlap<-tibble::tibble(SNP_ID=query,Peak_name=subjects)%>%left_join(peak_gene_info,by="Peak_name")

QCed_df<-dplyr::left_join(QCpassed_var_cell,ase_data,by=c("SNP_ID","Sample_ID"))%>%
  dplyr::select(SNP_ID,Sample_ID,POS,REF,ALT,cell_barcode,REFcount,ALTcount,OTHcount,Xa)%>%
  dplyr::inner_join(SNP_Peak_overlap,by="SNP_ID")%>%
  left_join(XCI_ref,by="Gene")
QCed_df$XCI_status[is.na(QCed_df$XCI_status)]<-"Unknown"

QCed_df<-QCed_df%>%mutate(CHR=case_when(POS<=2781479~"PAR1",POS>=155701383~"PAR2",TRUE~"nonPAR"))%>%
  mutate(Gene_class=case_when(
    CHR=="PAR1"~"PAR1",CHR=="PAR2"~"PAR2",
    XCI_status=="escape"~"nonPAR_escape",
    XCI_status=="variable"~"nonPAR_variable",
    XCI_status=="inactive"~"nonPAR_inactive",
    XCI_status=="Unknown"~"nonPAR_unknown"))%>%
  mutate(REFcount=ifelse(REFcount>1,1,REFcount),ALTcount=ifelse(ALTcount>1,1,ALTcount))

QCed_df$Gene_class<-factor(QCed_df$Gene_class,levels=c("PAR1","nonPAR_escape","nonPAR_variable","nonPAR_inactive","nonPAR_unknown","PAR2"))


phased_df<-QCed_df%>%
  mutate(A_allele_count=ifelse(Xa=="Allele_A",REFcount,ALTcount),B_allele_count=ifelse(Xa=="Allele_A",ALTcount,REFcount))

return(phased_df)
}





#' Summarize scLinaX-multi results
#' @param phased_df A dataframe generated by the `run_scLinaX_multi` function.
#' @inheritParams summarize_scLinaX
#' @export

summarize_scLinaX_multi<-function(phased_df,QC_total_allele_THR=10,Annotation=NULL){

  summary<-phased_df%>%
    dplyr::group_by(SNP_ID,Peak_name,Sample_ID,POS,REF,ALT,Gene,Gene_class)%>%
    dplyr::summarize(A_allele_count=sum(A_allele_count),B_allele_count=sum(B_allele_count))%>%
    dplyr::mutate(Total_allele_count=A_allele_count+B_allele_count)%>%
    dplyr::mutate(minor_allele_ratio=ifelse(A_allele_count>B_allele_count,B_allele_count/Total_allele_count,A_allele_count/Total_allele_count),
           Major_allele=ifelse(A_allele_count>B_allele_count,"A","B"))%>%
    dplyr::ungroup()%>%
    dplyr::filter(Total_allele_count>=QC_total_allele_THR)%>%
    dplyr::group_by(Peak_name)%>%
    dplyr::arrange(-Total_allele_count)%>%slice(1)%>%ungroup()

  if(is.null(Annotation)){
    return(summary)
  }else{
    major_allele_info<-dplyr::select(summary,SNP_ID,Major_allele)

    per_anno_summary<-phased_df%>%
      dplyr::inner_join(Annotation,by="cell_barcode")%>%
      dplyr::left_join(major_allele_info,by="SNP_ID")%>%
      dplyr::filter(SNP_ID %in% summary$SNP_ID)%>%
      dplyr::group_by(SNP_ID,Sample_ID,Annotation,POS,REF,ALT,Gene,Gene_class,Major_allele)%>%
      dplyr::summarize(A_allele_count=sum(A_allele_count),B_allele_count=sum(B_allele_count))%>%
      dplyr::mutate(Total_allele_count=A_allele_count+B_allele_count)%>%
      dplyr::mutate(minor_allele_ratio=ifelse(Major_allele=="A",B_allele_count/Total_allele_count,A_allele_count/Total_allele_count))%>%
      dplyr::ungroup()%>%
      dplyr::filter(Total_allele_count>=QC_total_allele_THR)

    return(per_anno_summary)
  }
}




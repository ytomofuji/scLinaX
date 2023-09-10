#prepare test data which will be stored in the scLinaX R package

library(tidyverse)
data<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/test_mutiome_RNA_scASE.tsv.gz")
data$cell_barcode<-gsub("-RNA","-1",data$cell_barcode)

metadata<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/Multiome.rm.plt.ery.db.metadata.tsv.gz")

multiome_ASE_df<-filter(data,cell_barcode %in% metadata$cell_barcode)

XCI_ref<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/scLinaX_XCI_status_ref.tsv")

AIDA_QCREF<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/test_reference_Gene_QC_AIDA.tsv")

multiome_Annotation<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/Multiome.rm.plt.ery.db.metadata.tsv.gz")%>%
  dplyr::select(cell_barcode,predicted.id)
colnames(multiome_Annotation)<-c("cell_barcode","Annotation")

multiome_ATAC_ASE_data<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/Multiome_ATAC_QC_passed_SNP_df.tsv.gz")%>%
  mutate(Sample_ID="Multiome")
multiome_ATAC_ASE_data$cell_barcode<-gsub("-ATAC","-1",multiome_ATAC_ASE_data$cell_barcode)
multiome_ATAC_ASE_data<-filter(multiome_ATAC_ASE_data,cell_barcode %in% metadata$cell_barcode)
multiome_ATAC_ASE_data<-dplyr::select(multiome_ATAC_ASE_data,-Region,-Gene,-OKG_All_Freq,-CHR)


multiome_peak_data<-read_tsv("C:/Users/YoshihikoTomofuji/Desktop/Analysis/scLinaX/Original_Data/Multiome_peak_data.tsv.gz")%>%dplyr::select(-seqname)

usethis::use_data(multiome_ASE_df)
usethis::use_data(XCI_ref)
usethis::use_data(AIDA_QCREF)
usethis::use_data(multiome_Annotation,overwrite = TRUE)
usethis::use_data(multiome_ATAC_ASE_data)
usethis::use_data(multiome_peak_data)

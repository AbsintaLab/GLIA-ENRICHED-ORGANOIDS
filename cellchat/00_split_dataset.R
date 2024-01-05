# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
test_all <- readRDS("../../../BS_drop_RR25CTRL/data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
DimPlot(test_all)

# add the splitting variable
table(test_all$orig.ident)
test_all$GroupCCC <- case_when(test_all$orig.ident %in% c("CTRL4_SOX10_CTRL","CTRL8_SOX10_CTRL","RR16_SOX10_CTRL")~"GroupCCC_CTRL",
          test_all$orig.ident %in% c("CTRL8_SOX10_CSF","RR16_SOX10_CSF","RR25_SOX10_CSF")~"GroupCCC_CSF",
          T~"exclude")

table(test_all$orig.ident,test_all$GroupCCC)

# split he object into individula one besed on the pathology variable
list <- SplitObject(test_all, split.by = "GroupCCC")

# save each object
pmap(list(list,names(list)),function(x,y){
  saveRDS(x,file = paste0("../../out/object/subset_",y,"_data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds"))
})

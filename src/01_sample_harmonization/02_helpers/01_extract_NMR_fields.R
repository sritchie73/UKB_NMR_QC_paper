library(data.table)
library(ukbnmr)

fields <- data.table(field=c(
  na.omit(ukbnmr::nmr_info$UKB.Field.ID),
  na.omit(ukbnmr::nmr_info$QC.Flag.Field.ID),
  ukbnmr:::sample_qc_fields$UKB.Field.ID
))

fwrite(fields, quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/nmr_fields.txt")

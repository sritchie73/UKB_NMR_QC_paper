library(data.table)
library(ukbnmr) # remotes::install_github("sritchie73/ukbnmr", ref="development")
library(MASS)

# Create output directories
if (!dir.exists("data/tech_qc")) dir.create("data/tech_qc")

################################################
# Load raw data that has undergone sample QC
################################################

raw <- fread("data/raw/processed/UKB_phase1_results.txt")
sinfo <- fread("data/raw/processed/sample_information.txt")

# melt to long
raw <- melt(raw, id.vars=c("sample_id", "plate_id", "plate_position", "visit", "spectrometer"))
raw <- raw[!is.na(value)] # simplifies downstream code to drop NAs now

################################################
# Filter to non-derived biomarkers
################################################

# Filter raw data to non-derived biomarkers
raw <- raw[variable %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]

################################################
# Log transform
################################################

# Determine offset for log transformation (required for variables with measurements == 0):
# Acetate, Acetoacetate, Albumin, bOHbutyrate, Clinical_LDL_C, Gly, Ile, L_LDL_CE, L_LDL_FC, M_LDL_CE
log_offset <- raw[!is.na(value),.(min=min(value), non_zero=min(value[value != 0])),by=variable]
log_offset[, log_offset := ifelse(min == 0, non_zero / 2, 0)]

# Get log transformed raw value
raw[log_offset, on = .(variable), log_value := log(value + log_offset)]

################################################
# Merge in relevant technical covariates
################################################

raw[sinfo[!(removed)], on = .(sample_id, plate_id, plate_position, visit, spectrometer), 
    c("sample_degredation", "plate_row", "plate_column") := .(prep_to_measured, well_row, well_column)]

#################################################################
# Step 1 - Adjust for sample degredation
#################################################################

raw[, log_adj1 := rlm(log_value ~ sample_degredation)$residuals, by=.(variable)]

#################################################################
# Step 2 - adjust for within plate structure across rows A-H
#################################################################

# Function for converting a set of numbers to a factor, ordering levels by group size
factor_by_size <- function(x) {
  group_sizes <- as.data.table(table(x))
  group_sizes <- group_sizes[order(-N)]
  factor(x, levels=group_sizes$x)
}

raw[, log_adj2 := rlm(log_adj1 ~ factor_by_size(plate_row))$residuals, by=.(variable)]

#################################################################
# Step 3 - adjust for within plate structure across rows A-H
#################################################################

raw[, log_adj3 := rlm(log_adj2 ~ factor_by_size(plate_column))$residuals, by=.(variable)]

#################################################################
# Step 4 - Adjust for drift over time within each spectrometer
#################################################################

# Function to cut vector x into N bins of equal duration
bin_duration <- function(x, n=10) {
  offset_x <- x + 1 # offset by 1 so we can apply ceiling function
  bin_size <- max(offset_x)/n
  as.integer(ceiling(offset_x/bin_size))
}

# Convert plate measure dates into a sequence of 1:N distinct dates per spectrometer, then cut into 10 equal size bins
sinfo[, plate_measure_bin := bin_duration(as.integer(factor(plate_MEASURED_DAYS))), by=spectrometer]
raw[sinfo, on = .(sample_id, plate_id, plate_position, visit, spectrometer), plate_measure_bin := plate_measure_bin]

# Adjust
raw[, log_adj4 := rlm(log_adj3 ~ factor_by_size(plate_measure_bin))$residuals, by=.(variable, spectrometer)]

################################################
# Rescale to absolute units
################################################

# Rescale to absolute units
rescale <- function(log_adj, log_raw) {
  # The residuals from any regression are defined as the difference between the
  # observed independent variable (e.g. biomarker concentrations) and the parameter
  # estimated by the regression. A resulting key property is that the distribution
  # of the residuals is centred on 0 for this estimated parameter. In the case of
  # robust linear regression, this parameter is an estimate of the mean that is
  # robust to outliers. As a consequence of the way residuals are defined, their
  # distribution can be scaled to match the distribution of the independent variable
  # by giving it the same estimated mean. For robust linear regression, residuals
  # can be returned to the same scale as the observed independent variable by
  # estimating the mean of the observed independent variable using robust linear
  # regression and adding it to the residuals.
  adj <- log_adj + as.vector(coef(rlm(log_raw ~ 1)))[1]

  # Now we can reverse the log transformation
  adj <- exp(adj) 

  return(adj)
}

# Put values back on absolute units scale
raw[, adj1 := rescale(log_adj1, log_value), by=.(variable)]
raw[, adj2 := rescale(log_adj2, log_value), by=.(variable)]
raw[, adj3 := rescale(log_adj3, log_value), by=.(variable)]
raw[, adj4 := rescale(log_adj4, log_value), by=.(variable)]

# Remove the small offset applied for biomarkers with 0 values
raw[log_offset, on = .(variable), adj1 := adj1 - log_offset]
raw[log_offset, on = .(variable), adj2 := adj2 - log_offset]
raw[log_offset, on = .(variable), adj3 := adj3 - log_offset]
raw[log_offset, on = .(variable), adj4 := adj4 - log_offset]

# Some values that were 0 are now < 0, apply small right shift
# for these biomarkers (shift is very small, i.e. the impact on 
# the distribution's median is essentially numeric error)
shift <- raw[,.(adj1=-pmin(0, min(adj1)),
                adj2=-pmin(0, min(adj2)),
                adj3=-pmin(0, min(adj3)),
                adj4=-pmin(0, min(adj4))
              ), by=variable]

# Shift
raw[shift, on = .(variable), adj1 := adj1 + i.adj1]
raw[shift, on = .(variable), adj2 := adj2 + i.adj2]
raw[shift, on = .(variable), adj3 := adj3 + i.adj3]
raw[shift, on = .(variable), adj4 := adj4 + i.adj4]

# Output log offset information
log_offset[shift, on = .(variable), right_shift := adj4]
fwrite(log_offset, sep="\t", file="data/tech_qc/log_offset.txt")

################################################
# Identify and remove outlier plates
################################################

# Model plate medians as a normal distribution (across all 1,352 plates), then
# we can flag outlier plates as those > 3.374446 standard deviations from the
# mean where 3.374446 is derived from the range of the theoretical normal
# distribution for 1,352 samples
n_plates <- sinfo[!(removed), length(unique(plate_id))] # 1,352
sdlim <- max(qnorm(ppoints(n_plates))) # 3.374446
plate_medians <- raw[!is.na(value),.(value = median(adj4)), by=.(variable, plate_id)]
outlier_lim <- plate_medians[,.(lower_lim = mean(value) - sd(value) * sdlim, upper_lim = mean(value) + sd(value) * sdlim), by=variable]
plate_medians[, outlier := "no"]
plate_medians[outlier_lim, on = .(variable, value < lower_lim), outlier := "low"]
plate_medians[outlier_lim, on = .(variable, value > upper_lim), outlier := "high"]

# Remove outlier plates
raw[, adj5 := adj4]
raw[plate_medians[outlier != "no"], on = .(variable, plate_id), adj5 := NA]

# Write out info
fwrite(plate_medians, sep="\t", quote=FALSE, file="data/tech_qc/plate_medians_outlier_tagged.txt")
fwrite(outlier_lim, sep="\t", quote=FALSE, file="data/tech_qc/plate_outlier_limits.txt")

################################################
# Compute composite and derived biomarkers
################################################

# Cast adjusted values to wide
adj1 <- dcast(raw, sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj1", fill=NA)
adj2 <- dcast(raw, sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj2", fill=NA)
adj3 <- dcast(raw, sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj3", fill=NA)
adj4 <- dcast(raw, sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj4", fill=NA)
adj5 <- dcast(raw, sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj5", fill=NA)

# Recompute composite and derived biomarkers
adj1 <- recompute_derived_biomarkers(adj1)
adj2 <- recompute_derived_biomarkers(adj2)
adj3 <- recompute_derived_biomarkers(adj3)
adj4 <- recompute_derived_biomarkers(adj4)
adj5 <- recompute_derived_biomarkers(adj5)

# Melt to long
adj1 <- melt(adj1, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="adj1")
adj2 <- melt(adj2, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="adj2")
adj3 <- melt(adj3, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="adj3")
adj4 <- melt(adj4, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="adj4")
adj5 <- melt(adj5, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="adj5")

# Combine into single data.table
adj <- merge(adj1, adj2, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"))
adj <- merge(adj, adj3, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"))
adj <- merge(adj, adj4, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"))
adj <- merge(adj, adj5, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"))

# Add in values prior to rescaling:
log_adj <- raw[,.(sample_id, visit, plate_id, plate_position, spectrometer, variable=as.vector(variable), log_adj1, log_adj2, log_adj3, log_adj4)]
adj <- merge(log_adj, adj, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"), all.y=TRUE)

# Also get version of raw data with rederived biomarkers
raw_rederived <- fread("data/raw/processed/UKB_phase1_results.txt")
raw_rederived <- recompute_derived_biomarkers(raw_rederived)
raw_rederived <- melt(raw_rederived, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="raw_rederived", variable.factor=FALSE)
adj <- merge(raw_rederived, adj, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"), all.y=TRUE)

# And log transformed raw data for non-derived biomarkers
raw <- raw[,.(sample_id, visit, plate_id, plate_position, spectrometer, variable=as.vector(variable), log_raw=log_value)]
adj <- merge(raw, adj, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"), all.y=TRUE)

# Reload raw data so we can add in those values (current 'raw' data.table only has non-derived biomarkers)
raw <- fread("data/raw/processed/UKB_phase1_results.txt")
raw <- melt(raw, id.vars=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer"), value.name="raw", variable.factor=FALSE)
adj <- merge(raw, adj, by=c("sample_id", "visit", "plate_id", "plate_position", "spectrometer", "variable"), all.y=TRUE)

# Add back in relevant technical covariates for plotting
tech_covar <- sinfo[!(removed), .(sample_id, plate_id, plate_position, visit, spectrometer, 
                    sample_degredation=prep_to_measured, plate_row=well_row, plate_column=well_column, plate_MEASURED_DAYS, plate_measure_bin)]
adj <- tech_covar[adj, on = .(sample_id, plate_id, plate_position, visit, spectrometer)]

# Link to participant IDs
ukblink <- fread("data/raw/nightingale/Cambridge_Nightingale_bridge.tsv")
adj <- merge(ukblink, adj, by=c("sample_id", "plate_id", "plate_position"), all.y=TRUE)

# Write out version with all steps of adjustment
fwrite(adj, sep="\t", quote=FALSE, file="data/tech_qc/multistep_adjusted_values.txt")

##########################################################
# Write out post-qc dataset(s)
##########################################################

# Write out post-QC data including outlier plates
adj4 <- dcast(adj, eid_7439 + sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj4")
fwrite(adj4, sep="\t", quote=FALSE, file="data/tech_qc/postqc_with_outlier_plates.txt")

# Write out post-QC data excluding outlier plates
adj5 <- dcast(adj, eid_7439 + sample_id + visit + plate_id + plate_position + spectrometer ~ variable, value.var="adj5")
fwrite(adj5, sep="\t", quote=FALSE, file="data/tech_qc/postqc_excluding_outlier_plates.txt")

# Add participant IDs to sinfo
sinfo[ukblink, on = .(sample_id, plate_id, plate_position), eid_7439 := eid_7439]

# Collate QC Flags for each biomarker:
# - Need to add flags for outlier plates
# - For derived biomarkers, collate flags from parts
flags <- fread("data/raw/processed/biomarker_QC_tags.txt")
outlier_plates <- plate_medians[outlier != "no"]
outlier_plates[, QC_tag := sprintf("%s_outlier_plate", outlier)]
outlier_plates <- merge(sinfo[,.(sample_id, visit, plate_id)], outlier_plates, by="plate_id", all.y=TRUE)
outlier_plates <- outlier_plates[, .(sample_id, visit, biomarker=variable, QC_tag)]
flags <- rbind(outlier_plates, flags)

flags[, QC_tag := gsub("_", " ", QC_tag)] # Harmonize with flags in UKB
flags[QC_tag == "ethanol", QC_tag := "Ethanol"]
flags[QC_tag == "high outlier plate", QC_tag := "High outlier plate"]
flags[QC_tag == "low outlier plate", QC_tag := "Low outlier plate"]

flags <- flags[, .(QC_tag = paste(sort(QC_tag), collapse="; ")), by=.(sample_id, visit, biomarker)]
flags[sinfo, on = .(sample_id, visit), eid_7439 := eid_7439]

flags <- flags[biomarker %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]
flags <- dcast(flags, eid_7439 + sample_id + visit ~ biomarker, value.var="QC_tag")
flags <- ukbnmr::recompute_derived_biomarker_qc_flags(flags)
fwrite(flags, sep="\t", quote=FALSE, file="data/tech_qc/biomarker_QC_flags.txt")

# Collate harmonized sample information. Keep:
# - standard sample identifier information (sample_id, plate_id, plate_position, spectrometer)
# - row number in original raw data (as the above is not always unique)
# - linked participant ID on project 7439
# - additional covariates used to adjust for technical variation (derived variables)
# - Additional sample flags, e.g. whether the sample was removed during our QC (and why),
#   sample QC flags, whether the sample belonged to a blind duplicate (eid_7439 visit combination),
#   and whether the sample is present in the released UKB raw data
sinfo <- fread("data/raw/processed/harmonized_sample_information.txt")
sinfo[, plate_measured_bin := bin_duration(as.integer(factor(plate_MEASURED_DAYS))), by=spectrometer]
sinfo <- sinfo[,.(raw_row_id=uid, eid_7439, Gender, visit, sample_id, plate_id, plate_position, spectrometer, 
                  plate_row=well_row, plate_column=well_column,
                  sample_degredation_time=prep_to_measured, plate_MEASURED_DATE, plate_measured_bin,
                  QC_flag=tag, sample_removed=removed, removal_reason,
                  blind_duplicate, in_ukb_raw=ifelse(in_ukb_raw == "No biomarkers", FALSE, TRUE))]
sinfo[, QC_flag := gsub(";", "; ", gsub("_", " ", QC_flag))]
fwrite(sinfo, sep="\t", quote=FALSE, file="data/tech_qc/sample_information.txt")


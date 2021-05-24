library(data.table)
library(ukbnmr) # remotes::install_github("sritchie73/ukbnmr", ref="development")
library(MASS)
library(ggplot2)
library(ggthemes)
library(palettetown)

## Warning: generating all the diagnostic plots takes ~10 hours.

# Create output directories
if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")
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

################################################
# Diagnostic plots
################################################

if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")

# Make sure biomarkers always shown in alphabetical order
adj <- adj[order(tolower(variable))]
adj[, variable := factor(variable, levels=unique(variable))]

# Analytical validity of biomarker rederivation
gg_dt <- adj[variable %in% ukbnmr::nmr_info[Type != "Non-derived" & (Nightingale), Biomarker]]
gg_dt <- gg_dt[!is.na(raw) & !is.na(raw_rederived)]

g <- ggplot(gg_dt) +
  aes(x=raw, y=raw_rederived) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_y_continuous(name="Derived and composite biomarkers computed from raw parts") +
  scale_x_continuous(name="Biomarker in raw data") +
  facet_wrap(~ variable, ncol=12, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=16.8, height=16.8, units="in", file="diagnostic_plots/rederivation_validity.png")

g <- ggplot(gg_dt) +
  aes(x=raw, y=raw - raw_rederived) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_hline(yintercept=0, linetype=2, color="red", size=0.4) +
  scale_y_continuous(name="Difference between raw and derived") +
  scale_x_continuous(name="Biomarker in raw data") +
  facet_wrap(~ variable, ncol=12, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=16.8, height=16.8, units="in", file="diagnostic_plots/rederivation_validity_diff.png")

# Rescale residuals to match log raw data for comparison
rescale <- function(log_adj, log_raw) {
  log_adj + as.vector(coef(rlm(log_raw ~ 1)))[1]
}

# Put values back on absolute units scale
adj[!is.na(log_raw), log_adj1 := rescale(log_adj1, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj2 := rescale(log_adj2, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj3 := rescale(log_adj3, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj4 := rescale(log_adj4, log_raw), by=.(variable)]

# Need to create version for adj5 (dropping outlier plates)
adj[!is.na(log_adj4) & !is.na(adj5), log_adj5 := log_adj4]

# Compare concentrations before and after adjustment
conc_comp <- function(x_col, y_col, xlab, ylab, fname) {
  setnames(adj, c(x_col, y_col), c("x_col", "y_col"))
  on.exit(setnames(adj, c("x_col", "y_col"), c(x_col, y_col)))
  n_bio <- adj[!is.na(x_col) & !is.na(y_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  g <- ggplot(adj[!is.na(x_col) & !is.na(y_col)], aes(x = x_col, y = y_col)) + 
    geom_point(shape=21, size=0.45, stroke=0.2, color="black", fill="white") +
    geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_x_continuous(name=xlab) +
    scale_y_continuous(name=ylab) +
		theme_bw() +
		theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
					strip.background=element_blank(), strip.text=element_text(size=6),
					panel.grid.major=element_blank(), panel.grid.minor=element_blank()
    )

  ggsave(g, width=width, height=width, units="in", file=fname)
}

conc_comp("raw", "adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_concentrations.png")
conc_comp("raw", "adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_concentrations.png")
conc_comp("raw", "adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_concentrations.png")
conc_comp("raw", "adj4", "Raw biomarker concentrations", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/raw_vs_adj4_concentrations.png")
conc_comp("raw", "adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_concentrations.png")
conc_comp("adj1", "adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_concentrations.png")
conc_comp("adj2", "adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_concentrations.png")
conc_comp("adj3", "adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_concentrations.png")
conc_comp("adj4", "adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_concentrations.png")

conc_comp("log_raw", "log_adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_log_concentrations.png")
conc_comp("log_raw", "log_adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_log_concentrations.png")
conc_comp("log_raw", "log_adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_log_concentrations.png")
conc_comp("log_raw", "log_adj4", "Raw biomarker concentrations", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/raw_vs_adj4_log_concentrations.png")
conc_comp("log_raw", "log_adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_log_concentrations.png")
conc_comp("log_adj1", "log_adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_log_concentrations.png")
conc_comp("log_adj2", "log_adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_log_concentrations.png")
conc_comp("log_adj3", "log_adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_log_concentrations.png")
conc_comp("log_adj4", "log_adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_log_concentrations.png")

# Compare distributions before and after adjustment
dens_comp <- function(raw_col, adj_col, raw_lab, adj_lab, fname) {
  gg_dt <- melt(adj, id.vars=c("sample_id", "visit", "variable"), measure.vars=c(raw_col, adj_col), variable.name="type")
  gg_dt <- gg_dt[!is.na(value)]
  n_bio <- gg_dt[, length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  g <- ggplot(gg_dt, aes(x = value, color=type)) + 
    geom_density(trim=TRUE, size=0.4) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_manual(name="", labels=structure(c(raw_lab, adj_lab), names=c(raw_col, adj_col)), 
                        values=structure(c("#e08214", "#542788"), names=c(raw_col, adj_col))) +
    xlab("Biomarker concentrations") +
    ylab("Density") +
		theme_bw() +
		theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
					strip.background=element_blank(), strip.text=element_text(size=6),
					panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    )

  ggsave(g, width=width, height=width, units="in", file=fname)
}

dens_comp("raw", "adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_densities.png")
dens_comp("raw", "adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_densities.png")
dens_comp("raw", "adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_densities.png")
dens_comp("raw", "adj4", "Raw biomarker concentrations", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/raw_vs_adj4_densities.png")
dens_comp("raw", "adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_densities.png")
dens_comp("adj1", "adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_densities.png")
dens_comp("adj2", "adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_densities.png")
dens_comp("adj3", "adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_densities.png")
dens_comp("adj4", "adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_densities.png")

dens_comp("log_raw", "log_adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_log_densities.png")
dens_comp("log_raw", "log_adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_log_densities.png")
dens_comp("log_raw", "log_adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_log_densities.png")
dens_comp("log_raw", "log_adj4", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj4_log_densities.png")
dens_comp("log_raw", "log_adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_log_densities.png")
dens_comp("log_adj1", "log_adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_log_densities.png")
dens_comp("log_adj2", "log_adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_log_densities.png")
dens_comp("log_adj3", "log_adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_log_densities.png")
dens_comp("log_adj4", "log_adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_log_densities.png")

# How does sample degradation time change at each step?
sample_deg <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  g <- ggplot(adj[!is.na(conc_col) & is.finite(conc_col)], aes(x=sample_degredation, y=conc_col)) +
    geom_hex() +
    scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10") +
    geom_smooth(color="red", size=0.2, method=MASS::rlm) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(conc_name) + 
    xlab("Sample degradation time (hours from sample prep to measurement") +
    theme_bw() +
		theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
					strip.background=element_blank(), strip.text=element_text(size=6),
					panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom"
    )

  ggsave(g, width=width, height=width, units="in", file=fname)
}

sample_deg("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_sample_degradation_time.png")
sample_deg("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_sample_degradation_time.png")
sample_deg("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_sample_degradation_time.png")
sample_deg("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_sample_degradation_time.png")
sample_deg("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_sample_degradation_time.png")
sample_deg("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_sample_degradation_time.png")

sample_deg("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_sample_degradation_time.png")
sample_deg("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_sample_degradation_time.png")
sample_deg("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_sample_degradation_time.png")
sample_deg("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_sample_degradation_time.png")
sample_deg("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_sample_degradation_time.png")
sample_deg("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_sample_degradation_time.png")

# Show concentrations vs. plate row
plate_row <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by row)
  well_summary <- well_summary[order(plate_column)][order(plate_row)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_tableau(name="Row") +
    scale_fill_tableau(name="Row") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by column (1-12) within each row (A-H))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_row("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_row.png")
plate_row("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_row.png")
plate_row("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_row.png")
plate_row("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_row.png")
plate_row("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_row.png")
plate_row("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_row.png")

plate_row("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_row.png")
plate_row("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_row.png")
plate_row("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_row.png")
plate_row("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_row.png")
plate_row("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_row.png")
plate_row("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_row.png")

# Zoom in on interquartile range
plate_row_iqr <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by row)
  well_summary <- well_summary[order(plate_column)][order(plate_row)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_tableau(name="Row") +
    scale_fill_tableau(name="Row") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by column (1-12) within each row (A-H))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_row_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_row_iqr_zoom.png")

plate_row_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_row_iqr_zoom.png")

# Show concentrations vs. plate column
plate_column <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by column)
  well_summary <- well_summary[order(plate_row)][order(plate_column)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_poke(name="Column", pokemon="typhlosion", spread=12) +
    scale_fill_poke(name="Column", pokemon="typhlosion", spread=12) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by row (A-H) within each column (1-12))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_column("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_column.png")
plate_column("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_column.png")
plate_column("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_column.png")
plate_column("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_column.png")
plate_column("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_column.png")
plate_column("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_column.png")

plate_column("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_column.png")
plate_column("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_column.png")
plate_column("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_column.png")
plate_column("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_column.png")
plate_column("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_column.png")
plate_column("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_column.png")

# Zoom in on interquartile range
plate_column_iqr <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by column)
  well_summary <- well_summary[order(plate_row)][order(plate_column)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_poke(name="Column", pokemon="typhlosion", spread=12) +
    scale_fill_poke(name="Column", pokemon="typhlosion", spread=12) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by row (A-H) within each column (1-12))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_column_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_column_iqr_zoom.png")

plate_column_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_column_iqr_zoom.png")

# Show concentrations vs. inter-plate variation (drift over time) within each spectrometer
spec_drift <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date within spectrometer
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)][order(spectrometer)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order within spectrometer)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

spec_drift("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_spectrometer_drift.png")
spec_drift("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_spectrometer_drift.png")
spec_drift("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_spectrometer_drift.png")
spec_drift("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_spectrometer_drift.png")
spec_drift("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_spectrometer_drift.png")
spec_drift("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_spectrometer_drift.png")

spec_drift("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_spectrometer_drift.png")
spec_drift("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_spectrometer_drift.png")
spec_drift("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_spectrometer_drift.png")
spec_drift("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_spectrometer_drift.png")
spec_drift("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_spectrometer_drift.png")
spec_drift("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_spectrometer_drift.png")

# Zoom in on interquartile range
spec_drift_iqr <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date within spectrometer
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)][order(spectrometer)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order within spectrometer)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

spec_drift_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_spectrometer_drift_iqr_zoom.png")

spec_drift_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_spectrometer_drift_iqr_zoom.png")

# Show concentrations vs. inter-plate variation without ordering by spectrometer
plate_by_date <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}


plate_by_date("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_by_date.png")
plate_by_date("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_by_date.png")
plate_by_date("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_by_date.png")
plate_by_date("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_by_date.png")
plate_by_date("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_by_date.png")
plate_by_date("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_by_date.png")

plate_by_date("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_by_date.png")
plate_by_date("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_by_date.png")
plate_by_date("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_by_date.png")
plate_by_date("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_by_date.png")
plate_by_date("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_by_date.png")
plate_by_date("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_by_date.png")

# Zoom in on interquartile range
plate_by_date_iqr <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_by_date_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_by_date_iqr_zoom.png")

plate_by_date_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_by_date_iqr_zoom.png")

# Show concentrations vs. inter-plate variation without ordering by spectrometer
outlier_plates <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))
  
  n_bio <- adj[!is.na(conc_col) & !is.na(log_raw), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4
  
  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col) & !is.na(log_raw),
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))),
                       by=.(plate_id, plate_MEASURED_DAYS, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  
  # Order plates by measurement date 
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]
  
  # Add outlier info
  plate_summary[plate_medians, on = .(variable, plate_id), outlier := i.outlier]
  
  # Split into three color groups (ggnewscale was not behaving)
  gg_no_outlier <- plate_summary[outlier == "no"]
  gg_high_outlier <- plate_summary[outlier == "high"]
  gg_low_outlier <- plate_summary[outlier == "low"]
  
  g <- ggplot(plate_summary, aes(x=factor(plate_order))) +
    geom_point(data=gg_no_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252") +
    geom_point(data=gg_no_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252") +
    geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252") +
    geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252") +

    geom_point(data=gg_high_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15") +
    geom_point(data=gg_high_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15") +
    geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15") +
    geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15") +
    
    geom_point(data=gg_low_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c") +
    geom_point(data=gg_low_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c") +
    geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c") +
    geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c") +
    
    geom_hline(data=outlier_lim, aes(yintercept=lower_lim), color="#dd1c77", linetype=2, size=0.2) +
    geom_hline(data=outlier_lim, aes(yintercept=upper_lim), color="#dd1c77", linetype=2, size=0.2) +

    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) +
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))
  
  ggsave(g, width=width, height=width, units="in", file=fname)
}

outlier_plates("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_outlier_plates.png")

# Zoom in on interquartile range
outlier_plates_iqr <- function(conc_col, conc_name, fname) {
  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))
  
  n_bio <- adj[!is.na(conc_col) & !is.na(log_raw), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4
  
  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col) & !is.na(log_raw),
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))),
                       by=.(plate_id, plate_MEASURED_DAYS, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  
  # Order plates by measurement date 
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]
  
  # Add outlier info
  plate_summary[plate_medians, on = .(variable, plate_id), outlier := i.outlier]
  
  # Split into three color groups (ggnewscale was not behaving)
  gg_no_outlier <- plate_summary[outlier == "no"]
  gg_high_outlier <- plate_summary[outlier == "high"]
  gg_low_outlier <- plate_summary[outlier == "low"]
  
  g <- ggplot(plate_summary, aes(x=factor(plate_order))) +
    geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252") +
    geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252") +
    
    geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15") +
    geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15") +
    
    geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c") +
    geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c") +
    
    geom_hline(data=outlier_lim, aes(yintercept=lower_lim), color="#dd1c77", linetype=2, size=0.2) +
    geom_hline(data=outlier_lim, aes(yintercept=upper_lim), color="#dd1c77", linetype=2, size=0.2) +
    
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) +
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))
  
  ggsave(g, width=width, height=width, units="in", file=fname)
}

outlier_plates_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_outlier_plates_iqr_zoom.png")

# Now show also on log scale
outlier_lim[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj4_2 := log(adj4 - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj4_2 ~ 1))[1]) - as.vector(coef(glm(log_adj4 ~ 1))[1])), by=variable]
outlier_lim[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_outlier_plates.png")
outlier_plates_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_outlier_plates_iqr_zoom.png")

# Show internal controls
cntrl <- fread("data/raw/processed/NGH_control_samples.txt")
cntrl_info <- fread("data/raw/processed/controls_information.txt")

cntrl <- melt(cntrl, id.vars=c("sample_id", "plate_id", "plate_position", "spectrometer"))
cntrl <- cntrl[variable %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]
cntrl[cntrl_info, on = .(sample_id, plate_id, plate_position, spectrometer), plate_MEASURED_DAYS := plate_MEASURED_DAYS]
cntrl <- cntrl[order(tolower(variable))]
cntrl[, variable := factor(variable, levels=unique(variable))]
cntrl[plate_medians, on = .(variable, plate_id), outlier := i.outlier] 
cntrl[, control_sample := sprintf("%s (Well %s)", sample_id, plate_position)]
cntrl[, control_sample := factor(control_sample, levels=c( # Control samples are paired
        "190404 (Well A01)", "190328 (Well H12)",
        "190425 (Well A01)", "190508 (Well H12)",
        "180829 (Well A01)", "180827 (Well H12)",
        "180830 (Well A01)", "180831 (Well H12)"
      ))]

internal_controls <- function(fname) {
  n_bio <- cntrl[, length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4
  
  # Order plates by measurement date 
  cntrl <- cntrl[order(plate_MEASURED_DAYS)]
  cntrl[, plate_order := 1:.N]
  
  # Split into three color groups (ggnewscale was not behaving)
  gg_no_outlier <- cntrl[outlier == "no"]
  gg_high_outlier <- cntrl[outlier == "high"]
  gg_low_outlier <- cntrl[outlier == "low"]
  
  g <- ggplot(cntrl, aes(x=factor(plate_order), shape=factor(control_sample))) +
    scale_shape_manual(name = "Control Sample (4-pairs used across all plates)",
      values=c( # Control samples are paired
				"190404 (Well A01)"=0, "190328 (Well H12)"=22,
				"190425 (Well A01)"=1, "190508 (Well H12)"=21,
				"180829 (Well A01)"=2, "180827 (Well H12)"=24,
				"180830 (Well A01)"=5, "180831 (Well H12)"=23
			)) +

    geom_point(data=gg_no_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#525252") +
    geom_point(data=gg_no_outlier[plate_position == "A01"], aes(y=value), color="#525252", size=0.45, stroke=0.2) +

    geom_point(data=gg_high_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#a50f15") +
    geom_point(data=gg_high_outlier[plate_position == "A01"], aes(y=value), color="#a50f15", size=0.45, stroke=0.2) +
    
    geom_point(data=gg_low_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#08519c") +
    geom_point(data=gg_low_outlier[plate_position == "A01"], aes(y=value), color="#08519c", size=0.45, stroke=0.2) +

    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(paste0("Concentration of control sample")) +
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(shape = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1, colour="black", fill="#525252")))
  
  ggsave(g, width=width, height=width, units="in", file=fname)
}

internal_controls("diagnostic_plots/internal_control_samples_outlier_plates.png")

# Show on log scale too
cntrl[log_offset, on = .(variable), value := log(value - right_shift + log_offset)]
cntrl[outlier_shift, on = .(variable), value := value - shift]

internal_controls("diagnostic_plots/internal_control_samples_log_outlier_plates.png")


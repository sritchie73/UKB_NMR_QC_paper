library(data.table)
library(ukbnmr)
library(foreach)
library(doMC)
library(ggplot2)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Filter to samples present in UKB raw data:
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(eid_7439, sample_id, visit)], on = .(eid_7439, sample_id, visit), nomatch=0]

# Filter to standard set of 249 biomarkers and ratios
dat <- dat[variable %in% ukbnmr::nmr_info[(Nightingale), Biomarker]]

# Load extended technical information
tech_info <- fread("data/raw/processed/sample_information.txt")
tech_info <- tech_info[!(removed)]

# Drop columns that are either (1) not technical covariates, (2) only have a single value 
# across the dataset, or (3) are date/time stamps with other columns already created that 
# have appropriate numeric representations.
tech_info[, c("uid", "removed", "removal_reason", "Collaborator", "Collection.Date", "Sample.Type", 
              "Sample.Source", "Gender", "Dispense.Sequence", "tag", "plate_MEASURED_DATE",
              "hours_on_machine", "hours_on_machine2", 
               grep("date_", names(tech_info), value=TRUE), grep("time_", names(tech_info), value=TRUE)) := NULL]

# Filter to samples present in UKB raw data:
tech_info <- tech_info[sinfo[,.(sample_id, visit)], on = .(sample_id, visit), nomatch=0]

# Function for converting a set of numbers to a factor, ordering levels by group size
factor_by_size <- function(x) {
  group_sizes <- as.data.table(table(x))
  group_sizes <- group_sizes[order(-N)]
  factor(x, levels=group_sizes$x)
}

# Function to cut vector x into N bins of equal duration
bin_duration <- function(x, n=10) {
  offset_x <- x + 1 # offset by 1 so we can apply ceiling function
  bin_size <- max(offset_x)/n
  as.integer(ceiling(offset_x/bin_size))
}

# Set column types for modelling
tech_info[, plate_position := factor_by_size(plate_position)]
tech_info[, well_row := factor_by_size(well_row)]
tech_info[, well_column := factor_by_size(well_column)]
tech_info[, spectrometer := factor_by_size(spectrometer)]
tech_info[, batch := factor_by_size(batch)]
tech_info[, Tecan.Name := factor_by_size(Tecan.Name)]

for (cn in grep("days_", names(tech_info), value=TRUE)) {
  tech_info[, c(cn) := bin_duration(tech_info[[cn]])]
}

for (cn in grep("hour_of_day_", names(tech_info), value=TRUE)) {
  tech_info[, c(cn) := bin_duration(tech_info[[cn]])]
}

# For plate, split into bins by spectrometer (running as factor for each plate
# takes too long, >18 hours with 10 cores).
tech_info[, plate_bin := bin_duration(as.integer(factor(plate_MEASURED_DAYS))),  by=spectrometer]
tech_info[, plate_bin := plate_bin + (as.integer(spectrometer)-1)*10]
tech_info[, plate_bin := factor_by_size(plate_bin)]
tech_info[, c("plate_id", "plate_MEASURED_DAYS") := NULL]

# Missing data for 609 samples on Tecan.Tip.
# We want these NAs to be treated as a separate group, not 
# dropped when examining variance explained.
tech_info[, Tecan.Tip := as.character(Tecan.Tip)]
tech_info[is.na(Tecan.Tip), Tecan.Tip := ""]
tech_info[, Tecan.Tip := factor_by_size(Tecan.Tip)]

# Fit linear regression and extract r2
tech_r2 <- foreach(biom = nmr_info[(Nightingale), Biomarker], .combine=rbind) %:% 
  foreach(techv = names(tech_info), .combine=rbind) %dopar% {
    if (techv %in% c("sample_id", "visit")) return(NULL)
    this_dat <- dat[variable == biom] 
    this_tech <- tech_info[, .SD, .SDcols=c("sample_id", "visit", techv)]
    setnames(this_tech, techv, "covar")
    this_dat <- this_dat[this_tech, on = .(sample_id, visit)]

	  raw_r2 = summary(lm(raw ~ covar, data=this_dat))[["r.squared"]]
    adj4_r2 = summary(lm(adj4 ~ covar, data=this_dat))[["r.squared"]]
    adj5_r2 = summary(lm(adj5 ~ covar, data=this_dat))[["r.squared"]]
    
    data.table(biomarker=biom, covariate=techv, raw_r2 = raw_r2, adj4_r2 = adj4_r2, adj5_r2 = adj5_r2)
}

# Add nice variable names for plotting
tech_r2[covariate == "batch", print_name := "Shipping batch"]
tech_r2[covariate == "plate_bin", print_name := "96-well plate (drift over time)"]
tech_r2[covariate == "plate_position", print_name := "Well position (A2-H11)"]
tech_r2[covariate == "well_row", print_name := "Well row (A-H)"]
tech_r2[covariate == "well_column", print_name := "Well column (1-12)"]
tech_r2[covariate == "Tecan.Name", print_name := "Tecan aliquoting robot"]
tech_r2[covariate == "Tecan.Tip", print_name := "Tecan aliquot tip"]
tech_r2[covariate == "spectrometer", print_name := "Spectrometer"]

tech_r2[covariate == "days_DISPATCHED_AT", print_name := "Date sample dispatched from UK Biobank"]
tech_r2[covariate == "days_ARRIVED_AT", print_name := "Date sample arrived at Nightingale"]
tech_r2[covariate == "days_FROZEN_AT", print_name := "Date sample frozen"]
tech_r2[covariate == "days_DEFROSTED_AT", print_name := "Date sample defrosted"]
tech_r2[covariate == "days_CENTRIFUGED_AT", print_name := "Date sample centrifuged"]
tech_r2[covariate == "days_PREPARED_AT", print_name := "Date sample prepared"]
tech_r2[covariate == "days_MEASURED_AT", print_name := "Date sample measured"]

tech_r2[covariate == "hour_of_day_DISPATCHED_AT", print_name := "Time of day sample dispatched from UK Biobank"]
tech_r2[covariate == "hour_of_day_ARRIVED_AT", print_name := "Time of day sample arrived at Nightingale"]
tech_r2[covariate == "hour_of_day_FROZEN_AT", print_name := "Time of day sample frozen"]
tech_r2[covariate == "hour_of_day_DEFROSTED_AT", print_name := "Time of day sample defrosted"]
tech_r2[covariate == "hour_of_day_CENTRIFUGED_AT", print_name := "Time of day sample centrifuged"]
tech_r2[covariate == "hour_of_day_PREPARED_AT", print_name := "Time of day sample prepared"]
tech_r2[covariate == "hour_of_day_MEASURED_AT", print_name := "Time of day sample measured"]

tech_r2[covariate == "dispatch_to_arrival", print_name := "Time taken to ship from UK Biobank to Nightingale Health"]
tech_r2[covariate == "arrival_to_freeze", print_name := "Time between sample arrival and freezing"]
tech_r2[covariate == "freeze_duration", print_name := "Time sample frozen for"]
tech_r2[covariate == "defrost_to_centrifuge", print_name := "Time between sample defrosting and centrifugation"]
tech_r2[covariate == "centrifuge_to_prep", print_name := "Time between sample centrifugation and sample preparation"]
tech_r2[covariate == "prep_to_measured", print_name := "Time between sample preparation and measurement"]

# Write out
fwrite(tech_r2, sep="\t", quote=FALSE, file="paper_output/tech_r2.txt")

# Does the r2 change after excluding outlier plates?
g <- ggplot(tech_r2) +
  aes(x=adj4_r2*100, y=adj5_r2*100) +
  geom_abline(intercept=0, slope=1, color="red", linetype=2) +
  geom_point(shape=19, size=0.5) + 
  facet_wrap(~ print_name, scales="free", labeller = labeller(print_name = label_wrap_gen(18))) +
  xlab("Variance Explained (%) in Post-QC data") +
  ylab("Variance Explained (%) after excluding outlier plates") +
  theme_bw() +
	theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
			axis.text.x = element_text(size=6), axis.text.y=element_text(size=7),
			axis.title.x = element_text(size=7), legend.position="none",
			strip.background=element_blank(), strip.text=element_text(size=6, face=2))
ggsave(g, width=7.2, height=7.2, file="paper_output/tech_r2_diff_outlier_plates.pdf") 

# Prepare for plotting
tech_r2 <- melt(tech_r2, id.vars=c("biomarker", "covariate", "print_name"), value.name="r2", variable.name="data")
tech_r2[, data := fcase(
  data == "raw_r2", "Raw data", 
  data == "adj4_r2", "Post-QC data",
  data == "adj5_r2", "Post-QC data excluding outlier plates"
)]
tech_r2[, pct := r2 * 100]

# Get row and column order
tech_r2[, data := factor(data, levels=c("Raw data", "Post-QC data", "Post-QC data excluding outlier plates"))]
tech_r2[, print_name := factor(print_name, levels=rev(c(
    "Tecan aliquoting robot", "Tecan aliquot tip",
    "96-well plate (drift over time)", "Well position (A2-H11)", "Well row (A-H)", "Well column (1-12)",
    "Shipping batch", "Date sample dispatched from UK Biobank", "Time of day sample dispatched from UK Biobank",
    "Date sample arrived at Nightingale", "Time of day sample arrived at Nightingale", "Time taken to ship from UK Biobank to Nightingale Health",
    "Date sample frozen", "Time of day sample frozen", "Time between sample arrival and freezing",
    "Date sample defrosted", "Time of day sample defrosted", "Time sample frozen for",
    "Date sample centrifuged", "Time of day sample centrifuged", "Time between sample defrosting and centrifugation",
    "Date sample prepared", "Time of day sample prepared", "Time between sample centrifugation and sample preparation",
    "Date sample measured", "Time of day sample measured", "Time between sample preparation and measurement",
    "Spectrometer"
  )))]

# Plot variance explained in raw data by technical print_names
g <- ggplot(tech_r2[data != "Post-QC data excluding outlier plates"]) +  # redundant
  aes(x = pct, y = print_name, fill = print_name) +
  facet_wrap(~ data) +
  xlab("Variance explained (%)") + ylab("") +
  geom_boxplot(outlier.size=0.6, color="black", size=0.3, outlier.stroke=0, outlier.color="#525252") +
  geom_vline(xintercept=1, linetype=2, color="red") +
  theme_bw() +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
        axis.text.x = element_text(size=6), axis.text.y=element_text(size=7),
        axis.title.x = element_text(size=7), legend.position="none",
        strip.background=element_blank(), strip.text=element_text(size=7, face=2))
ggsave(g, width=7.2, height=3.6, file="paper_output/tech_r2_boxplot.pdf")

# Get the maximum variance explained by any factor in the raw data for each biomarker and make a density plot
r2_max <- tech_r2[data == "Raw data", .SD[which.max(r2)], by=.(biomarker)]

g <- ggplot(r2_max, aes(x=r2)) +
  geom_density() + geom_rug() +
  scale_x_continuous(name="Maximum Variance Explained (%) by any technical covariate in raw data", expand=expansion(mult=0.01)) +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
        axis.text.x = element_text(size=6), axis.text.y=element_text(size=6),
        axis.title = element_text(size=7))
ggsave(g, width=7.2, height=1.5, file="paper_output/max_r2_density.pdf")


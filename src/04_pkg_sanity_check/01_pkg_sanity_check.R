library(data.table)
library(ukbnmr)
library(ggplot2)
library(ggthemes)
library(palettetown)

# Make output directory
system("mkdir -p pkg_check", wait=TRUE)

# Load in full pre-release raw data and extended sample information
raw <- fread("data/raw/processed/UKB_phase1_results.txt")
sinfo <- fread("data/raw/processed/harmonized_sample_information.txt")

# Load in and extract raw data available directly from UK Biobank
ukb_raw <- fread("data/raw/ukbiobank/extracted/nightingale.csv")
ukb_raw <- extract_biomarkers(ukb_raw)

# The pre-release raw data contains blind duplicates and samples that
# may no longer be in the UK Biobank data, so we need to identify the
# matching samples
raw[sinfo, on = .(sample_id, visit, plate_id, plate_position, spectrometer), eid_7439 := eid_7439] 
raw[, visit_index := ifelse(visit == "Main Phase", 0L, 1L)]
raw <- melt(raw, id.vars=c("eid_7439", "visit_index", "visit", "sample_id", "plate_id", "plate_position", "spectrometer"),
            variable.name="biomarker", value.name="prerelease")

ukb_raw <- melt(ukb_raw, id.vars=c("eid", "visit_index"), variable.name="biomarker", value.name="ukb")
raw[ukb_raw, on = .(eid_7439=eid, visit_index, biomarker), ukb := ukb]

# Get best matching sample ([sample_id, visit] combination for each [eid, visit] combination
# with the most matching biomarker values.
matching_values <- raw[abs(prerelease - ukb) < 1e-5 | (is.na(prerelease) & is.na(ukb)), .N, by=.(eid_7439, visit, sample_id)]
best_match <- matching_values[,.SD[which.max(N)], by=.(eid_7439, visit)]
raw[best_match, on=.(eid_7439, visit, sample_id), best_matching_sample := TRUE]
raw[is.na(best_matching_sample), best_matching_sample := FALSE]

# Get mismatching values
mismatch <- raw[(best_matching_sample) & (
  abs(prerelease - ukb) > 1e-5 |
  (is.na(prerelease) & !is.na(ukb)) |
  (!is.na(prerelease) & is.na(ukb))
)]

if(mismatch[,.N] > 0) {
	fwrite(mismatch, sep="\t", quote=FALSE, file="pkg_check/raw_concentration_mismatch_among_matched_samples.txt")
}

# At the moment this is just one participant whose XXL_VLDL_*_pct are NA in the prerelease
# data but non-missing in the data now on UK Biobank. This [eid, visit] combination have
# two samples in the pre-release data, but the values clearly mismatch for the other sample
# not flagged as the best match.

# Add in these to the extended sample information sheet
best_match <- unique(raw[(best_matching_sample),.(eid_7439, visit, sample_id, plate_id, plate_position, spectrometer)])
sinfo[, c("in_ukb_raw", "withdrawn") := NULL]
sinfo[best_match, on=.(eid_7439, visit, sample_id, plate_id, plate_position, spectrometer), in_ukb_raw := TRUE]
sinfo[is.na(in_ukb_raw), in_ukb_raw := FALSE]

# Check for samples that are in the UKB raw, but failed QC in our prerelease data
fail_qc <- sinfo[(in_ukb_raw) & (removed)]
fail_qc <- sinfo[fail_qc[,.(sample_id, visit)], on = .(sample_id, visit)]
fail_qc <- fail_qc[,.(uid, in_ukb_raw, sample_id, plate_id, plate_position, visit, spectrometer, removed, removal_reason)]

if (fail_qc[,.N] > 0) { 
  fwrite(fail_qc, sep="\t", quote=FALSE, file="pkg_check/samples_failing_qc_but_kept_by_ukb.txt")
}

# Find samples that are in the pre-release data, but not in the latest UKB data
# Obtain list of all sample withdrawals
withdrawal_files <- list.files(recursive=TRUE, path="data/ceu_curated_phenotypes/Withdrawals/", pattern="*.csv", full.names=TRUE)
withdrawals <- unique(rbindlist(lapply(withdrawal_files, fread, header=FALSE)))
withdrawals <- withdrawals[V1 %like% "^[0-9]*$"]
setnames(withdrawals, "eid_7439")
withdrawals[, eid_7439 := as.integer(eid_7439)]
sinfo[, withdrawn := FALSE]
sinfo[withdrawals, on = .(eid_7439), withdrawn := TRUE]

missing <- sinfo[!(in_ukb_raw) & !(removed) & !(blind_duplicate) & !(withdrawn)]
missing <- missing[,.(uid, sample_id, plate_id, plate_position, visit, spectrometer)]
if (missing[,.N] > 0) {
  fwrite(missing, sep="\t", quote=FALSE, file="pkg_check/samples_missing_not_withdrawals_from_ukb.txt")
}

# Now, Load the data we adjusted for technical variation in this project
adj <- fread("data/tech_qc/postqc_excluding_outlier_plates.txt")

# Get output from package applied to most recent raw data released by UKB
ukb <- fread("data/my_curated_phenotypes/NMR_metabolomics/output/nmr_techadj.txt")
adj_sinfo <- fread("data/my_curated_phenotypes/NMR_metabolomics/output/sample_processing_information.txt")

# Harmonize and compare key sample processing/QC variables:
sinfo[, visit_index := ifelse(visit == "Main Phase", 0L, 1L)]
adj_bin <- fread("data/tech_qc/sample_information.txt")
sinfo[adj_bin, on = .(uid=raw_row_id), plate_measured_bin := plate_measured_bin]

adj_sinfo[sinfo[(in_ukb_raw)], on = .(eid=eid_7439, visit_index), 
  c("Prerelease.Sample.ID", "Prerelease.Plate.ID", "Prerelease.Plate.Position", "Prerelease.Spectrometer", 
    "Prerelease.Prep.to.Measure.Duration", "Prerelease.Spectrometer.Date.Bin") :=
  .(sample_id, plate_id, plate_position, spectrometer, prep_to_measured, plate_measured_bin)]

fwrite(adj_sinfo, sep="\t", quote=FALSE, file="pkg_check/sample_qc_info_compare.txt")

# Check spectrometers - should be six all up when matching UKB Spectrometer field to 
# prerelease spectrometer number
spec <- adj_sinfo[,.N,by=.(Spectrometer, Prerelease.Spectrometer)][order(Spectrometer)]
fwrite(spec, sep="\t", quote=FALSE, file="pkg_check/spectrometer_map.txt")

# Check for plate mismatch
plate_mismatch <- adj_sinfo[Shipment.Plate != paste("Plate", Prerelease.Plate.ID)]
plate_mismatch <- plate_mismatch[,.(Prerelease.Sample.ID, visit_index, Shipment.Plate, Prerelease.Plate.ID, Prerelease.Spectrometer)]
if (plate_mismatch[,.N] > 0) {
  fwrite(plate_mismatch, sep="\t", quote=FALSE, file="pkg_check/plate_mismatch.txt")
}

# Check for well mismatch
well_mismatch <- adj_sinfo[Well.Position.Within.Plate != Prerelease.Plate.Position]
well_mismatch <- well_mismatch[,.(Prerelease.Sample.ID, visit_index, Prerelease.Plate.ID, Well.Position.Within.Plate, Prerelease.Plate.Position, Prerelease.Spectrometer)]
if (well_mismatch[,.N] > 0) {
  fwrite(well_mismatch, sep="\t", quote=FALSE, file="pkg_check/well_mismatch.txt")
}

# Check for differences in time between sample prep and sample measurement
duration_check <- adj_sinfo[abs(Prep.to.Measure.Duration - Prerelease.Prep.to.Measure.Duration) > 0.01]
duration_check <- duration_check[,.(Prerelease.Sample.ID, visit_index, Sample.Prepared.Date.and.Time, Sample.Measured.Date.and.Time, Prep.to.Measure.Duration, Prerelease.Prep.to.Measure.Duration)]
duration_check[sinfo[(in_ukb_raw)], on = .(Prerelease.Sample.ID=sample_id, visit_index), Prerelease.Prepared := paste(date_PREPARED_AT, time_PREPARED_AT)]
duration_check[sinfo[(in_ukb_raw)], on = .(Prerelease.Sample.ID=sample_id, visit_index), Prerelease.Measured := paste(date_MEASURED_AT, time_MEASURED_AT)]
if (duration_check[,.N] > 0) {
  fwrite(duration_check, sep="\t", quote=FALSE, file="pkg_check/duration_mismatch.txt")
}

# Check for differences in spectrometer bin over time
bin_check <- adj_sinfo[Prerelease.Spectrometer.Date.Bin != Spectrometer.Date.Bin %% 10 + 1]
bin_check <- bin_check[, .(eid, visit_index, Spectrometer, Prerelease.Spectrometer.Date.Bin, Spectrometer.Date.Bin=Spectrometer.Date.Bin %% 10 + 1)]
bin_check[, bin_diff := Spectrometer.Date.Bin - Prerelease.Spectrometer.Date.Bin]
if (bin_check[abs(bin_diff) == 1,.N] > 0) {
  fwrite(bin_check[abs(bin_diff) == 1], sep="\t", quote=FALSE, file="pkg_check/spectrometer_date_bin_boundary_changes.txt")
}
if (bin_check[abs(bin_diff) > 1,.N] > 0) {
  fwrite(bin_check[abs(bin_diff) > 1], sep="\t", quote=FALSE, file="pkg_check/spectrometer_date_bin_mismatch.txt")
}

# Compare adjusted values
adj[, visit_index := ifelse(visit == "Main Phase", 0L, 1L)]
adj <- melt(adj, id.vars=c("eid_7439", "visit_index", "visit", "sample_id", "plate_id", "plate_position", "spectrometer"),
            variable.name="biomarker", value.name="prerelease")

ukb[sinfo[(in_ukb_raw)], on = .(eid=eid_7439, visit_index), c("visit", "sample_id", "plate_id", "plate_position", "spectrometer") :=
    .(visit, sample_id, plate_id, plate_position, spectrometer)]
ukb <- melt(ukb, id.vars=c("eid", "visit_index", "visit", "sample_id", "plate_id", "plate_position", "spectrometer"),
            variable.name="biomarker", value.name="ukb")

adj[ukb, on = .(eid_7439=eid, visit_index, visit, sample_id, plate_id, plate_position, spectrometer, biomarker), c("ukb") := i.ukb]

# Compute correlations
corr <- adj[is.finite(prerelease) & is.finite(ukb), .(pearson=cor(prerelease, ukb), spearman=cor(prerelease, ukb, method="spearman")), by=biomarker]
corr <- corr[order(spearman)]
fwrite(corr, sep="\t", quote=FALSE, file="pkg_check/prerelease_vs_pkg_adjusted_values_correlation.txt")

# Plot differences
n_bio <- adj[,.N,by=biomarker][,.N]
n_col <- ceiling(sqrt(n_bio))
width <- n_col * 1.4

g <- ggplot(adj) +
  aes(x=prerelease, y=ukb) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_x_continuous(name="Adjusted biomarker in pre-release data") +
  scale_y_continuous(name="Adjusted biomarker from UKB release data using ukbnmr R package") +
  facet_wrap(~ biomarker, ncol=n_col, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=width, height=width, units="in", file="pkg_check/prerelease_vs_pkg_adjusted_values.png")




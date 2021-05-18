library(data.table)
library(ukbnmr) # remotes::install_github("sritchie73/ukbnmr")

# Load raw data provided by Nightingale that we have applied sample QC to:
proc <- fread("data/raw/processed/UKB_phase1_results.txt", colClasses=c("sample_id"="character", "plate_id"="character"))
sinfo <- fread("data/raw/processed/sample_information.txt", colClasses=c("sample_id"="character", "plate_id"="character"))
tags <- fread("data/raw/processed/biomarker_QC_tags.txt", colClasses=c("sample_id"="character"))

# Link to P7439 sample ids
ukblink <- fread("data/raw/nightingale/Cambridge_Nightingale_bridge.tsv", colClasses="character")
sinfo[ukblink, on = .(sample_id, plate_id, plate_position), eid_7439 := eid_7439]
proc[ukblink, on = .(sample_id, plate_id, plate_position), eid_7439 := eid_7439]
tags[sinfo[!(removed)], on = .(sample_id, visit), eid_7439 := eid_7439]

# Flag samples that are blind duplicates
sinfo[, blind_duplicate := FALSE]
sinfo[sinfo[,.N,by=.(eid_7439, visit)][N > 1], on = .(eid_7439, visit), blind_duplicate := TRUE]

# Load raw data available through UKB
ukb <- fread("data/raw/ukbiobank/extracted/nightingale.csv", colClasses=c("eid"="character")) # includes biomarker and sample QC Flags
ukb_nmr <- extract_biomarkers(ukb) # extract biomarker fields
ukb_nmr[, eid := as.character(eid)]
ukb_nmr[, visit_index := ifelse(visit_index == 0L, "Main Phase", "Repeat Assessment")]

# Find which samples have been released by UKB
biomarkers <- intersect(names(ukb_nmr), names(proc))
proc <- melt(proc, id.vars=c("eid_7439", "sample_id", "visit"), measure.vars=biomarkers)
ukb_nmr <- melt(ukb_nmr, id.vars=c("eid", "visit_index"), measure.vars=biomarkers)

proc[, in_ukb := FALSE]
proc[ukb_nmr, on = .(eid_7439=eid, visit=visit_index, variable, value), in_ukb := TRUE]
ukb_nmr[, in_proc := FALSE]
ukb_nmr[proc, on = .(eid=eid_7439, visit_index=visit, variable, value), in_proc := TRUE]

# Make note in sinfo of the types of overlap:
sinfo[, in_ukb_raw := "No biomarkers"]
overlap_type <- proc[,.(all_in_ukb=all(in_ukb), any_in_ukb=any(in_ukb)), by=.(eid_7439, sample_id, visit)]
sinfo[overlap_type[(all_in_ukb)], on=.(eid_7439, sample_id, visit), in_ukb_raw := "All 168 biomarkers"]
sinfo[overlap_type[!(all_in_ukb) & (any_in_ukb)], on=.(eid_7439, sample_id, visit), in_ukb_raw := "Some biomarkers"]
sinfo[(removed), in_ukb_raw := "No biomarkers"]

# What samples in UKB do we not have?
miss_proc <- ukb_nmr[!(in_proc), .N, by=.(eid, visit_index)]
miss_proc <- miss_proc[!unique(proc[,.(eid_7439, visit)]), on = .(eid=eid_7439, visit_index=visit)]
miss_proc <- sinfo[miss_proc, on = .(eid_7439=eid, visit=visit_index)]
miss_proc <- miss_proc[, .(eid_7439, sample_id, visit, plate_id, plate_position, spectrometer, removal_reason)]
fwrite(miss_proc, sep="\t", quote=FALSE, file="data/raw/processed/samples_in_ukb_failing_qc.txt")

# What measurements are in UKB but cannot be matched?
miss_meas <- ukb_nmr[!(in_proc)][!miss_proc, on = .(eid=eid_7439, visit_index=visit)]

# Are any of them just numeric error? Note we need to do one-many matches here due to 
# presence of blind duplicates in raw data
miss_meas <- miss_meas[proc, on = .(eid=eid_7439, visit_index=visit, variable), nomatch=0]
setnames(miss_meas, c("value", "i.value"), c("ukb.value", "prerelease.value"))
miss_meas[, row := .I]
numeric_match <- miss_meas[,.(numeric_equal=isTRUE(all.equal(ukb.value, prerelease.value)), ukb.value, prerelease.value), by=.(row, eid, visit_index, variable)]
numeric_match <- numeric_match[(numeric_equal)]
ukb_nmr[numeric_match, on = .(eid, visit_index, variable), in_proc := TRUE]

# Pull out remaining measurements that cannot be matched?
miss_meas <- ukb_nmr[!(in_proc)][!miss_proc, on = .(eid=eid_7439, visit_index=visit)]
miss_meas <- miss_meas[proc, on = .(eid=eid_7439, visit_index=visit, variable), nomatch=0]
setnames(miss_meas, c("value", "i.value"), c("ukb.value", "prerelease.value"))

# Pull out ones which are wrong order of magnitude for some reason
magnitude <- miss_meas[ukb.value == prerelease.value / 10]
magnitude[sinfo[!(removed)], on = .(eid=eid_7439, visit_index=visit), 
 c("sample_id", "plate_id", "plate_position", "spectrometer") :=
 .(sample_id, plate_id, plate_position, spectrometer)]
magnitude <- magnitude[, .(eid_7439=eid, sample_id, visit=visit_index, plate_id, plate_position, spectrometer, 
                           biomarker=variable, ukb.value, prerelease.value)]
fwrite(magnitude, sep="\t", quote=FALSE, file="data/raw/processed/values_div10_in_ukb_release.txt")

# None left in UKB we can't match to our sample QC data

# Now update our sample QC data to include those that are numeric, but not equal, matches
proc[ukb_nmr, on = .(eid_7439=eid, visit=visit_index, variable), ukb.value := i.value]
proc[, row := .I]
proc[(in_ukb), numeric_equal := TRUE]
proc[!(in_ukb) & !is.na(ukb.value), numeric_equal := isTRUE(all.equal(value, ukb.value)), by=.(row)]
proc[!(in_ukb) & (numeric_equal), in_ukb := TRUE]

# And also those divided by 10 for some reason (we just want to know which samples are
# present in the UKB release in this instance)
proc[!(in_ukb) & value == ukb.value * 10, in_ukb := TRUE]

# Any missing in both datasets should be counted as a match
proc[is.na(value) & is.na(ukb.value), in_ukb := TRUE]

# For blind duplicates, there are some cases where a biomarker has identical values across
# samples for some biomarkers. Identify and drop these, keeping the primary sample
dupes <- sinfo[(blind_duplicate), .(eid_7439, sample_id, visit)]
dupes <- proc[dupes, on = .(eid_7439, sample_id, visit)]
primary <- dupes[(in_ukb), .N, by=.(eid_7439, sample_id, visit)][N == 168]
alt <- dupes[(in_ukb)][!primary, on = .(eid_7439, sample_id, visit)]
proc[alt, on = .(eid_7439, sample_id, visit, variable), in_ukb := FALSE

# Collate overlap information again
sinfo[, in_ukb_raw := "No biomarkers"]
overlap_type <- proc[,.(all_in_ukb=all(in_ukb), any_in_ukb=any(in_ukb & !is.na(value) & !is.na(ukb.value))), by=.(eid_7439, sample_id, visit)]
sinfo[overlap_type[(all_in_ukb)], on=.(eid_7439, sample_id, visit), in_ukb_raw := "All 168 biomarkers"]
sinfo[overlap_type[!(all_in_ukb) & (any_in_ukb)], on=.(eid_7439, sample_id, visit), in_ukb_raw := "Some biomarkers"]
sinfo[(removed), in_ukb_raw := "No biomarkers"]

# Pull out data on samples present in our sample QC version, but not in the released UKB data,
# where those samples also are not blind duplicates
dropped <- sinfo[!(removed) & !(blind_duplicate) & in_ukb_raw == "No biomarkers"]
dropped <- dropped[,.(eid_7439, sample_id, visit, plate_id, plate_position, spectrometer, QC_flag=tag)]

# Are they the latest withdrawals? (No, only 5 / 100)
withdrawals <- fread("data/ceu_curated_phenotypes/Withdrawals/20210221_withdrawals.csv", colClasses="character")
sinfo[, withdrawn := FALSE]
sinfo[withdrawals, on = .(eid_7439=n_eid), withdrawn := TRUE]
dropped <- sinfo[!(removed) & !(blind_duplicate) & !(withdrawn) & in_ukb_raw == "No biomarkers"]
dropped <- dropped[,.(eid_7439, sample_id, visit, plate_id, plate_position, spectrometer, QC_flag=tag)]

# What biomarker QC tags do they have?
dropped <- proc[dropped[, .(eid_7439, visit)], on = .(eid_7439, visit), nomatch=0]
dropped <- dropped[tags, on = .(eid_7439, sample_id, visit, variable=biomarker), nomatch=0]

# 95 dropped samples with no obvious explanation - no obvious sample QC flags or
# biomarker QC flags.
dropped <- dropped[,.(eid_7439, sample_id, visit, biomarker=variable, prerelease.value=value, QC_flag=QC_tag)]
fwrite(dropped, sep="\t", quote=FALSE, file="data/raw/processed/dropped_in_ukb_biomarker_QC_flags.txt")

dropped <- sinfo[!(removed) & !(blind_duplicate) & !(withdrawn) & in_ukb_raw == "No biomarkers"]
dropped <- dropped[,.(eid_7439, sample_id, visit, plate_id, plate_position, spectrometer, QC_flag=tag)]
fwrite(dropped, sep="\t", quote=FALSE, file="data/raw/processed/dropped_in_ukb_sampless.txt")

# Write out selected sample for each blind duplicate
dupes <- sinfo[(blind_duplicate) & in_ukb_raw == "All 168 biomarkers", .(eid_7439, sample_id, visit)]
fwrite(dupes, sep="\t", quote=FALSE, file="data/raw/processed/selected_blind_duplicate_samples.txt")

# Write out updated sample information sheet
fwrite(sinfo, sep="\t", quote=FALSE, file="data/raw/processed/harmonized_sample_information.txt")



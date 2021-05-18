library(data.table)

# Load wide-format raw data containing all quantified biomarkers
raw <- fread("data/raw/nightingale/UKB-phase1__results.tsv", colClasses=list("character"=1:5, "numeric"=6:254))

# Load and add internal control samples
cntrl <- fread("data/raw/nightingale/NGH_control_samples.tsv", colClasses=list("character"=1:9, "numeric"=10:258), na.strings=c("\"NA\"", "NA"))
cntrl[,visit := "Internal Control"]
cntrl[, names(cntrl)[5:9] := NULL]
raw <- rbind(raw, cntrl)

# Add well position info
raw[, well_row := substr(plate_position, 1, 1)]
raw[, well_column := as.integer(substr(plate_position, 2, 3))]

# Add unique row identifier (so we can identify and remove duplicates)
raw[, uid := .I]

# split out sample info, and melt measurements:
sinfo <- raw[, .(uid, sample_id, plate_id, plate_position, well_row, well_column, visit, spectrometer)]
raw[, c("sample_id", "plate_id", "plate_position", "well_row", "well_column", "visit", "spectrometer") := NULL]
raw <- melt(raw, id.vars=c("uid"))

# Add columns to sinfo for flagging sample removal
sinfo[, removed := FALSE]
sinfo[, removal_reason := NA_character_]

# Remove 40 samples that were accidentally measured despite having insufficient sample material:
remove <- fread("data/raw/nightingale/ukb_phase1_samples_to_be_removed.tsv", colClasses="character")
remove[, plate_id := paste0("0", plate_id)] # to match other files
sinfo[remove, on = .(sample_id, plate_id, plate_position, visit), c("removed", "removal_reason") := .(TRUE, "Insufficient sample material")]

# Flag rows with all missing data for removal
all_NA <- raw[sinfo[!(removed), .(uid)], on = .(uid)][, .(all_NA=all(is.na(value))), by=.(uid)][(all_NA)]
sinfo[all_NA, on = .(uid), c("removed", "removal_reason") := .(TRUE, "All measurements missing")]

# Load and add information about aliquoting robot and shipping batch and other batch information
ff <- list.files(path="data/raw/nightingale/UKB_technical_variables/", full.names=TRUE)
ll <- lapply(ff, fread, check.names=TRUE, colClasses="character")
names(ll) <- basename(ff)
names(ll) <- gsub("Nightingale_", "", names(ll))
names(ll) <- gsub(".csv", "", names(ll))
tech <- rbindlist(idcol="batch", ll)
rm(ll)
sinfo <- merge(sinfo, tech, by.x=c("sample_id", "plate_id", "plate_position", "visit"), by.y=c("Patient.ID", "Sample.ID", "Plate.Position", "Visit"), all.x=TRUE)

# Load and add sample QC tags from the Nightingale QC pipeline
stags <- fread("data/raw/nightingale/UKB-phase1__sample_tags_and_notes.tsv", colClasses="character")
stags[, uid := .I] # same row order as raw data file
stags[, names(stags)[1:5] := NULL]
stags <- melt(stags, id.vars="uid")
stags <- stags[!is.na(value) & value != ""]
stags[value == "1", value := variable]
stags <- stags[,.(uid, tag=value)]

sinfo[, tag := NA_character_]
sinfo[stags[,.(tag = paste(tag, collapse=";")), by=.(uid)], on = .(uid), tag := i.tag]

# Load extended date/time based technical covariates
tech <- fread("data/raw/nightingale/NGH_sample_handling.tsv", colClasses="character")
tech[sinfo, on = .(sample_id, visit), uid := uid]

# Fix some duplicates - assume same order as raw data file
tech[sample_id == "1219449" & visit == "Main Phase", uid := c(22004, 22005)]
tech[sample_id == "1556057" & visit == "Main Phase", uid := c(54813, 54814)]
tech[sample_id == "2553581" & visit == "Main Phase", uid := c(102853, 102854)]
tech[sample_id == "4478891" & visit == "Main Phase", uid := c(70674, 70675)]

# Melt
tech <- melt(tech, id.vars=c("uid", "sample_id", "visit"))

# Add Internal Controls
cntrl <- fread("data/raw/nightingale/NGH_control_samples.tsv", colClasses=list("character"=1:9, "numeric"=10:258), na.strings=c("\"NA\"", "NA"))
cntrl[,visit := "Internal Control"]
cntrl <- cntrl[, .(sample_id, visit, plate_id, plate_position, spectrometer, FROZEN_AT, DEFROSTED_AT, CENTRIFUGED_AT, PREPARED_AT, MEASURED_AT)]
cntrl[sinfo, on = .(sample_id, visit, plate_id, plate_position, spectrometer), uid := uid]

# Fix some duplicates - assume same order as raw data file
cntrl[plate_id == "0490000006069" & plate_position == "H12", uid := c(127203, 127207)]
cntrl[plate_id == "0490000006069" & plate_position == "A01", uid := c(127620, 127630)]

# Melt
cntrl <- cntrl[,.(uid, sample_id, visit, FROZEN_AT, DEFROSTED_AT, CENTRIFUGED_AT, PREPARED_AT, MEASURED_AT)]
cntrl <- melt(cntrl, id.vars=c("uid", "sample_id", "visit"))

# Combine and drop missing values
tech <- rbind(tech, cntrl)
tech <- tech[!is.na(value)]

# Split out dates and times
tech[, date := as.IDate(gsub("T*", "", value))]
tech[, time := as.ITime(gsub("Z$", "", gsub("*T", "", value)))]

# Also get numeric representation for later removal of technical variation
# For date, get # days relative to first observed date
# For time, get hours relative to midnight
earliest_date <- tech[,min(date)]
tech[, days := date - earliest_date]
tech[, hour_of_day := hour(time) + minute(time)/60 + second(time)/3600]

# Expand to wide so we can compute durations
tech <- dcast(tech, uid ~ variable, value.var=c("date", "time", "days", "hour_of_day"))

# Compute durations between events (hours)
duration <- function(days1, hour_of_day1, days2, hour_of_day2) {
  (days2 - days1)*24 + (hour_of_day2 - hour_of_day1)
}

tech[, dispatch_to_arrival := duration(days_DISPATCHED_AT, hour_of_day_DISPATCHED_AT, days_ARRIVED_AT, hour_of_day_ARRIVED_AT)]
tech[, arrival_to_freeze := duration(days_ARRIVED_AT, hour_of_day_ARRIVED_AT, days_FROZEN_AT, hour_of_day_FROZEN_AT)]
tech[, freeze_duration := duration(days_FROZEN_AT, hour_of_day_FROZEN_AT, days_DEFROSTED_AT, hour_of_day_DEFROSTED_AT)]
tech[, defrost_to_centrifuge := duration(days_DEFROSTED_AT, hour_of_day_DEFROSTED_AT, days_CENTRIFUGED_AT, hour_of_day_CENTRIFUGED_AT)]
tech[, centrifuge_to_prep := duration(days_CENTRIFUGED_AT, hour_of_day_CENTRIFUGED_AT, days_PREPARED_AT, hour_of_day_PREPARED_AT)]
tech[, prep_to_measured := duration(days_PREPARED_AT, hour_of_day_PREPARED_AT, days_MEASURED_AT, hour_of_day_MEASURED_AT)]

# Add to sinfo
sinfo <- merge(sinfo, tech, by="uid", all.x=TRUE)

# Compute time on machine
first_measure <- sinfo[,.SD[date_MEASURED_AT == min(date_MEASURED_AT)],by=.(plate_id, spectrometer)]
first_measure <- first_measure[,.SD[which.min(time_MEASURED_AT)], by=.(plate_id, spectrometer)]
sinfo[first_measure, on = .(plate_id, spectrometer), hours_on_machine := duration(i.days_MEASURED_AT, i.hour_of_day_MEASURED_AT, days_MEASURED_AT, hour_of_day_MEASURED_AT)]

# Get majority date for each plate measurement
plate_date <- sinfo[,.N, by=.(date_MEASURED_AT, days_MEASURED_AT, plate_id, spectrometer)][,.SD[which.max(N)], by=.(plate_id, spectrometer)]
sinfo[plate_date, on = .(plate_id, spectrometer), c("plate_MEASURED_DATE", "plate_MEASURED_DAYS") := .(i.date_MEASURED_AT, i.days_MEASURED_AT)]

# Resolve duplicate plate+position measurements and dupliate sample_id + visit measurements (for non-controls), 
# dropping one from each pair for reasons listed below.
sinfo[sample_id == "191127" & visit == "Internal Control", c("removed", "removal_reason") := .(TRUE, "Orphaned control pair")]
sinfo[sample_id == "190125" & visit == "Internal Control", c("removed", "removal_reason") := .(TRUE, "Orphaned control pair")]
sinfo[sample_id == "180817" & visit == "Internal Control", c("removed", "removal_reason") := .(TRUE, "Orphaned control pair")]
sinfo[sample_id == "190404" & plate_position == "H12", c("removed", "removal_reason") := .(TRUE, "Wrong control in well")]
sinfo[plate_id == "0490000006107" & spectrometer == "10176949", c("removed", "removal_reason") := .(TRUE, "Different spectrometer to rest of plate")]
sinfo[uid %in% c(127203, 127620), c("removed", "removal_reason") := .(TRUE, "Duplicate control with largest time difference to rest of plate")]
sinfo[uid %in% c(54814, 102854), c("removed", "removal_reason") := .(TRUE, "Duplicate sample with largest time difference to rest of plate")]

# Drop sample with tag "technical_error" on the advice of Nightingale (1 / 189 samples not already removed,
# the other 188 had all missing values)
sinfo[!(removed) & tag %like% "Technical_error", c("removed", "removal_reason") := .(TRUE, "Sample tagged technical error")]

# Recompute time on machine now that we have removed abberant duplicates
first_measure <- sinfo[!(removed),.SD[date_MEASURED_AT == min(date_MEASURED_AT)],by=.(plate_id, spectrometer)]
first_measure <- first_measure[,.SD[which.min(time_MEASURED_AT)], by=.(plate_id, spectrometer)]
sinfo[first_measure, on = .(plate_id, spectrometer), hours_on_machine2 := duration(i.days_MEASURED_AT, i.hour_of_day_MEASURED_AT, days_MEASURED_AT, hour_of_day_MEASURED_AT)]

# Drop all removed samples from raw
raw <- raw[!sinfo[(removed)], on = .(uid)]

# load and filter biomarker tags
tags <- fread("data/raw/nightingale/UKB-phase1__biomarker_tags_and_notes.tsv", colClasses="character")
tags[, uid := .I] # same row order as raw data file
tags[, names(tags)[1:5] := NULL]
tags <- melt(tags, id.vars="uid")
tags <- tags[!is.na(value) & value != ""]
tags <- tags[,.(value=strsplit(value, ";")[[1]]), by=.(uid, variable)]
tags[value == "Low_glutamine_/_high_glutamate", value := "Low_glutamine_or_high_glutamate"] # 1 tag with typo
tags <- tags[!(uid %in% sinfo[(removed), uid])]

# Add to raw
raw[, tag := NA_character_]
raw[tags[,.(tag=paste(value,collapse=";")),by=.(uid, variable)], on = .(uid, variable), tag := i.tag]

# On the advice of nightingale, set technical_error measurements to missing (if they aren't already)
raw[tag %like% "error" & !is.na(value), value := NA]

# Remove 4 internal control samples failing quantification
cntrl_fail_quant <- raw[uid %in% sinfo[visit == "Internal Control", uid] & variable == "DHA" & value > 1, uid]
raw <- raw[!(uid %in% cntrl_fail_quant)]
sinfo[uid %in% cntrl_fail_quant, removed := TRUE]
sinfo[uid %in% cntrl_fail_quant, removal_reason := "Failed Quantification"]

# Save data on samples passing qc
if (!dir.exists("data/raw/processed")) dir.create("data/raw/processed")

# Cast to wide and get original column and row order
passqc <- dcast(raw, uid ~ variable, value.var="value")
passqc <- passqc[sinfo, on = .(uid), nomatch=0]
passqc <- passqc[order(uid)]

cntrl_passqc <- passqc[visit == "Internal Control"]
passqc <- passqc[visit != "Internal Control"]

raw <- fread("data/raw/nightingale/UKB-phase1__results.tsv", colClasses=list("character"=1:5, "numeric"=6:254))
cntrl <- fread("data/raw/nightingale/NGH_control_samples.tsv", colClasses=list("character"=1:9, "numeric"=10:258), na.strings=c("\"NA\"", "NA"))

passqc <- passqc[, .SD, .SDcols=names(raw)]
cntrl_passqc <- cntrl_passqc[, .SD, .SDcols=intersect(names(cntrl), names(cntrl_passqc))]

tags <- tags[sinfo, on = .(uid), c("sample_id", "visit") := .(sample_id, visit)]
tags <- tags[order(uid), .(sample_id, visit, biomarker=variable, QC_tag=value)]

fwrite(passqc, sep="\t", quote=FALSE, file="data/raw/processed/UKB_phase1_results.txt")
fwrite(cntrl_passqc, sep="\t", quote=FALSE, file="data/raw/processed/NGH_control_samples.txt")
fwrite(sinfo[visit != "Internal Control"], sep="\t", quote=FALSE, file="data/raw/processed/sample_information.txt")
fwrite(sinfo[visit == "Internal Control"], sep="\t", quote=FALSE, file="data/raw/processed/controls_information.txt")
fwrite(tags, sep="\t", quote=FALSE, file="data/raw/processed/biomarker_QC_tags.txt")


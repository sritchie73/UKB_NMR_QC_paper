library(data.table)
library(ggplot2)
library(ggrastr)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in
# UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Extract albumin
alb <- dat[variable == "Albumin"]

# Get plate outlier limits
outlier_lim <- fread("data/tech_qc/plate_outlier_limits.txt")
outlier_lim <- outlier_lim[variable == "Albumin"]

# Get plate outlier tags
plate_tags <- fread("data/tech_qc/plate_medians_outlier_tagged.txt")
plate_tags <- plate_tags[variable == "Albumin"]

# Compute summary per plate
plate_summary <- alb[!is.na(adj4), .(metric=names(summary(adj4)), value=as.vector(summary(adj4))),
										 by=.(plate_id, plate_MEASURED_DAYS)]
plate_summary <- dcast(plate_summary, plate_id + plate_MEASURED_DAYS ~ metric, value.var="value")
setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

# Order plates by measurement date
plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
plate_summary[, plate_order := 1:.N]

# Add outlier info
plate_summary[plate_tags, on = .(plate_id), outlier := i.outlier]

# Split into three color groups
gg_no_outlier <- plate_summary[outlier == "no"]
gg_high_outlier <- plate_summary[outlier == "high"]
gg_low_outlier <- plate_summary[outlier == "low"]

# Plot plates vs outlier tags:
g1 <- ggplot(plate_summary, aes(x=factor(plate_order))) +
	rasterize(geom_point(data=gg_no_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252")) +
	rasterize(geom_point(data=gg_no_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252")) +
	rasterize(geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252")) +
	rasterize(geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252")) +

	rasterize(geom_point(data=gg_high_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15")) +
	rasterize(geom_point(data=gg_high_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15")) +
	rasterize(geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15")) +
	rasterize(geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15")) +

	rasterize(geom_point(data=gg_low_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c")) +
	rasterize(geom_point(data=gg_low_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c")) +
	rasterize(geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c")) +
	rasterize(geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c")) +

	geom_hline(data=outlier_lim, aes(yintercept=lower_lim), color="#dd1c77", linetype=2, size=0.2) +
	geom_hline(data=outlier_lim, aes(yintercept=upper_lim), color="#dd1c77", linetype=2, size=0.2) +

	ylab("Albumin (g/L)") + xlab("Plate\n(chronological order)") +
	theme_bw() +
	theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				strip.background=element_blank(), strip.text=element_text(size=6),
				panel.grid.major=element_blank(), panel.grid.minor=element_blank()
	)

# Zoom in on interquartile range
g2 <- ggplot(plate_summary, aes(x=factor(plate_order))) +
  rasterize(geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252")) +
  rasterize(geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252")) +

  rasterize(geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15")) +
  rasterize(geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15")) +

  rasterize(geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c")) +
  rasterize(geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c")) +

  geom_hline(data=outlier_lim, aes(yintercept=lower_lim), color="#dd1c77", linetype=2, size=0.2) +
  geom_hline(data=outlier_lim, aes(yintercept=upper_lim), color="#dd1c77", linetype=2, size=0.2) +

  ylab("Albumin (g/L)") + xlab("Plate\n(chronological order)") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# Show clinical biochemistry albumin organised into these same plates
clin <- fread("data/my_curated_phenotypes/biomarkers/output/biomarkers_fixed_limits.txt")
clin <- clin[variable == "alb" & !is.na(value)]
clin[, visit := ifelse(visit == 0, "Main Phase", "Repeat Assessment")]
clin[sinfo[!(sample_removed) & (in_ukb_raw)], on = .(eid=eid_7439, visit), plate_id := plate_id]
clin <- clin[!is.na(plate_id)]
clin[unique(alb[,.(plate_id, plate_MEASURED_DAYS)]), on = .(plate_id), plate_MEASURED_DAYS := plate_MEASURED_DAYS]

# Compute summary per plate
plate_summary <- clin[, .(metric=names(summary(value)), value=as.vector(summary(value))), by=.(plate_id, plate_MEASURED_DAYS)]
plate_summary <- dcast(plate_summary, plate_id + plate_MEASURED_DAYS ~ metric, value.var="value")
setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

# Order plates by measurement date
plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
plate_summary[, plate_order := 1:.N]

# Add outlier info
plate_summary[plate_tags, on = .(plate_id), outlier := i.outlier]

# Split into three color groups
gg_no_outlier <- plate_summary[outlier == "no"]
gg_high_outlier <- plate_summary[outlier == "high"]
gg_low_outlier <- plate_summary[outlier == "low"]

# Zoom in on interquartile range
g3 <- ggplot(plate_summary, aes(x=factor(plate_order))) +
  rasterize(geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252")) +
  rasterize(geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252")) +

  rasterize(geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15")) +
  rasterize(geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15")) +

  rasterize(geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c")) +
  rasterize(geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c")) +

  ylab("Clin.Albumin (g/L)") + xlab("Plate\n(chronological order)") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# Show Nightingale internal controls
cntrl <- fread("data/raw/processed/NGH_control_samples.txt")
cntrl_info <- fread("data/raw/processed/controls_information.txt")

cntrl <- melt(cntrl, id.vars=c("sample_id", "plate_id", "plate_position", "spectrometer"))
cntrl <- cntrl[variable == "Albumin"]
cntrl[cntrl_info, on = .(sample_id, plate_id, plate_position, spectrometer), plate_MEASURED_DAYS := plate_MEASURED_DAYS]
cntrl[plate_tags, on = .(plate_id), outlier := i.outlier]
cntrl[, control_sample := sprintf("%s (Well %s)", sample_id, plate_position)]
cntrl[, control_sample := factor(control_sample, levels=c( # Control samples are paired
        "190404 (Well A01)", "190328 (Well H12)",
        "190425 (Well A01)", "190508 (Well H12)",
        "180829 (Well A01)", "180827 (Well H12)",
        "180830 (Well A01)", "180831 (Well H12)"
      ))]

# Order plates by measurement date
cntrl <- cntrl[order(plate_MEASURED_DAYS)]
cntrl[, plate_order := 1:.N]

# Split into three color groups
gg_no_outlier <- cntrl[outlier == "no"]
gg_high_outlier <- cntrl[outlier == "high"]
gg_low_outlier <- cntrl[outlier == "low"]

g4 <- ggplot(cntrl, aes(x=factor(plate_order), shape=factor(control_sample))) +
	scale_shape_manual(name = "Control Sample (4-pairs used across all plates)",
		values=c( # Control samples are paired
			"190404 (Well A01)"=0, "190328 (Well H12)"=22,
			"190425 (Well A01)"=1, "190508 (Well H12)"=21,
			"180829 (Well A01)"=2, "180827 (Well H12)"=24,
			"180830 (Well A01)"=5, "180831 (Well H12)"=23
		)) +

	rasterize(geom_point(data=gg_no_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#525252")) +
	rasterize(geom_point(data=gg_no_outlier[plate_position == "A01"], aes(y=value), color="#525252", size=0.45, stroke=0.2)) +

	rasterize(geom_point(data=gg_high_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#a50f15")) +
	rasterize(geom_point(data=gg_high_outlier[plate_position == "A01"], aes(y=value), color="#a50f15", size=0.45, stroke=0.2)) +

	rasterize(geom_point(data=gg_low_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#08519c")) +
	rasterize(geom_point(data=gg_low_outlier[plate_position == "A01"], aes(y=value), color="#08519c", size=0.45, stroke=0.2)) +

	ylab(paste0("Albumin (g/L)")) +
	xlab("Plate\n(chronological order)") +
	theme_bw() +
	theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				strip.background=element_blank(), strip.text=element_text(size=6),
				panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
				legend.position="none"
	)

# Collate and output
g <- plot_grid(g1, g2, g3, g4, align="hv", nrow=1)
ggsave(g, width=7.2, height=2, file="paper_output/albumin_outlier_plates.pdf")





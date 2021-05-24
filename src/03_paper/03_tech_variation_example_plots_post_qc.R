library(data.table)
library(ukbnmr)
library(ggplot2)
library(ggrastr)
library(ggthemes)
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

# Plot Alanine levels as a function of plate over time within spectrometer
ala <- dat[variable == "Ala", .(sample_id, visit, adj5)]
ala[sinfo, on = .(sample_id, visit), c("spectrometer", "plate_id", "plate_MEASURED_DATE") := .(spectrometer, plate_id, plate_MEASURED_DATE)]

# Compute per-plate summary
ala <- ala[!is.na(adj5), .(metric=names(summary(adj5)), value=as.vector(summary(adj5))), by=.(spectrometer, plate_id, plate_MEASURED_DATE)]
ala <- dcast(ala, spectrometer + plate_id + plate_MEASURED_DATE ~ metric, value.var="value")
setnames(ala, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

# Order plates by measurement date within spectrometer
ala <- ala[order(plate_MEASURED_DATE)][order(spectrometer)]
ala[, plate_order := 1:.N]

# Plot median and IQR along with min and max values
g1 <- ggplot(ala, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) + 
  rasterize(geom_point(aes(y=Min), shape=1, size=0.2, stroke=0.15, alpha=0.5)) +
  rasterize(geom_point(aes(y=Max), shape=1, size=0.2, stroke=0.15, alpha=0.5)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
  scale_color_calc() +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.3, stroke=0.15)) +
  scale_fill_calc() +
  scale_y_continuous(name="Alanine (mmol/L)", breaks=seq(0.25, 1.25, length=5)) + xlab("Plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7), 
        axis.ticks.x=element_blank(), panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        legend.position="none"
  )

# Now plot without grouping by spectrometer
ala <- ala[order(plate_MEASURED_DATE)]
ala[, plate_order := 1:.N]

g2 <- ggplot(ala, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) + 
  rasterize(geom_point(aes(y=Min), shape=1, size=0.2, stroke=0.15, alpha=0.5)) +
  rasterize(geom_point(aes(y=Max), shape=1, size=0.2, stroke=0.15, alpha=0.5)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
  scale_color_calc() +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.3, stroke=0.15)) +
  scale_fill_calc() +
  scale_y_continuous(name="Alanine (mmol/L)", breaks=seq(0.25, 1.25, length=5)) + xlab("Plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7), 
        axis.ticks.x=element_blank(), panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), axis.title.y=element_blank(),
        legend.position="none"
  )

g <- plot_grid(g1, g2, nrow=1, rel_widths=c(1.2,1))
ggsave(g, width=2.8, height=1.4, file="paper_output/ala_inter_plate_variation_postqc.pdf")

# Plot histidne vs. sample degredation time
his <- dat[variable == "His", .(sample_id, visit, adj5)]
his[sinfo, on = .(sample_id, visit), sample_degredation_time := sample_degredation_time]

g <- ggplot(his, aes(x=sample_degredation_time, y=adj5)) +
	geom_hex() +
	scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
	geom_smooth(color="red", size=0.2, method=MASS::rlm) +
	ylab("Histidine (mmol/L)") +
	xlab("Sample degradation time (hours from sample prep to measurement") +
	theme_bw() +
	theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
				strip.background=element_blank(), strip.text=element_text(size=6),
				panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
				legend.position="bottom", legend.text=element_text(size=6),
        legend.title=element_text(size=7)
	)
ggsave(g, width=1.55, height=2.2, file="paper_output/his_vs_sample_deg_postqc.pdf")


# Plot glycine by well position
gln <- dat[variable == "Gln", .(sample_id, visit, plate_position, plate_row, plate_column, adj5)]

# Compute per-well summary
gln <- gln[!is.na(adj5), .(metric=names(summary(adj5)), value=as.vector(summary(adj5))), by=.(plate_position, plate_row, plate_column)]
gln <- dcast(gln, plate_position + plate_row + plate_column ~ metric, value.var="value")
setnames(gln, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

# Order wells by column within each row
gln <- gln[order(plate_column)][order(plate_row)]
gln[, well_order := 1:.N]

g <- ggplot(gln, aes(x=factor(well_order), color=plate_row, fill=plate_row)) + 
  rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2)) +
  scale_colour_tableau() +
  scale_fill_tableau() +
  ylab("Glycine (mmol/L)") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7), 
        axis.ticks.x=element_blank(), panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        legend.position="none"
  )
ggsave(g, width=1.3, height=1.4, file="paper_output/gln_well_variation_postqc.pdf")

g <- ggplot(gln, aes(x=factor(well_order), color=plate_row, fill=plate_row)) + 
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2)) +
  scale_colour_tableau() +
  scale_fill_tableau() +
  ylab("Glycine (mmol/L)") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7), 
        axis.ticks.x=element_blank(), panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=6), axis.title.y=element_blank(),
        legend.position="none"
  )
ggsave(g, width=1.3, height=1.4, file="paper_output/gln_well_variation_zoom_postqc.pdf")


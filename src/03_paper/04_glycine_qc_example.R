library(data.table)
library(ggplot2)
library(ggrastr)
library(ggthemes)
library(palettetown)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step

# Load technical information and filter to samples in
# UKB adj5 data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Extract glycine
gln <- dat[variable == "Gln" & !is.na(raw)]

# Function to generate plots after each step:
make_plots <- function(column, filename) {
  # Manually rename column of interest for easy extraction
  setnames(gln, column, "value")
  on.exit({ setnames(gln, "value", column) })
  
  # Plot value vs. sample degradation time:
  g1 <- ggplot(gln[!is.na(value)]) +
    aes(x=sample_degredation, y=value) + 
		geom_hex() +
		scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
    geom_smooth(color="red", size=0.2, method=MASS::rlm) + 
    ylab("") + xlab("Hours from prep") + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
          axis.title.x=element_text(size=7), axis.title.y=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="none"
    )

  # Compute summary per well position
	gln_well_summary <- gln[!is.na(value), .(metric=names(summary(value)), value=as.vector(summary(value))), by=.(plate_position, plate_row, plate_column)]
	gln_well_summary <- dcast(gln_well_summary, plate_position + plate_row + plate_column ~ metric, value.var="value")
	setnames(gln_well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  gln_well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by row)
  gln_well_summary <- gln_well_summary[order(plate_column)][order(plate_row)]
  gln_well_summary[, well_order := 1:.N]

  g2 <- ggplot(gln_well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
		rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
		rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
		rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
		rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2)) +
		scale_colour_tableau() +
		scale_fill_tableau() +
		ylab("") + xlab("Well on plate") +
		theme_bw() +
		theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
					axis.ticks.x=element_blank(), panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(), legend.position="none",
					axis.text.y=element_text(size=6), axis.title.y=element_blank()
		)

  # Zoom in to just show interquartile range
	g3 <- ggplot(gln_well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
		rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
		rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2)) +
		scale_colour_tableau() +
		scale_fill_tableau() +
		ylab("") + xlab("Well on plate") +
		theme_bw() +
		theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
					axis.ticks.x=element_blank(), panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(), legend.position="none",
					axis.text.y=element_text(size=6), axis.title.y=element_blank()
		)

  # Show also by plate column
  gln_well_summary <- gln_well_summary[order(plate_row)][order(plate_column)]
  gln_well_summary[, well_order := 1:.N]
	g4 <- ggplot(gln_well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
		rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
		rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2)) +
		scale_colour_poke(pokemon="typhlosion", spread=12) +
		scale_fill_poke(pokemon="typhlosion", spread=12) +
		ylab("") + xlab("Well on plate") +
		theme_bw() +
		theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
					axis.ticks.x=element_blank(), panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(), legend.position="none",
					axis.text.y=element_text(size=6), axis.title.y=element_blank()
		)

  # Compute summary per plate
	gln_plate_summary <- gln[!is.na(value), .(metric=names(summary(value)), value=as.vector(summary(value))), by=.(plate_id, plate_MEASURED_DAYS, spectrometer)]
	gln_plate_summary <- dcast(gln_plate_summary, plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
	setnames(gln_plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date within spectrometer
	gln_plate_summary <- gln_plate_summary[order(plate_MEASURED_DAYS)][order(spectrometer)]
	gln_plate_summary[, plate_order := 1:.N]

	g5 <- ggplot(gln_plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
		rasterize(geom_point(aes(y=Min), shape=1, size=0.2, stroke=0.15, alpha=0.5)) +
		rasterize(geom_point(aes(y=Max), shape=1, size=0.2, stroke=0.15, alpha=0.5)) +
		rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
		scale_color_calc() +
		rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.3, stroke=0.15)) +
		scale_fill_calc() +
		ylab("") + xlab("Plate") +
		theme_bw() +
		theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
					axis.ticks.x=element_blank(), panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(), legend.position="none",
					axis.text.y=element_text(size=6), axis.title.y=element_blank()
    )

  # Zoom in
	g6 <- ggplot(gln_plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
		rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1)) +
		scale_color_calc() +
		rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.3, stroke=0.15)) +
		scale_fill_calc() +
		ylab("") + xlab("Plate") +
		theme_bw() +
		theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
					axis.ticks.x=element_blank(), panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(), legend.position="none",
					axis.text.y=element_text(size=6), axis.title.y=element_blank()
    )

  g <- plot_grid(g1, g2, g3, g4, g5, g6, nrow=1, align="h")
  ggsave(g, width=7.2, height=1, file=filename)
}

make_plots("raw", "paper_output/glycine_techqc_step1_raw.pdf")
make_plots("log_raw", "paper_output/glycine_techqc_step2_log.pdf")
make_plots("log_adj1", "paper_output/glycine_techqc_step3_sample_degredation.pdf")
make_plots("log_adj2", "paper_output/glycine_techqc_step4_plate_row.pdf")
make_plots("log_adj3", "paper_output/glycine_techqc_step5_plate_column.pdf")
make_plots("log_adj4", "paper_output/glycine_techqc_step6_plate_drift.pdf")
make_plots("adj4", "paper_output/glycine_techqc_step7_abs_units.pdf")
make_plots("adj5", "paper_output/glycine_techqc_step8_no_outlier_plates.pdf")


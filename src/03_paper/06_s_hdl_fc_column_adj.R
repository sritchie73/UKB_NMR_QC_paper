library(data.table)
library(ggplot2)
library(ggrastr)
library(palettetown)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step and extract S_HDL_FC
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")
dat <- dat[variable == "S_HDL_FC"]

# Compute summary per well position
well_summary <- dat[!is.na(log_adj2), .(metric=names(summary(log_adj2)), log_adj2=as.vector(summary(log_adj2))), by=.(plate_position, plate_row, plate_column)]
well_summary <- dcast(well_summary, plate_position + plate_row + plate_column ~ metric, value.var="log_adj2")
setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
well_summary[, plate_column := as.factor(plate_column)]

# Plot value vs. plate position (organised by column)
well_summary <- well_summary[order(plate_row)][order(plate_column)]
well_summary[, well_order := 1:.N]

g1 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
	rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
	rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
	rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
	rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_poke(pokemon="typhlosion", spread=12) +
	scale_fill_poke(pokemon="typhlosion", spread=12) +
	ylab("S_HDL_FC") + xlab("Well on plate") +
	theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g2 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
	rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
	rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_poke(pokemon="typhlosion", spread=12) +
	scale_fill_poke(pokemon="typhlosion", spread=12) +
	ylab("S_HDL_FC") + xlab("Well on plate") +
	theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# Show also after removal of inter-column effects
well_summary <- dat[!is.na(log_adj3), .(metric=names(summary(log_adj3)), log_adj3=as.vector(summary(log_adj3))), by=.(plate_position, plate_row, plate_column)]
well_summary <- dcast(well_summary, plate_position + plate_row + plate_column ~ metric, value.var="log_adj3")
setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
well_summary[, plate_column := as.factor(plate_column)]

# Plot value vs. plate position (organised by column)
well_summary <- well_summary[order(plate_row)][order(plate_column)]
well_summary[, well_order := 1:.N]

g3 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
	rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
	rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
	rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
	rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_poke(pokemon="typhlosion", spread=12) +
	scale_fill_poke(pokemon="typhlosion", spread=12) +
	ylab("S_HDL_FC") + xlab("Well on plate") +
	theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g4 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
	rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
	rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_poke(pokemon="typhlosion", spread=12) +
	scale_fill_poke(pokemon="typhlosion", spread=12) +
	ylab("S_HDL_FC") + xlab("Well on plate") +
	theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g <- plot_grid(g1, g2, g3, g4, nrow=2, align="hv")
ggsave(g, width=3.6, height=3.4, file="paper_output/S_HDL_FC_column_variation.pdf")

library(data.table)
library(ggplot2)
library(ggrastr)
library(ggthemes)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step and extract Glycine
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")
dat <- dat[variable == "Gly"]
dat <- dat[!is.na(raw)]

# First, build plots showing the log raw data
well_summary <- dat[, .(metric=names(summary(log_raw)), log_raw=as.vector(summary(log_raw))), by=.(plate_position, plate_row, plate_column)]
well_summary <- dcast(well_summary, plate_position + plate_row + plate_column ~ metric, value.var="log_raw")
setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
well_summary[, plate_column := as.factor(plate_column)]
well_summary <- well_summary[order(plate_column)][order(plate_row)]
well_summary[, well_order := 1:.N]

g1 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
  rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_tableau() +
	scale_fill_tableau() +
  ylab("Glycine (log mmol/L)") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g2 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_tableau() +
	scale_fill_tableau() +
  ylab("Glycine") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# Show results of robust linear regression to remove sample degradation time
well_summary <- dat[, .(metric=names(summary(log_adj1)), log_adj1=as.vector(summary(log_adj1))), by=.(plate_position, plate_row, plate_column)]
well_summary <- dcast(well_summary, plate_position + plate_row + plate_column ~ metric, value.var="log_adj1")
setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
well_summary[, plate_column := as.factor(plate_column)]
well_summary <- well_summary[order(plate_column)][order(plate_row)]
well_summary[, well_order := 1:.N]

g3 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
  rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_tableau() +
	scale_fill_tableau() +
  ylab("Glycine (residuals)") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g4 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_tableau() +
	scale_fill_tableau() +
  ylab("Glycine") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )


# Now show what happens when we use linear regression instead 
dat[, log_adj1_lm := lm(log_raw ~ sample_degredation)$residuals]
well_summary <- dat[, .(metric=names(summary(log_adj1_lm)), log_adj1_lm=as.vector(summary(log_adj1_lm))), by=.(plate_position, plate_row, plate_column)]
well_summary <- dcast(well_summary, plate_position + plate_row + plate_column ~ metric, value.var="log_adj1_lm")
setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
well_summary[, plate_column := as.factor(plate_column)]
well_summary <- well_summary[order(plate_column)][order(plate_row)]
well_summary[, well_order := 1:.N]

g5 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
  rasterize(geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_tableau() +
	scale_fill_tableau() +
  ylab("Glycine (residuals)") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g6 <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
  rasterize(geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.25)) +
  rasterize(geom_point(aes(y=Median), color="black", shape=21, size=0.6, stroke=0.3)) +
	scale_colour_tableau() +
	scale_fill_tableau() +
  ylab("Glycine") + xlab("Well on plate") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        axis.ticks.x=element_blank(),  legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

# build plot
g <- plot_grid(g1, g2, g3, g4, NULL, NULL, g5, g6, nrow=2, align="hv")
ggsave(g, width=5.7, height=3.6, file="paper_output/glycine_lm_vs_rlm.pdf")

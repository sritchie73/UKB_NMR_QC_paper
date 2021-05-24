library(data.table)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(ochRe) # remotes::install_github("ropenscilabs/ochRe")

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load technical covariates
tech_info <- fread("data/tech_qc/sample_information.txt")

# Extract plate information
plate_info <- tech_info[!(sample_removed), .(samples=.N), by=.(plate_id, plate_MEASURED_DATE, plate_measured_bin, spectrometer)]

# Plot number of plates by date, coloured by bin
plate_by_date <- plate_info[,.N,by=.(plate_MEASURED_DATE, plate_measured_bin, spectrometer)]

g1 <- ggplot(plate_by_date) +
  aes(x = plate_MEASURED_DATE, y = N, fill = factor(plate_measured_bin)) +
  rasterize(geom_col(width=0.6)) +
  facet_wrap(~ spectrometer, nrow=3) +
  scale_fill_ochre(name="Drift over time bin", palette = "olsen_seq") +
  xlab("Plate measurement date") +
  ylab("Number of plates") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7), 
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.text=element_text(size=7), strip.background=element_blank(),
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(size=0.3), panel.grid.minor.y=element_line(size=0.2),
        legend.position="bottom", legend.title=element_text(size=7), legend.text=element_text(size=6)
  )

# How many samples in each group?
samples_by_bin <- plate_info[, .(samples=sum(samples)), by=.(spectrometer, plate_measured_bin)]

g2 <- ggplot(samples_by_bin) +
  aes(x = factor(plate_measured_bin), y = samples, fill = factor(plate_measured_bin)) +
  geom_col(width=0.6) +
  facet_wrap(~ spectrometer, nrow=3) +
  scale_fill_ochre(name="Drift over time bin", palette = "olsen_seq") +
  xlab("Drift over time bin") +
  ylab("Number of samples") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7), 
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.text=element_text(size=7), strip.background=element_blank(),
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(size=0.3), panel.grid.minor.y=element_line(size=0.2),
        legend.position="bottom", legend.title=element_text(size=7), legend.text=element_text(size=6)
  )

g <- plot_grid(g1, g2, nrow=1, align="hv")
ggsave(g, width=7.2, height=3.6, file="paper_output/plate_bins_drift_over_time.pdf")


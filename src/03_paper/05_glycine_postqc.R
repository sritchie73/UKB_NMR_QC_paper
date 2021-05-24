library(data.table)
library(ggplot2)
library(ggrastr)
library(ggthemes)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Filter to glycine
gln <- dat[variable == "Gln" & !is.na(raw)]

# Compare raw vs. postqc data
g1 <- ggplot(gln) +
  aes(x = raw, y = adj5) +
  rasterize(geom_point(shape=21, size=0.45, stroke=0.2, color="black", fill="white")) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_x_continuous(name="Raw glycine (mmol/L)", limits=c(0, 1.5)) +
  scale_y_continuous(name="Post-QC glycine (mmol/L)", limits=c(0, 1.5)) +
  theme_bw() +
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
				axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  )

# Now on log scale
g2 <- ggplot(gln) +
  aes(x = raw, y = adj5) +
  rasterize(geom_point(shape=21, size=0.45, stroke=0.2, color="black", fill="white")) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_x_log10(name="Raw (log scale)") +
  scale_y_log10(name="Post-QC (log scale)") +
  theme_bw() +
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
				axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  )

# Overlay distributions
g3 <- ggplot(gln) +
  geom_density(aes(x=raw), color="#e08214", size=0.4, trim=TRUE) +
  geom_density(aes(x=adj5), color="#542788", size=0.4, trim=TRUE) +
  scale_x_continuous(name="Glycine (mmol/L)", limits=c(0, 1.5), breaks=c(0, 0.5, 1, 1.5)) + 
  scale_y_continuous(name="Density") +
  theme_bw() +
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
				axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  )

# Now on log scale
g4 <- ggplot(gln) +
  geom_density(aes(x=raw), color="#e08214", size=0.4, trim=TRUE) +
  geom_density(aes(x=adj5), color="#542788", size=0.4, trim=TRUE) +
  scale_x_log10(name="Glycine (log scale)") + 
  scale_y_continuous(name="Density") +
  theme_bw() +
	theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
				axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.2), panel.grid.minor=element_line(size=0.1)
  )

g <- plot_grid(g1, g2, g3, g4, align="hv", nrow=1)
ggsave(g, height=1, width=4.4, file="paper_output/raw_vs_postqc_glycine.pdf")

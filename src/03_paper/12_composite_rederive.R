library(data.table)
library(ukbnmr)
library(ggplot2)
library(ggrastr)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in
# UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Show analytic validity of rederiving biomarkers
gg_dt <- dat[variable %in% nmr_info[Type != "Non-derived" & (Nightingale), Biomarker]]
gg_dt <- gg_dt[!is.na(raw) & !is.na(raw_rederived)]
diff <- gg_dt[, .(Pearson = cor(raw, raw_rederived), Spearman = cor(raw, raw_rederived, method="spearman")), by=.(variable)]
fwrite(diff, sep="\t", quote=FALSE, file="diagnostic_plots/analytical_validity_correlations.txt")

# Filter to top 6 most different composite biomarkers
diff <- diff[variable %in% nmr_info[Type == "Composite", Biomarker]]
diff <- diff[order(Pearson)][1:6]
gg_dt <- gg_dt[diff[,.(variable)], on = .(variable)]
gg_dt[, variable := factor(variable, levels=diff$variable)]

g1 <- ggplot(gg_dt) + 
  aes(x=raw, y=raw_rederived) +
  rasterize(geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2)) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.2) +  
  scale_x_continuous(name="Biomarker in raw data (mmol/L or %)") +
  scale_y_continuous(name="Rederived from parts") +
  facet_wrap(~ variable, nrow=1, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g2 <- ggplot(gg_dt) + 
  aes(x=raw, y=raw - raw_rederived) +
  rasterize(geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2)) +
  geom_hline(yintercept=0, linetype=2, color="red", size=0.2) +  
  scale_x_continuous(name="Biomarker in raw data (mmol/L or %)") +
  scale_y_continuous(name="Raw - Rederived") +
  facet_wrap(~ variable, nrow=1, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

g <- plot_grid(g1, g2, align="hv", nrow=2)
ggsave(g, width=7.2, height=3, file="paper_output/composite_rederive_validity.pdf")



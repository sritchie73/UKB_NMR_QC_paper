library(data.table)
library(ukbnmr) # remotes::install_github("sritchie73/ukbnmr", ref="development")
library(MASS)
library(ggplot2)
library(hexbin)
library(ggthemes)
library(palettetown)

## Warning: generating all the diagnostic plots takes ~20 hours.
ggsave <- function(file, ...) {
  if (!file.exists(file)) {
    ggplot2::ggsave(file=file, ...)
  } else {
    message(sprintf("%s already exists, skipping", file))
  }
}

# Create output directories
if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")

# Load data from previous script
adj <- fread("data/tech_qc/multistep_adjusted_values.txt")
sinfo <- fread("data/tech_qc/sample_information.txt")
flags <- fread("data/tech_qc/biomarker_QC_flags.txt")
log_offset <- fread("data/tech_qc/log_offset.txt")
plate_medians <- fread("data/tech_qc/plate_medians_outlier_tagged.txt")
outlier_lim <- fread("data/tech_qc/plate_outlier_limits.txt")
plate_medians_comb <- fread("data/tech_qc/plate_medians_outlier_tagged_comb_regression.txt")
outlier_lim_comb <- fread("data/tech_qc/plate_outlier_limits_comb_regression.txt")
plate_medians_comb_strata <- fread("data/tech_qc/plate_medians_outlier_tagged_comb_strata_regression.txt")
outlier_lim_comb_strata <- fread("data/tech_qc/plate_outlier_limits_comb_strata_regression.txt")
plate_medians_ts <- fread("data/tech_qc/plate_medians_outlier_tagged_ts_regression.txt")
outlier_lim_ts <- fread("data/tech_qc/plate_outlier_limits_ts_regression.txt")
plate_medians_4321 <- fread("data/tech_qc/plate_medians_outlier_tagged_4321_regression.txt")
outlier_lim_4321 <- fread("data/tech_qc/plate_outlier_limits_4321_regression.txt")
plate_medians_3214 <- fread("data/tech_qc/plate_medians_outlier_tagged_3214_regression.txt")
outlier_lim_3214 <- fread("data/tech_qc/plate_outlier_limits_3214_regression.txt")

# Make sure biomarkers always shown in alphabetical order
adj <- adj[order(tolower(variable))]
adj[, variable := factor(variable, levels=unique(variable))]

# Analytical validity of biomarker rederivation
gg_dt <- adj[variable %in% ukbnmr::nmr_info[Type != "Non-derived" & (Nightingale), Biomarker]]
gg_dt <- gg_dt[!is.na(raw) & !is.na(raw_rederived)]

g <- ggplot(gg_dt) +
  aes(x=raw, y=raw_rederived) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_y_continuous(name="Derived and composite biomarkers computed from raw parts") +
  scale_x_continuous(name="Biomarker in raw data") +
  facet_wrap(~ variable, ncol=12, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=16.8, height=16.8, units="in", file="diagnostic_plots/rederivation_validity.png")

g <- ggplot(gg_dt) +
  aes(x=raw, y=raw - raw_rederived) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_hline(yintercept=0, linetype=2, color="red", size=0.4) +
  scale_y_continuous(name="Difference between raw and derived") +
  scale_x_continuous(name="Biomarker in raw data") +
  facet_wrap(~ variable, ncol=12, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=16.8, height=16.8, units="in", file="diagnostic_plots/rederivation_validity_diff.png")

# Rescale residuals to match log raw data for comparison
rescale <- function(log_adj, log_raw) {
  log_adj + as.vector(coef(rlm(log_raw ~ 1)))[1]
}

# Put values back on absolute units scale
adj[!is.na(log_raw), log_adj1 := rescale(log_adj1, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj2 := rescale(log_adj2, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj3 := rescale(log_adj3, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj4 := rescale(log_adj4, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj_comb := rescale(log_adj_comb, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj_comb_strata := rescale(log_adj_comb_strata, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj_ts := rescale(log_adj_ts, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj4_4321 := rescale(log_adj4_4321, log_raw), by=.(variable)]
adj[!is.na(log_raw), log_adj4_3214 := rescale(log_adj4_3214, log_raw), by=.(variable)]

# Need to create version for adj5 (dropping outlier plates)
adj[!is.na(log_adj4) & !is.na(adj5), log_adj5 := log_adj4]
adj[!is.na(log_adj_comb) & !is.na(adj5_comb), log_adj5_comb := log_adj_comb]
adj[!is.na(log_adj_comb_strata) & !is.na(adj5_comb_strata), log_adj5_comb_strata := log_adj_comb_strata]
adj[!is.na(log_adj_ts) & !is.na(adj5_ts), log_adj5_ts := log_adj_ts]
adj[!is.na(log_adj4_4321) & !is.na(adj5_4321), log_adj5_4321 := log_adj4_4321]
adj[!is.na(log_adj4_3214) & !is.na(adj5_3214), log_adj5_3214 := log_adj4_3214]

# Compare concentrations before and after adjustment
conc_comp <- function(x_col, y_col, xlab, ylab, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, c(x_col, y_col), c("x_col", "y_col"))
  on.exit(setnames(adj, c("x_col", "y_col"), c(x_col, y_col)))
  n_bio <- adj[!is.na(x_col) & !is.na(y_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  g <- ggplot(adj[!is.na(x_col) & !is.na(y_col)], aes(x = x_col, y = y_col)) + 
    geom_point(shape=21, size=0.45, stroke=0.2, color="black", fill="white") +
    geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_x_continuous(name=xlab) +
    scale_y_continuous(name=ylab) +
		theme_bw() +
		theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
					strip.background=element_blank(), strip.text=element_text(size=6),
					panel.grid.major=element_blank(), panel.grid.minor=element_blank()
    )

  ggsave(g, width=width, height=width, units="in", file=fname)
}

conc_comp("raw", "adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_concentrations.png")
conc_comp("raw", "adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_concentrations.png")
conc_comp("raw", "adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_concentrations.png")
conc_comp("raw", "adj4", "Raw biomarker concentrations", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/raw_vs_adj4_concentrations.png")
conc_comp("raw", "adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_concentrations.png")
conc_comp("adj1", "adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_concentrations.png")
conc_comp("adj2", "adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_concentrations.png")
conc_comp("adj3", "adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_concentrations.png")
conc_comp("adj4", "adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_concentrations.png")

conc_comp("log_raw", "log_adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_log_concentrations.png")
conc_comp("log_raw", "log_adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_log_concentrations.png")
conc_comp("log_raw", "log_adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_log_concentrations.png")
conc_comp("log_raw", "log_adj4", "Raw biomarker concentrations", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/raw_vs_adj4_log_concentrations.png")
conc_comp("log_raw", "log_adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_log_concentrations.png")
conc_comp("log_adj1", "log_adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_log_concentrations.png")
conc_comp("log_adj2", "log_adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_log_concentrations.png")
conc_comp("log_adj3", "log_adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_log_concentrations.png")
conc_comp("log_adj4", "log_adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_log_concentrations.png")

conc_comp("raw", "adj_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/raw_vs_adj_comb_concentrations.png")
conc_comp("raw", "adj_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/raw_vs_adj_comb_strata_concentrations.png")
conc_comp("raw", "adj_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/raw_vs_adj_ts_concentrations.png")
conc_comp("raw", "adj4_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/raw_vs_adj4_4321_concentrations.png")
conc_comp("raw", "adj4_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/raw_vs_adj4_3214_concentrations.png")

conc_comp("raw", "adj5_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_concentrations.png")
conc_comp("raw", "adj5_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_strata_concentrations.png")
conc_comp("raw", "adj5_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_ts_concentrations.png")
conc_comp("raw", "adj5_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_4321_concentrations.png")
conc_comp("raw", "adj5_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_3214_concentrations.png")

conc_comp("adj4", "adj_comb", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step adjustment", "diagnostic_plots/adj4_vs_adj_comb_concentrations.png")
conc_comp("adj4", "adj_comb_strata", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step spectrometer stratified adjustment", "diagnostic_plots/adj4_vs_adj_comb_strata_concentrations.png")
conc_comp("adj4", "adj_ts", "(Post-QC) Multi-step adjustment", "(Post-QC) Two-step adjustment", "diagnostic_plots/adj4_vs_adj_ts_concentrations.png")
conc_comp("adj4", "adj4_4321", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_vs_adj4_4321_concentrations.png")
conc_comp("adj4", "adj4_3214", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_vs_adj4_3214_concentrations.png")

conc_comp("adj5", "adj5_comb", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step adjustment, outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_comb_concentrations.png")
conc_comp("adj5", "adj5_comb_strata", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step spectrometer stratified adjustment, outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_comb_strata_concentrations.png")
conc_comp("adj5", "adj5_ts", "(Post-QC) Two-step adjustment, outlier plates excluded", "(Post-QC) Single step adjustment, outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_ts_concentrations.png")
conc_comp("adj5", "adj5_4321", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_4321_concentrations.png")
conc_comp("adj5", "adj5_3214", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_3214_concentrations.png")

conc_comp("log_raw", "log_adj_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/raw_vs_adj_comb_log_concentrations.png")
conc_comp("log_raw", "log_adj_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/raw_vs_adj_comb_strata_log_concentrations.png")
conc_comp("log_raw", "log_adj_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/raw_vs_adj_ts_log_concentrations.png")
conc_comp("log_raw", "log_adj4_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/raw_vs_adj4_4321_log_concentrations.png")
conc_comp("log_raw", "log_adj4_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/raw_vs_adj4_3214_log_concentrations.png")

conc_comp("log_raw", "log_adj5_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_log_concentrations.png")
conc_comp("log_raw", "log_adj5_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_strata_log_concentrations.png")
conc_comp("log_raw", "log_adj5_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_ts_log_concentrations.png")
conc_comp("log_raw", "log_adj5_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_4321_log_concentrations.png")
conc_comp("log_raw", "log_adj5_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_3214_log_concentrations.png")

conc_comp("log_adj4", "log_adj_comb", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step adjustment", "diagnostic_plots/adj4_vs_adj_comb_log_concentrations.png")
conc_comp("log_adj4", "log_adj_comb_strata", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step spectrometer stratified adjustment", "diagnostic_plots/adj4_vs_adj_comb_strata_log_concentrations.png")
conc_comp("log_adj4", "log_adj_ts", "(Post-QC) Multi-step adjustment", "(Post-QC) Two-stage adjustment", "diagnostic_plots/adj4_vs_adj_ts_log_concentrations.png")
conc_comp("log_adj4", "log_adj4_4321", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_vs_adj4_4321_log_concentrations.png")
conc_comp("log_adj4", "log_adj4_3214", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_vs_adj4_3214_log_concentrations.png")

conc_comp("log_adj5", "log_adj5_comb", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step adjustment, outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_comb_log_concentrations.png")
conc_comp("log_adj5", "log_adj5_comb_strata", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step spectrometer stratified adjustment, outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_comb_strata_log_concentrations.png")
conc_comp("log_adj5", "log_adj5_ts", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Two-stage adjustment, outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_ts_log_concentrations.png")
conc_comp("log_adj5", "log_adj5_4321", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_4321_log_concentrations.png")
conc_comp("log_adj5", "log_adj5_3214", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_vs_adj5_3214_log_concentrations.png")

# Compare distributions before and after adjustment
dens_comp <- function(raw_col, adj_col, raw_lab, adj_lab, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  gg_dt <- melt(adj, id.vars=c("sample_id", "visit", "variable"), measure.vars=c(raw_col, adj_col), variable.name="type")
  gg_dt <- gg_dt[!is.na(value)]
  n_bio <- gg_dt[, length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  g <- ggplot(gg_dt, aes(x = value, color=type)) + 
    geom_density(trim=TRUE, size=0.4) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_manual(name="", labels=structure(c(raw_lab, adj_lab), names=c(raw_col, adj_col)), 
                        values=structure(c("#e08214", "#542788"), names=c(raw_col, adj_col))) +
    xlab("Biomarker concentrations") +
    ylab("Density") +
		theme_bw() +
		theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
					strip.background=element_blank(), strip.text=element_text(size=6),
					panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    )

  ggsave(g, width=width, height=width, units="in", file=fname)
}

dens_comp("raw", "adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_densities.png")
dens_comp("raw", "adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_densities.png")
dens_comp("raw", "adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_densities.png")
dens_comp("raw", "adj4", "Raw biomarker concentrations", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/raw_vs_adj4_densities.png")
dens_comp("raw", "adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_densities.png")
dens_comp("adj1", "adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_densities.png")
dens_comp("adj2", "adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_densities.png")
dens_comp("adj3", "adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_densities.png")
dens_comp("adj4", "adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_densities.png")

dens_comp("log_raw", "log_adj1", "Raw biomarker concentrations", "Adjusted for sample degradation time", "diagnostic_plots/raw_vs_adj1_log_densities.png")
dens_comp("log_raw", "log_adj2", "Raw biomarker concentrations", "Adjusted for plate row", "diagnostic_plots/raw_vs_adj2_log_densities.png")
dens_comp("log_raw", "log_adj3", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj3_log_densities.png")
dens_comp("log_raw", "log_adj4", "Raw biomarker concentrations", "Adjusted for plate column", "diagnostic_plots/raw_vs_adj4_log_densities.png")
dens_comp("log_raw", "log_adj5", "Raw biomarker concentrations", "(Post-QC) Excluding outlier plates", "diagnostic_plots/raw_vs_adj5_log_densities.png")
dens_comp("log_adj1", "log_adj2", "Adjusted for sample degradation time", "Adjusted for plate row", "diagnostic_plots/adj1_vs_adj2_log_densities.png")
dens_comp("log_adj2", "log_adj3", "Adjusted for plate row", "Adjusted for plate column", "diagnostic_plots/adj2_vs_adj3_log_densities.png")
dens_comp("log_adj3", "log_adj4", "Adjusted for plate column", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj3_vs_adj4_log_densities.png")
dens_comp("log_adj4", "log_adj5", "(Post-QC) Adjusted for spectrometer drift over time", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj4_vs_adj5_log_densities.png")

dens_comp("raw", "adj_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/raw_vs_adj_comb_densities.png")
dens_comp("raw", "adj_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/raw_vs_adj_comb_strata_densities.png")
dens_comp("raw", "adj_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/raw_vs_adj_ts_densities.png")
dens_comp("raw", "adj4_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/raw_vs_adj4_4321_densities.png")
dens_comp("raw", "adj4_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/raw_vs_adj4_3214_densities.png")

dens_comp("raw", "adj5_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_densities.png")
dens_comp("raw", "adj5_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_strata_densities.png")
dens_comp("raw", "adj5_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_ts_densities.png")
dens_comp("raw", "adj5_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_4321_densities.png")
dens_comp("raw", "adj5_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_3214_densities.png")

dens_comp("adj4", "adj_comb", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step adjustment", "diagnostic_plots/adj4_vs_adj_comb_densities.png")
dens_comp("adj4", "adj_comb_strata", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step spectrometer stratified adjustment", "diagnostic_plots/adj4_vs_adj_comb_strata_densities.png")
dens_comp("adj4", "adj_ts", "(Post-QC) Multi-step adjustment", "(Post-QC) Two-stage adjustment", "diagnostic_plots/adj4_vs_adj_ts_densities.png")
dens_comp("adj4", "adj4_4321", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_vs_adj4_4321_densities.png")
dens_comp("adj4", "adj4_3214", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_vs_adj4_3214_densities.png")

dens_comp("adj5", "adj5_comb", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step adjustment, outlier plates excluded", "diagnostic_plots/adj4_vs_adj_comb_densities.png")
dens_comp("adj5", "adj5_comb_strata", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step spectrometer stratified adjustment, outlier plates excluded", "diagnostic_plots/adj4_vs_adj_comb_strata_densities.png")
dens_comp("adj5", "adj5_ts", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Two-stage adjustment, outlier plates excluded", "diagnostic_plots/adj4_vs_adj_ts_densities.png")
dens_comp("adj5", "adj5_4321", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_4321_densities.png")
dens_comp("adj5", "adj5_3214", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_3214_densities.png")

dens_comp("log_raw", "log_adj_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/raw_vs_adj_comb_log_densities.png")
dens_comp("log_raw", "log_adj_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/raw_vs_adj_comb_strata_log_densities.png")
dens_comp("log_raw", "log_adj_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/raw_vs_adj_ts_log_densities.png")
dens_comp("log_raw", "log_adj4_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/raw_vs_adj4_4321_log_densities.png")
dens_comp("log_raw", "log_adj4_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/raw_vs_adj4_3214_log_densities.png")

dens_comp("log_raw", "log_adj5_comb", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_log_densities.png")
dens_comp("log_raw", "log_adj5_comb_strata", "Raw biomarker concentrations", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_comb_strata_log_densities.png")
dens_comp("log_raw", "log_adj5_ts", "Raw biomarker concentrations", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/raw_vs_adj5_ts_log_densities.png")
dens_comp("log_raw", "log_adj5_4321", "Raw biomarker concentrations", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_4321_log_densities.png")
dens_comp("log_raw", "log_adj5_3214", "Raw biomarker concentrations", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_3214_log_densities.png")

dens_comp("log_adj4", "log_adj_comb", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step adjustment", "diagnostic_plots/adj4_vs_adj_comb_log_densities.png")
dens_comp("log_adj4", "log_adj_comb_strata", "(Post-QC) Multi-step adjustment", "(Post-QC) Single step spectrometer stratified adjustment", "diagnostic_plots/adj4_vs_adj_comb_strata_log_densities.png")
dens_comp("log_adj4", "log_adj_ts", "(Post-QC) Multi-step adjustment", "(Post-QC) Two-stage adjustment", "diagnostic_plots/adj4_vs_adj_ts_log_densities.png")
dens_comp("log_adj4", "log_adj4_4321", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_vs_adj4_4321_log_densities.png")
dens_comp("log_adj4", "log_adj4_3214", "(Post-QC) Multi-step adjustment", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_vs_adj4_3214_log_densities.png")

dens_comp("log_adj5", "log_adj5_comb", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step adjustment, outlier plates excluded", "diagnostic_plots/adj4_vs_adj_comb_log_densities.png")
dens_comp("log_adj5", "log_adj5_comb_strata", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Single step spectrometer stratified adjustment, outlier plates excluded", "diagnostic_plots/adj4_vs_adj_comb_strata_log_densities.png")
dens_comp("log_adj5", "log_adj5_ts", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Two-stage adjustment, outlier plates excluded", "diagnostic_plots/adj4_vs_adj_ts_log_densities.png")
dens_comp("log_adj5", "log_adj5_4321", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_4321_log_densities.png")
dens_comp("log_adj5", "log_adj5_3214", "(Post-QC) Multi-step adjustment, outlier plates excluded", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/raw_vs_adj5_3214_log_densities.png")

# How does sample degradation time change at each step?
sample_deg <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  g <- ggplot(adj[!is.na(conc_col) & is.finite(conc_col)], aes(x=sample_degredation, y=conc_col)) +
    geom_hex() +
    scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10") +
    geom_smooth(color="red", size=0.2, method=MASS::rlm) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(conc_name) + 
    scale_x_log10("Sample degradation time (hours from sample prep to measurement") +
    theme_bw() +
		theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
					strip.background=element_blank(), strip.text=element_text(size=6),
					panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom"
    )

  ggsave(g, width=width, height=width, units="in", file=fname)
}

sample_deg("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_sample_degradation_time.png")
sample_deg("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_sample_degradation_time.png")
sample_deg("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_sample_degradation_time.png")
sample_deg("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_sample_degradation_time.png")
sample_deg("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_sample_degradation_time.png")
sample_deg("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_sample_degradation_time.png")

sample_deg("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_sample_degradation_time.png")
sample_deg("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_sample_degradation_time.png")
sample_deg("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_sample_degradation_time.png")
sample_deg("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_sample_degradation_time.png")
sample_deg("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_sample_degradation_time.png")
sample_deg("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_sample_degradation_time.png")

sample_deg("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_sample_degradation_time.png")
sample_deg("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_sample_degradation_time.png")
sample_deg("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_sample_degradation_time.png")
sample_deg("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_sample_degradation_time.png")
sample_deg("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_sample_degradation_time.png")

sample_deg("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_sample_degradation_time.png")
sample_deg("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_sample_degradation_time.png")
sample_deg("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_sample_degradation_time.png")
sample_deg("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_sample_degradation_time.png")
sample_deg("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_sample_degradation_time.png")

sample_deg("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_sample_degradation_time.png")
sample_deg("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_sample_degradation_time.png")
sample_deg("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_sample_degradation_time.png")
sample_deg("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_sample_degradation_time.png")
sample_deg("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_sample_degradation_time.png")

sample_deg("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_sample_degradation_time.png")
sample_deg("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_sample_degradation_time.png")
sample_deg("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_sample_degradation_time.png")
sample_deg("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_sample_degradation_time.png")
sample_deg("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_sample_degradation_time.png")

# Show concentrations vs. plate row
plate_row <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by row)
  well_summary <- well_summary[order(plate_column)][order(plate_row)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_tableau(name="Row") +
    scale_fill_tableau(name="Row") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by column (1-12) within each row (A-H))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_row("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_row.png")
plate_row("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_row.png")
plate_row("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_row.png")
plate_row("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_row.png")
plate_row("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_row.png")
plate_row("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_row.png")

plate_row("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_row.png")
plate_row("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_row.png")
plate_row("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_row.png")
plate_row("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_row.png")
plate_row("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_row.png")
plate_row("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_row.png")

plate_row("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_plate_row.png")
plate_row("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_plate_row.png")
plate_row("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_plate_row.png")
plate_row("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_plate_row.png")
plate_row("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_plate_row.png")
plate_row("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_plate_row.png")

plate_row("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_plate_row.png")
plate_row("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_plate_row.png")
plate_row("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_plate_row.png")
plate_row("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_plate_row.png")
plate_row("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_plate_row.png")
plate_row("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_plate_row.png")

plate_row("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_plate_row.png")
plate_row("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_plate_row.png")
plate_row("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_plate_row.png")
plate_row("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_plate_row.png")

plate_row("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_plate_row.png")
plate_row("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_plate_row.png")
plate_row("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_plate_row.png")
plate_row("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_plate_row.png")

# Zoom in on interquartile range
plate_row_iqr <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by row)
  well_summary <- well_summary[order(plate_column)][order(plate_row)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_row, fill=plate_row)) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_tableau(name="Row") +
    scale_fill_tableau(name="Row") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by column (1-12) within each row (A-H))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_row_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_row_iqr_zoom.png")

plate_row_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_row_iqr_zoom.png")

plate_row_iqr("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_plate_row_iqr_zoom.png")
plate_row_iqr("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_plate_row_iqr_zoom.png")

plate_row_iqr("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_plate_row_iqr_zoom.png")
plate_row_iqr("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_plate_row_iqr_zoom.png")

plate_row_iqr("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_plate_row_iqr.png")
plate_row_iqr("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_plate_row_iqr.png")
plate_row_iqr("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_plate_row_iqr.png")
plate_row_iqr("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_plate_row_iqr.png")

plate_row_iqr("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_plate_row_iqr.png")
plate_row_iqr("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_plate_row_iqr.png")
plate_row_iqr("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_plate_row_iqr.png")
plate_row_iqr("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_plate_row_iqr.png")

# Show concentrations vs. plate column
plate_column <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by column)
  well_summary <- well_summary[order(plate_row)][order(plate_column)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_poke(name="Column", pokemon="typhlosion", spread=12) +
    scale_fill_poke(name="Column", pokemon="typhlosion", spread=12) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by row (A-H) within each column (1-12))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_column("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_column.png")
plate_column("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_column.png")
plate_column("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_column.png")
plate_column("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_column.png")
plate_column("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_column.png")
plate_column("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_column.png")

plate_column("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_column.png")
plate_column("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_column.png")
plate_column("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_column.png")
plate_column("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_column.png")
plate_column("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_column.png")
plate_column("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_column.png")

plate_column("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_plate_column.png")
plate_column("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_plate_column.png")
plate_column("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_plate_column.png")
plate_column("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_plate_column.png")
plate_column("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_plate_column.png")
plate_column("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_plate_column.png")

plate_column("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_plate_column.png")
plate_column("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_plate_column.png")
plate_column("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_plate_column.png")
plate_column("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_plate_column.png")
plate_column("log_adj5_comb_strata", "(Post-QC) All adjustments in a spectrometer stratified single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_plate_column.png")
plate_column("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_plate_column.png")

plate_column("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_plate_column.png")
plate_column("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_plate_column.png")
plate_column("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_plate_column.png")
plate_column("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_plate_column.png")

plate_column("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_plate_column.png")
plate_column("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_plate_column.png")
plate_column("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_plate_column.png")
plate_column("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_plate_column.png")

# Zoom in on interquartile range
plate_column_iqr <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per well position
  well_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                      .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                      by=.(plate_position, plate_row, plate_column, variable)]
  well_summary <- dcast(well_summary, variable + plate_position + plate_row + plate_column ~ metric, value.var="value")
  setnames(well_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  well_summary[, plate_column := as.factor(plate_column)]

  # Plot value vs. plate position (organised by column)
  well_summary <- well_summary[order(plate_row)][order(plate_column)][order(variable)]
  well_summary[, well_order := 1:.N]

  g <- ggplot(well_summary, aes(x=factor(well_order), color=plate_column, fill=plate_column)) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_colour_poke(name="Column", pokemon="typhlosion", spread=12) +
    scale_fill_poke(name="Column", pokemon="typhlosion", spread=12) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Position on 96-well plate\n(ordered by row (A-H) within each column (1-12))") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_column_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_column_iqr_zoom.png")

plate_column_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_column_iqr_zoom.png")

plate_column_iqr("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_plate_column_iqr_zoom.png")
plate_column_iqr("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_plate_column_iqr_zoom.png")

plate_column_iqr("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_plate_column_iqr_zoom.png")
plate_column_iqr("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_plate_column_iqr_zoom.png")

plate_column_iqr("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_plate_column_iqr.png")
plate_column_iqr("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_plate_column_iqr.png")
plate_column_iqr("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_plate_column_iqr.png")
plate_column_iqr("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_plate_column_iqr.png")

plate_column_iqr("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_plate_column_iqr.png")
plate_column_iqr("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_plate_column_iqr.png")
plate_column_iqr("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_plate_column_iqr.png")
plate_column_iqr("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_plate_column_iqr.png")

# Show concentrations vs. inter-plate variation (drift over time) within each spectrometer
spec_drift <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date within spectrometer
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)][order(spectrometer)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order within spectrometer)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

spec_drift("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_spectrometer_drift.png")
spec_drift("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_spectrometer_drift.png")
spec_drift("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_spectrometer_drift.png")
spec_drift("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_spectrometer_drift.png")
spec_drift("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_spectrometer_drift.png")
spec_drift("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_spectrometer_drift.png")

spec_drift("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_spectrometer_drift.png")
spec_drift("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_spectrometer_drift.png")
spec_drift("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_spectrometer_drift.png")
spec_drift("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_spectrometer_drift.png")
spec_drift("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_spectrometer_drift.png")
spec_drift("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_spectrometer_drift.png")

spec_drift("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_spectrometer_drift.png")
spec_drift("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_spectrometer_drift.png")
spec_drift("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_spectrometer_drift.png")
spec_drift("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_spectrometer_drift.png")
spec_drift("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_spectrometer_drift.png")
spec_drift("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_spectrometer_drift.png")

spec_drift("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_spectrometer_drift.png")
spec_drift("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_spectrometer_drift.png")
spec_drift("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_spectrometer_drift.png")
spec_drift("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_spectrometer_drift.png")
spec_drift("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_spectrometer_drift.png")
spec_drift("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_spectrometer_drift.png")

spec_drift("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_spectrometer_drift.png")
spec_drift("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_spectrometer_drift.png")
spec_drift("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_spectrometer_drift.png")
spec_drift("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_spectrometer_drift.png")

spec_drift("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_spectrometer_drift.png")
spec_drift("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_spectrometer_drift.png")
spec_drift("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_spectrometer_drift.png")
spec_drift("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_spectrometer_drift.png")

# Zoom in on interquartile range
spec_drift_iqr <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date within spectrometer
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)][order(spectrometer)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order within spectrometer)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

spec_drift_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_spectrometer_drift_iqr_zoom.png")

spec_drift_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_spectrometer_drift_iqr_zoom.png")

spec_drift_iqr("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_spectrometer_drift_iqr_zoom.png")

spec_drift_iqr("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_spectrometer_drift_iqr_zoom.png")
spec_drift_iqr("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_spectrometer_drift_iqr_zoom.png")

spec_drift_iqr("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_spectrometer_drift_iqr.png")
spec_drift_iqr("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_spectrometer_drift_iqr.png")
spec_drift_iqr("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_spectrometer_drift_iqr.png")
spec_drift_iqr("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_spectrometer_drift_iqr.png")

spec_drift_iqr("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_spectrometer_drift_iqr.png")
spec_drift_iqr("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_spectrometer_drift_iqr.png")
spec_drift_iqr("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_spectrometer_drift_iqr.png")
spec_drift_iqr("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_spectrometer_drift_iqr.png")

# Show concentrations vs. inter-plate variation without ordering by spectrometer
plate_by_date <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_point(aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_point(aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}


plate_by_date("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_by_date.png")
plate_by_date("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_by_date.png")
plate_by_date("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_by_date.png")
plate_by_date("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_by_date.png")
plate_by_date("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_by_date.png")
plate_by_date("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_by_date.png")

plate_by_date("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_by_date.png")
plate_by_date("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_by_date.png")
plate_by_date("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_by_date.png")
plate_by_date("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_by_date.png")
plate_by_date("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_by_date.png")
plate_by_date("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_by_date.png")

plate_by_date("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_plate_by_date.png")
plate_by_date("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_plate_by_date.png")
plate_by_date("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_plate_by_date.png")
plate_by_date("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_plate_by_date.png")
plate_by_date("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_plate_by_date.png")
plate_by_date("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_plate_by_date.png")

plate_by_date("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_plate_by_date.png")
plate_by_date("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_plate_by_date.png")
plate_by_date("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_plate_by_date.png")
plate_by_date("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_plate_by_date.png")
plate_by_date("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_plate_by_date.png")
plate_by_date("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_plate_by_date.png")

plate_by_date("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_plate_by_date.png")
plate_by_date("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_plate_by_date.png")
plate_by_date("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_plate_by_date.png")
plate_by_date("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_plate_by_date.png")

plate_by_date("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_plate_by_date.png")
plate_by_date("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_plate_by_date.png")
plate_by_date("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_plate_by_date.png")
plate_by_date("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_plate_by_date.png")

# Zoom in on interquartile range
plate_by_date_iqr <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))

  n_bio <- adj[!is.na(conc_col), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4

  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col), 
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))), 
                       by=.(plate_id, plate_MEASURED_DAYS, spectrometer, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS + spectrometer ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

  # Order plates by measurement date
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]

  g <- ggplot(plate_summary, aes(x=factor(plate_order), color=factor(spectrometer), fill=factor(spectrometer))) +
    geom_errorbar(aes(ymin=Q25, ymax=Q75), width=0, size=0.1) +
    geom_point(aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2) +
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    scale_color_calc(name="Spectrometer") +
    scale_fill_calc(name="Spectrometer") +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) + 
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))

  ggsave(g, width=width, height=width, units="in", file=fname)
}

plate_by_date_iqr("raw", "Raw biomarker concentrations", "diagnostic_plots/raw_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj2", "Adjusted for plate row", "diagnostic_plots/adj2_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj3", "Adjusted for plate column", "diagnostic_plots/adj3_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_vs_plate_by_date_iqr_zoom.png")

plate_by_date_iqr("log_raw", "Raw biomarker concentrations", "diagnostic_plots/raw_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj1", "Adjusted for sample degradation time", "diagnostic_plots/adj1_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj2", "Adjusted for plate row", "diagnostic_plots/adj2_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj3", "Adjusted for plate column", "diagnostic_plots/adj3_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj5", "(Post-QC) Excluding outlier plates", "diagnostic_plots/adj5_log_vs_plate_by_date_iqr_zoom.png")

plate_by_date_iqr("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_vs_plate_by_date_iqr_zoom.png")

plate_by_date_iqr("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj5_comb", "(Post-QC) All adjustments in a single regression, outlier plates excluded", "diagnostic_plots/adj5_comb_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj5_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression, outlier plates excluded", "diagnostic_plots/adj5_comb_strata_log_vs_plate_by_date_iqr_zoom.png")
plate_by_date_iqr("log_adj5_ts", "(Post-QC) All adjustments in a two-stage regression, outlier plates excluded", "diagnostic_plots/adj5_ts_log_vs_plate_by_date_iqr_zoom.png")

plate_by_date_iqr("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_plate_by_date_iqr.png")
plate_by_date_iqr("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_plate_by_date_iqr.png")
plate_by_date_iqr("adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_vs_plate_by_date_iqr.png")
plate_by_date_iqr("adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_vs_plate_by_date_iqr.png")

plate_by_date_iqr("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_vs_plate_by_date_iqr.png")
plate_by_date_iqr("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_vs_plate_by_date_iqr.png")
plate_by_date_iqr("log_adj5_4321", "(Post-QC) Steps reversed (4-3-2-1), outlier plates excluded", "diagnostic_plots/adj5_4321_log_vs_plate_by_date_iqr.png")
plate_by_date_iqr("log_adj5_3214", "(Post-QC) Steps reordered (3-2-1-4), outlier plates excluded", "diagnostic_plots/adj5_3214_log_vs_plate_by_date_iqr.png")

# Show concentrations vs. inter-plate variation without ordering by spectrometer
outlier_plates <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))
  
  n_bio <- adj[!is.na(conc_col) & !is.na(log_raw), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4
  
  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col) & !is.na(log_raw),
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))),
                       by=.(plate_id, plate_MEASURED_DAYS, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  
  # Order plates by measurement date 
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]
  
  # Add outlier info
  if (conc_col %like% "comb_strata") {
    plate_summary[plate_medians_comb_strata, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_comb_strata
  } else if (conc_col %like% "_ts") {
    plate_summary[plate_medians_ts, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_ts
  } else if (conc_col %like% "comb") {
    plate_summary[plate_medians_comb, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_comb
  } else if (conc_col %like% "4321") {
    plate_summary[plate_medians_4321, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_4321
  } else if (conc_col %like% "3214") {
    plate_summary[plate_medians_3214, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_3214
  } else {
    plate_summary[plate_medians, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim
  }
  
  # Split into three color groups (ggnewscale was not behaving)
  gg_no_outlier <- plate_summary[outlier == "no"]
  gg_high_outlier <- plate_summary[outlier == "high"]
  gg_low_outlier <- plate_summary[outlier == "low"]
  
  g <- ggplot(plate_summary, aes(x=factor(plate_order))) +
    geom_point(data=gg_no_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252") +
    geom_point(data=gg_no_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252") +
    geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252") +
    geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252") +

    geom_point(data=gg_high_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15") +
    geom_point(data=gg_high_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15") +
    geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15") +
    geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15") +
    
    geom_point(data=gg_low_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c") +
    geom_point(data=gg_low_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c") +
    geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c") +
    geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c") +
    
    geom_hline(data=this_outlier_lim, aes(yintercept=lower_lim), color="#dd1c77", linetype=2, size=0.2) +
    geom_hline(data=this_outlier_lim, aes(yintercept=upper_lim), color="#dd1c77", linetype=2, size=0.2) +

    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) +
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))
  
  ggsave(g, width=width, height=width, units="in", file=fname)
}

outlier_plates("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_outlier_plates.png")
outlier_plates("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_outlier_plates.png")
outlier_plates("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_outlier_plates.png")
outlier_plates("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_outlier_plates.png")
outlier_plates("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_outlier_plates.png")
outlier_plates("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_outlier_plates.png")

# Zoom in on interquartile range
outlier_plates_iqr <- function(conc_col, conc_name, fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  setnames(adj, conc_col, "conc_col")
  on.exit(setnames(adj, "conc_col", conc_col))
  
  n_bio <- adj[!is.na(conc_col) & !is.na(log_raw), length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4
  
  # Compute summary per plate
  plate_summary <- adj[!is.na(conc_col) & is.finite(conc_col) & !is.na(log_raw),
                       .(metric=names(summary(conc_col)), value=as.vector(summary(conc_col))),
                       by=.(plate_id, plate_MEASURED_DAYS, variable)]
  plate_summary <- dcast(plate_summary, variable + plate_id + plate_MEASURED_DAYS ~ metric, value.var="value")
  setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))
  
  # Order plates by measurement date 
  plate_summary <- plate_summary[order(plate_MEASURED_DAYS)]
  plate_summary[, plate_order := 1:.N]
  
  # Add outlier info
  if (conc_col %like% "comb_strata") {
    plate_summary[plate_medians_comb_strata, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_comb_strata
  } else if (conc_col %like% "_ts") {
    plate_summary[plate_medians_ts, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_ts
  } else if (conc_col %like% "comb") {
    plate_summary[plate_medians_comb, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_comb
  } else if (conc_col %like% "4321") {
    plate_summary[plate_medians_4321, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_4321
  } else if (conc_col %like% "3214") {
    plate_summary[plate_medians_3214, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim_3214
  } else {
    plate_summary[plate_medians, on = .(variable, plate_id), outlier := i.outlier]
    this_outlier_lim <- outlier_lim
  }
  
  # Split into three color groups (ggnewscale was not behaving)
  gg_no_outlier <- plate_summary[outlier == "no"]
  gg_high_outlier <- plate_summary[outlier == "high"]
  gg_low_outlier <- plate_summary[outlier == "low"]
  
  g <- ggplot(plate_summary, aes(x=factor(plate_order))) +
    geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252") +
    geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252") +
    
    geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15") +
    geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15") +
    
    geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c") +
    geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c") +
    
    geom_hline(data=this_outlier_lim, aes(yintercept=lower_lim), color="#dd1c77", linetype=2, size=0.2) +
    geom_hline(data=this_outlier_lim, aes(yintercept=upper_lim), color="#dd1c77", linetype=2, size=0.2) +
    
    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(paste0(conc_name, "\n", "(Median, interquartile range, minimum, and maximum values)")) +
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(color = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1)))
  
  ggsave(g, width=width, height=width, units="in", file=fname)
}

outlier_plates_iqr("adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_outlier_plates_iqr_zoom.png")
outlier_plates_iqr("adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_outlier_plates_iqr_zoom.png")
outlier_plates_iqr("adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_outlier_plates_iqr_zoom.png")
outlier_plates_iqr("adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_outlier_plates_iqr_zoom.png")
outlier_plates_iqr("adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_vs_outlier_plates_iqr.png")
outlier_plates_iqr("adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_vs_outlier_plates_iqr.png")

# Now show also on log scale
outlier_lim[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj4_2 := log(adj4 - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj4_2 ~ 1))[1]) - as.vector(coef(glm(log_adj4 ~ 1))[1])), by=variable]
outlier_lim[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_outlier_plates.png")
outlier_plates_iqr("log_adj4", "(Post-QC) Adjusted for spectrometer drift over time", "diagnostic_plots/adj4_log_outlier_plates_iqr_zoom.png")

outlier_lim_comb[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim_comb[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj_comb_2 := log(adj_comb - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj_comb_2 ~ 1))[1]) - as.vector(coef(glm(log_adj_comb ~ 1))[1])), by=variable]
outlier_lim_comb[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim_comb[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_outlier_plates.png")
outlier_plates_iqr("log_adj_comb", "(Post-QC) All adjustments in a single regression", "diagnostic_plots/adj_comb_log_outlier_plates_iqr_zoom.png")

outlier_lim_comb_strata[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim_comb_strata[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj_comb_strata_2 := log(adj_comb_strata - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj_comb_strata_2 ~ 1))[1]) - as.vector(coef(glm(log_adj_comb_strata ~ 1))[1])), by=variable]
outlier_lim_comb_strata[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim_comb_strata[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_outlier_plates.png")
outlier_plates_iqr("log_adj_comb_strata", "(Post-QC) All adjustments in a single spectrometer stratified regression", "diagnostic_plots/adj_comb_strata_log_outlier_plates_iqr_zoom.png")

outlier_lim_ts[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim_ts[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj_ts_2 := log(adj_ts - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj_ts_2 ~ 1))[1]) - as.vector(coef(glm(log_adj_ts ~ 1))[1])), by=variable]
outlier_lim_ts[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim_ts[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_outlier_plates.png")
outlier_plates_iqr("log_adj_ts", "(Post-QC) All adjustments in a two-stage regression", "diagnostic_plots/adj_ts_log_outlier_plates_iqr_zoom.png")

outlier_lim_4321[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim_4321[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj4_4321_2 := log(adj4_4321 - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj4_4321_2 ~ 1))[1]) - as.vector(coef(glm(log_adj4_4321 ~ 1))[1])), by=variable]
outlier_lim_4321[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim_4321[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_outlier_plates.png")
outlier_plates_iqr("log_adj4_4321", "(Post-QC) Steps reversed (4-3-2-1)", "diagnostic_plots/adj4_4321_log_outlier_plates_iqr_zoom.png")

outlier_lim_3214[log_offset, on = .(variable), lower_lim := log(lower_lim - right_shift + log_offset)]
outlier_lim_3214[log_offset, on = .(variable), upper_lim := log(upper_lim - right_shift + log_offset)]
adj[log_offset, on = .(variable), log_adj4_3214_2 := log(adj4_3214 - right_shift + log_offset)]
outlier_shift <- adj[!is.na(log_raw),.(shift=as.vector(coef(glm(log_adj4_3214_2 ~ 1))[1]) - as.vector(coef(glm(log_adj4_3214 ~ 1))[1])), by=variable]
outlier_lim_3214[outlier_shift, on = .(variable), lower_lim := lower_lim - shift]
outlier_lim_3214[outlier_shift, on = .(variable), upper_lim := upper_lim - shift]

outlier_plates("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_outlier_plates.png")
outlier_plates_iqr("log_adj4_3214", "(Post-QC) Steps reordered (3-2-1-4)", "diagnostic_plots/adj4_3214_log_outlier_plates_iqr_zoom.png")

# Show internal controls
cntrl <- fread("data/raw/processed/NGH_control_samples.txt")
cntrl_info <- fread("data/raw/processed/controls_information.txt")

cntrl <- melt(cntrl, id.vars=c("sample_id", "plate_id", "plate_position", "spectrometer"))
cntrl <- cntrl[variable %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]
cntrl[cntrl_info, on = .(sample_id, plate_id, plate_position, spectrometer), plate_MEASURED_DAYS := plate_MEASURED_DAYS]
cntrl <- cntrl[order(tolower(variable))]
cntrl[, variable := factor(variable, levels=unique(variable))]
cntrl[plate_medians, on = .(variable, plate_id), outlier := i.outlier] 
cntrl[, control_sample := sprintf("%s (Well %s)", sample_id, plate_position)]
cntrl[, control_sample := factor(control_sample, levels=c( # Control samples are paired
        "190404 (Well A01)", "190328 (Well H12)",
        "190425 (Well A01)", "190508 (Well H12)",
        "180829 (Well A01)", "180827 (Well H12)",
        "180830 (Well A01)", "180831 (Well H12)"
      ))]

internal_controls <- function(fname) {
  if (file.exists(fname)) {
    message(sprintf("%s already exists, skipping", fname))
    invisible(return(NULL))
  }

  n_bio <- cntrl[, length(unique(variable))]
  n_col <- ceiling(sqrt(n_bio))
  width <- n_col * 1.4
  
  # Order plates by measurement date 
  cntrl <- cntrl[order(plate_MEASURED_DAYS)]
  cntrl[, plate_order := 1:.N]
  
  # Split into three color groups (ggnewscale was not behaving)
  gg_no_outlier <- cntrl[outlier == "no"]
  gg_high_outlier <- cntrl[outlier == "high"]
  gg_low_outlier <- cntrl[outlier == "low"]
  
  g <- ggplot(cntrl, aes(x=factor(plate_order), shape=factor(control_sample))) +
    scale_shape_manual(name = "Control Sample (4-pairs used across all plates)",
      values=c( # Control samples are paired
				"190404 (Well A01)"=0, "190328 (Well H12)"=22,
				"190425 (Well A01)"=1, "190508 (Well H12)"=21,
				"180829 (Well A01)"=2, "180827 (Well H12)"=24,
				"180830 (Well A01)"=5, "180831 (Well H12)"=23
			)) +

    geom_point(data=gg_no_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#525252") +
    geom_point(data=gg_no_outlier[plate_position == "A01"], aes(y=value), color="#525252", size=0.45, stroke=0.2) +

    geom_point(data=gg_high_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#a50f15") +
    geom_point(data=gg_high_outlier[plate_position == "A01"], aes(y=value), color="#a50f15", size=0.45, stroke=0.2) +
    
    geom_point(data=gg_low_outlier[plate_position == "H12"], aes(y=value), color="black", size=0.45, stroke=0.2, fill="#08519c") +
    geom_point(data=gg_low_outlier[plate_position == "A01"], aes(y=value), color="#08519c", size=0.45, stroke=0.2) +

    facet_wrap(~ variable, scales="free", ncol=n_col) +
    ylab(paste0("Concentration of control sample")) +
    xlab("Plate\n(chronological order)") +
    theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          strip.background=element_blank(), strip.text=element_text(size=6),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          legend.position="bottom", legend.text=element_text(size=10)
    ) +
    guides(shape = guide_legend(override.aes = list(size=1, stroke=0.8, alpha=1, colour="black", fill="#525252")))
  
  ggsave(g, width=width, height=width, units="in", file=fname)
}

internal_controls("diagnostic_plots/internal_control_samples_outlier_plates.png")

# Show on log scale too
cntrl[log_offset, on = .(variable), value := log(value - right_shift + log_offset)]
cntrl[outlier_shift, on = .(variable), value := value - shift]

internal_controls("diagnostic_plots/internal_control_samples_log_outlier_plates.png")

# For biomarkers also present in the clinical biochemistry, do the outlier plates
# appear as outliers?
clin <- fread("data/my_curated_phenotypes/biomarkers/output/biomarkers.txt")
clin <- melt(clin, id.vars=c("eid", "visit_index"))
setnames(clin, "visit_index", "visit")

# Create map between blood biochemistry and Nightingale biomarkers
map <- rbind(use.names=FALSE,
  data.table(clin_var="tchol", clin_name="Total cholesterol", clin_units="mmol/L",
             nightingale_var="Total_C", nightingale_name="Total cholesterol", nightingale_units="mmol/L"),
  data.table("ldl", "LDL cholesterol", "mmol/L", "Clinical_LDL_C", "Clinical LDL cholesterol", "mmol/L"),
  data.table("ldl", "LDL cholesterol", "mmol/L", "LDL_C", "LDL cholesterol", "mmol/L"),
  data.table("hdl", "HDL cholesterol", "mmol/L", "HDL_C", "HDL cholesterol", "mmol/L"),
  data.table("trig", "Total triglycerides", "mmol/L", "Total_TG", "Total triglycerides", "mmol/L"),
  data.table("apoa1", "Apolipoprotein A1", "g/L", "ApoA1", "Apolipoprotein A1", "g/L"),
  data.table("apob", "Apolipoprotein B", "g/L", "ApoB", "Apolipoprotein B", "g/L"),
  data.table("glucose1", "Glucose", "mmol/L", "Glucose", "Glucose", "mmol/L"),
  data.table("creat", "Creatinine", "umol/L", "Creatinine", "Creatinine", "mmol/L"),
  data.table("alb", "Albumin", "g/L", "Albumin", "Albumin", "g/L"),
  data.table("nonhdl", "Non HDL Cholesterol", "mmol/L", "non_HDL_C", "Non HDL Cholesterol", "mmol/L"),
  data.table("apobapoa1", "ApoB / ApoA1", "ratio", "ApoB_by_ApoA1", "ApoB / ApoA1", "ratio")
)

clin <- clin[variable %in% map$clin_var]
clin[, visit := ifelse(visit == 0, "Main Phase", "Repeat Assessment")]
clin[sinfo[!(sample_removed) & (in_ukb_raw)], on = .(eid=eid_7439, visit), c("plate_id", "plate_MEASURED_DATE") := .(plate_id, plate_MEASURED_DATE)]
clin <- clin[!is.na(plate_id)]
clin[map, on = .(variable=clin_var), c("nightingale_var", "clin_name") := .(nightingale_var, clin_name)]

plate_outliers <- copy(plate_medians)
plate_outliers[outlier == "no", outlier := NA_character_]
plate_outliers <- dcast(plate_outliers, plate_id ~ variable, value.var="outlier")
plate_outliers <- ukbnmr::recompute_derived_biomarker_qc_flags(plate_outliers)
plate_outliers <- melt(plate_outliers, id.vars="plate_id")
plate_outliers <- plate_outliers[variable %in% map$nightingale_var & !is.na(value)]
plate_outliers[,tag := paste(unique(gsub("^.*: ", "", strsplit(value, "\\.")[[1]])), collapse=";"), by=.(plate_id, variable)]
plate_outliers[tag == "high;low", tag := "high"] # Total_TG, 9 lipo-protein sub fractions are high outliers, 1 is low (XXL_VLDL_TG), 4 not outliers.

clin[, outlier := "no"]
clin[plate_outliers, on = .(nightingale_var=variable, plate_id), outlier := tag]

# Compute summary per plate
plate_summary <- clin[, .(metric=names(summary(value)), value=as.vector(summary(value))), by=.(clin_name, plate_id, plate_MEASURED_DATE, outlier)]
plate_summary <- dcast(plate_summary, clin_name + plate_id + plate_MEASURED_DATE + outlier ~ metric, value.var="value")
setnames(plate_summary, c("1st Qu.", "3rd Qu.", "Max.", "Min."), c("Q25", "Q75", "Max", "Min"))

# Order plates by measurement date
plate_summary <- plate_summary[order(plate_MEASURED_DATE)]
plate_summary[, plate_order := 1:.N]

# Split into three color groups
gg_no_outlier <- plate_summary[outlier == "no"]
gg_high_outlier <- plate_summary[outlier == "high"]
gg_low_outlier <- plate_summary[outlier == "low"]

# How many biomarkers?
n_bio <- clin[, length(unique(variable))]
n_col <- ceiling(sqrt(n_bio))
width <- n_col * 2

# Show with min and max
g <- ggplot(plate_summary, aes(x=factor(plate_order))) +
  geom_point(data=gg_no_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252") +
  geom_point(data=gg_no_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#525252") +
  geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252") +
  geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252") +

  geom_point(data=gg_high_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15") +
  geom_point(data=gg_high_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#a50f15") +
  geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15") +
  geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15") +

  geom_point(data=gg_low_outlier, aes(y=Min), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c") +
  geom_point(data=gg_low_outlier, aes(y=Max), shape=1, size=0.3, stroke=0.15, alpha=0.5, color="#08519c") +
  geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c") +
  geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c") +

  facet_wrap(~ clin_name, scales="free", ncol=n_col) +
	ylab("Clinical chemistry concentration") +
	xlab("Plate\n(chronological order)") +
	theme_bw() +
	theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				strip.background=element_blank(), strip.text=element_text(size=6),
				panel.grid.major=element_blank(), panel.grid.minor=element_blank()
	)
ggsave(g, width=width, height=width, units="in", file="diagnostic_plots/clinical_chemistry_vs_nmr_outlier_plates.png")


# Zoom in on IQR
g <- ggplot(plate_summary, aes(x=factor(plate_order))) +
  geom_errorbar(data=gg_no_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#525252") +
  geom_point(data=gg_no_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#525252") +

  geom_errorbar(data=gg_high_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#a50f15") +
  geom_point(data=gg_high_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#a50f15") +

  geom_errorbar(data=gg_low_outlier, aes(ymin=Q25, ymax=Q75), width=0, size=0.1, color="#08519c") +
  geom_point(data=gg_low_outlier, aes(y=Median), color="black", shape=21, size=0.45, stroke=0.2, fill="#08519c") +

  facet_wrap(~ clin_name, scales="free", ncol=n_col) +
	ylab("Clinical chemistry concentration") +
	xlab("Plate\n(chronological order)") +
	theme_bw() +
	theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				strip.background=element_blank(), strip.text=element_text(size=6),
				panel.grid.major=element_blank(), panel.grid.minor=element_blank()
	)
ggsave(g, width=width, height=width, units="in", file="diagnostic_plots/clinical_chemistry_vs_nmr_outlier_plates_iqr_zoom.png")

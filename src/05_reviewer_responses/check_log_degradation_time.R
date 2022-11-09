library(data.table)
library(ggplot2)
library(cowplot)

# make output directory
system("mkdir -p reviewer_responses/")

##' # Load multi-step adjusted concentrations
##' multi_adj <- fread("data/tech_qc/multistep_adjusted_values.txt")
##' 
##' # Add in concentrations from preprint QC procedure
##' preprint <- fread("data/preprint/nmr_techadj.txt")
##' preprint <- melt(preprint, id.vars=c("eid", "visit_index"))
##' preprint[, visit := ifelse(visit_index == 0, "Main Phase", "Repeat Assessment")]
##' multi_adj[preprint, on = .(eid_7439=eid, visit, variable), preprint := i.value]
##' 
##' # Filter to non-derived biomarkers
##' nmr_info <- fread("data/preprint/biomarker_information.txt")
##' multi_adj <- multi_adj[variable %in% nmr_info[Type == "Non-derived", Biomarker]]
##' 
##' # Compute correlations between old and new QC values
##' qc_cor <- multi_adj[,.(
##'   pearson=cor(adj4, preprint, use="pairwise.complete.obs"),
##'   spearman=cor(adj4, preprint, use="pairwise.complete.obs", method="spearman")
##' ), by=variable]
##' 
##' fwrite(qc_cor, sep="\t", quote=FALSE, file="reviewer_responses/cor_qc_preprint_vs_new.txt")

# Extract histidine for plotting
his <- multi_adj[variable == "His"]

# First, show raw histidine vs. sample degradation time
g1 <- ggplot(his, aes(x=sample_degredation, y=raw)) +
  geom_hex() +
  scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
  geom_smooth(color="red", size=0.2, method=MASS::rlm) +
  ylab("Histidine (mmol/L)") +
  scale_x_continuous("Sample degradation time\n(hours from sample prep to measurement)") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.text=element_text(size=6), legend.title=element_text(size=7)
  )

g2 <- ggplot(his, aes(x=sample_degredation, y=raw)) +
  geom_hex() +
  scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
  geom_smooth(color="red", size=0.2, method=MASS::rlm) +
  ylab("Histidine (mmol/L)") +
  scale_x_log10(bquote("On a" ~ log[10] ~ "scale")) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.text=element_text(size=6), legend.title=element_text(size=7)
  )

# Then after adjusting for log-degradation time
g3 <- ggplot(his, aes(x=sample_degredation, y=adj4)) +
  geom_hex() +
  scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
  geom_smooth(color="red", size=0.2, method=MASS::rlm) +
  ylab("Histidine (mmol/L)") +
  scale_x_continuous("Sample degradation time\n(hours from sample prep to measurement)") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.text=element_text(size=6), legend.title=element_text(size=7)
  )

g4 <- ggplot(his, aes(x=sample_degredation, y=adj4)) +
  geom_hex() +
  scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
  geom_smooth(color="red", size=0.2, method=MASS::rlm) +
  ylab("Histidine (mmol/L)") +
  scale_x_log10(bquote("On a" ~ log[10] ~ "scale")) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.text=element_text(size=6), legend.title=element_text(size=7)
  )


# Compared to preprint
g5 <- ggplot(his, aes(x=sample_degredation, y=preprint)) +
  geom_hex() +
  scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
  geom_smooth(color="red", size=0.2, method=MASS::rlm) +
  ylab("Histidine (mmol/L)") +
  scale_x_continuous("Sample degradation time\n(hours from sample prep to measurement)") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.text=element_text(size=6), legend.title=element_text(size=7)
  )

g6 <- ggplot(his, aes(x=sample_degredation, y=preprint)) +
  geom_hex() +
  scale_fill_gradient(name="Samples", low="#d9d9d9", high="#252525", trans="log10", limits=c(1,10000)) +
  geom_smooth(color="red", size=0.2, method=MASS::rlm) +
  ylab("Histidine (mmol/L)") +
  scale_x_log10(bquote("On a" ~ log[10] ~ "scale")) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.text=element_text(size=6), legend.title=element_text(size=7)
  )


g <- plot_grid(g1, g2, g3, g4, g5, g6, ncol=2, align="hv", label_size=8, hjust=0,
  labels=c("Original histidine", "", "Adjusted for log degradation time", "", "Adjusted for linear degradation time", ""))
ggsave(g, width=6, height=8, file="reviewer_responses/log_degradation_time_example.png")









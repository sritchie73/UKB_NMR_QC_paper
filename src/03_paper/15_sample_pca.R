library(data.table)
library(ukbnmr)
library(foreach)
library(ggplot2)
library(ggrastr)
library(ggnewscale)
library(cowplot)
library(palr)
library(RColorBrewer)
source("src/utilities/logit.R")

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

# For samples with data at both baseline and repeat assessment, take baseline data
baseline <- sinfo[!(sample_removed) & (in_ukb_raw) & visit == "Main Phase"]
repeat_visit <- sinfo[!(sample_removed) & (in_ukb_raw) & visit == "Repeat Assessment"]
repeat_visit <- repeat_visit[!baseline, on = .(eid_7439)]
first <- rbind(baseline, repeat_visit)
dat <- dat[first[,.(eid_7439, sample_id, visit)], on = .(eid_7439, sample_id, visit)]

# Extract target raw and post-qc data
raw <- dat[,.(eid_7439, variable, value=raw)][!is.na(value)]
adj <- dat[,.(eid_7439, variable, value=adj5)][!is.na(value)]

# Log / logit transform, adding small offset for biomarkers with concentrations of 0 (half the minimum non-zero value).
# For percentages, the maximum value is also accounted for when determining the offset to ensure that values don't become
# greater than 1 after applying the offset.
raw[!is.na(value) & !is.finite(value), value := NA]
raw <- raw[!is.na(value)]
raw[variable %in% nmr_info[Units == "%", Biomarker], value := value/100]
lim_value <- raw[, .(min=min(value), non_zero=min(value[value != 0]), non_one=max(value[value != 1]), max=max(value)), by=variable]
lim_value[variable %in% nmr_info[Units != "%", Biomarker], offset := ifelse(min == 0, non_zero / 2, 0)]
lim_value[variable %in% nmr_info[Units == "%", Biomarker], offset := fcase(
  min == 0 & max == 1, 0, # cannot offset
  min == 0, pmin(non_zero / 2, (1 - max)/2), # make sure offset doesnt make % > 100
  max == 1, -pmin(min / 2, (1 - non_one)/2), # make sure offset doesnt make % < 0
  default = 0)]
raw[lim_value, on = .(variable), value := value + offset]
raw[variable %in% nmr_info[Units == "%", Biomarker], value := logit(value)]
raw[variable %in% nmr_info[Units != "%", Biomarker], value := log(value)]

adj[!is.na(value) & !is.finite(value), value := NA]
adj <- adj[!is.na(value)]
adj[variable %in% nmr_info[Units == "%", Biomarker], value := value/100]
lim_value <- adj[, .(min=min(value), non_zero=min(value[value != 0]), non_one=max(value[value != 1]), max=max(value)), by=variable]
lim_value[variable %in% nmr_info[Units != "%", Biomarker], offset := ifelse(min == 0, non_zero / 2, 0)]
lim_value[variable %in% nmr_info[Units == "%", Biomarker], offset := fcase(
  min == 0 & max == 1, 0, # cannot offset
  min == 0, pmin(non_zero / 2, (1 - max)/2), # make sure offset doesnt make % > 100
  max == 1, -pmin(min / 2, (1 - non_one)/2), # make sure offset doesnt make % < 0
  default = 0)]
adj[lim_value, on = .(variable), value := value + offset]
adj[variable %in% nmr_info[Units == "%", Biomarker], value := logit(value)]
adj[variable %in% nmr_info[Units != "%", Biomarker], value := log(value)]

# Standardise to mean 0 sd of 1
raw[, value := scale(value), by=variable]
adj[, value := scale(value), by=variable]

# Cast to wide with participants as columns
raw <- dcast(raw, variable ~ eid_7439, value.var="value", fill=0) # missing values imputed as mean for PCA
raw <- as.matrix(raw, rownames="variable")

adj <- dcast(adj, variable ~ eid_7439, value.var="value", fill=0) # missing values imputed as mean for PCA
adj <- as.matrix(adj, rownames="variable")

# Run PCA
raw_pcs <- prcomp(raw, scale.=TRUE)
adj_pcs <- prcomp(adj, scale.=TRUE)

# Extract summary of PCs
raw_pc_summary <- as.data.table(summary(raw_pcs)$importance, keep.rownames="metric")
raw_pc_summary <- melt(raw_pc_summary, id.vars="metric", variable.name="PC")

adj_pc_summary <- as.data.table(summary(adj_pcs)$importance, keep.rownames="metric")
adj_pc_summary <- melt(adj_pc_summary, id.vars="metric", variable.name="PC")

##
# Main fig plots
## 

# Show variance explained by PCs explaining > 1% of the variance
pc_gt1 <- adj_pc_summary[metric == "Proportion of Variance"][value > 0.01]
pc_gt1[, PC := as.integer(gsub("PC", "", PC))]
g1 <- ggplot(pc_gt1) + 
  aes(x=PC, y=value * 100) +
  geom_line(size=0.3, colour="#209058") +
  geom_point(shape=21, fill="#209058", colour="black", stroke=0.2, size=0.5) +
  scale_x_continuous(name="Principal component", breaks=seq(1, 12, by=2)) +
  ylab("Variance explained (%)") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.25), panel.grid.minor=element_line(size=0.1)
  )

# Extract PCs and show separation by sex
pc_gt1 <- as.data.table(adj_pcs$rotation[, pc_gt1$PC], keep.rownames="eid_7439")
pc_gt1[, eid_7439 := as.integer(eid_7439)]
pc_gt1[sinfo, on = .(eid_7439), Sex := i.Gender]
g2 <- ggplot(pc_gt1) + aes(x = PC1, color=Sex) +
  geom_density(size=0.3) +
  scale_color_manual(values=c("Male"="#D86020", "Female"="#184068")) +
  scale_x_continuous(name="PC1", limits=c(-0.0051, 0.0051), breaks=c(-0.005, 0, 0.005)) +
  ylab("Density") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none"
  )

g3 <- ggplot(pc_gt1) + aes(x = PC3, color=Sex) +
  geom_density(size=0.3) +
  scale_color_manual(values=c("Male"="#D86020", "Female"="#184068")) +
  scale_x_continuous(name="PC3", limits=c(-0.0075, 0.0075), breaks=c(-0.0075, 0, 0.0075)) +
  ylab("Density") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none"
  )

# Plot PC1 vs. PC3, the ones which most stratify by sex
g4 <- ggplot(pc_gt1) +
  aes(x = PC1, y = PC3) +
  geom_hex(data=pc_gt1[Sex == "Male"]) +
  scale_fill_gradient(name="Males", low="#D860201a", high="#D86020ff") + 
  new_scale_fill() +
  geom_hex(data=pc_gt1[Sex == "Female"]) +
  scale_fill_gradient(name="Females", low="#1840681a", high="#184068ff") + 
  scale_x_continuous(name="PC1", limits=c(-0.006, 0.006), breaks=c(-0.0050, 0, 0.0050)) +
  scale_y_continuous(name="PC3", limits=c(-0.008, 0.008), breaks=c(-0.0075, 0, 0.0075)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.title=element_text(size=7), legend.text=element_text(size=6), 
        legend.box="horizontal", legend.box.margin=margin(0, 0, 0, 0, unit = "pt"),
        legend.key.height=unit(0.1, "in"), legend.key.width=unit(0.1, "in")
  )

# Arrange plots
g23 <- plot_grid(g2, g3, nrow=2)
g <- plot_grid(g1, g23, g4, rel_widths=c(1,1,2), nrow=1)
ggsave(g, width=5, height=1, file="paper_output/postqc_pcs_by_sex.pdf")

##
# Replication in raw data
## 

# Show variance explained by PCs explaining > 1% of the variance
pc_gt1 <- raw_pc_summary[metric == "Proportion of Variance"][value > 0.01]
pc_gt1[, PC := as.integer(gsub("PC", "", PC))]
g1 <- ggplot(pc_gt1) + 
  aes(x=PC, y=value * 100) +
  geom_line(size=0.3, colour="#209058") +
  geom_point(shape=21, fill="#209058", colour="black", stroke=0.2, size=0.5) +
  scale_x_continuous(name="Principal component", breaks=seq(1, 12, by=2)) +
  ylab("Variance explained (%)") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_line(size=0.25), panel.grid.minor=element_line(size=0.1)
  )

# Extract PCs and show separation by sex
pc_gt1 <- as.data.table(raw_pcs$rotation[, pc_gt1$PC], keep.rownames="eid_7439")
pc_gt1[, eid_7439 := as.integer(eid_7439)]
pc_gt1[sinfo, on = .(eid_7439), Sex := i.Gender]
g2 <- ggplot(pc_gt1) + aes(x = PC1, color=Sex) +
  geom_density(size=0.3) +
  scale_color_manual(values=c("Male"="#D86020", "Female"="#184068")) +
  scale_x_continuous(name="PC1", limits=c(-0.0051, 0.0051), breaks=c(-0.005, 0, 0.005)) +
  ylab("Density") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none"
  )

g3 <- ggplot(pc_gt1) + aes(x = PC3, color=Sex) +
  geom_density(size=0.3) +
  scale_color_manual(values=c("Male"="#D86020", "Female"="#184068")) +
  scale_x_continuous(name="PC3", limits=c(-0.0075, 0.0075), breaks=c(-0.0075, 0, 0.0075)) +
  ylab("Density") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none"
  )

# Plot PC1 vs. PC3, the ones which most stratify by sex
g4 <- ggplot(pc_gt1) +
  aes(x = PC1, y = PC3) +
  geom_hex(data=pc_gt1[Sex == "Male"]) +
  scale_fill_gradient(name="Males", low="#D860201a", high="#D86020ff") + 
  new_scale_fill() +
  geom_hex(data=pc_gt1[Sex == "Female"]) +
  scale_fill_gradient(name="Females", low="#1840681a", high="#184068ff") + 
  scale_x_continuous(name="PC1", limits=c(-0.006, 0.006), breaks=c(-0.0050, 0, 0.0050)) +
  scale_y_continuous(name="PC3", limits=c(-0.008, 0.008), breaks=c(-0.0075, 0, 0.0075)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.title=element_text(size=7), legend.text=element_text(size=6), 
        legend.box="horizontal", legend.box.margin=margin(0, 0, 0, 0, unit = "pt"),
        legend.key.height=unit(0.1, "in"), legend.key.width=unit(0.1, "in")
  )

# Arrange plots
g23 <- plot_grid(g2, g3, nrow=2)
g <- plot_grid(g1, g23, g4, rel_widths=c(1,1,2), nrow=1)
ggsave(g, width=5, height=1, file="paper_output/raw_pcs_by_sex.pdf")

# Build pairwise plots akin to GGally::ggpairs. Top half shows raw data,
# bottom shows post-qc.
npcs <- max(adj_pc_summary[metric == "Proportion of Variance", sum(value > 0.01)],
            raw_pc_summary[metric == "Proportion of Variance", sum(value > 0.01)])
pc_gt1 <- rbind(idcol="dataset",
  raw=as.data.table(raw_pcs$rotation[, 1:npcs], keep.rownames="eid_7439"),
  adj=as.data.table(adj_pcs$rotation[, 1:npcs], keep.rownames="eid_7439")
)
gl <- foreach(pcii = 1:npcs, .combine=c) %:% foreach(pcjj = 1:npcs) %do% {
  if (pcjj == pcii) {
		ggplot(pc_gt1) + 
      aes_string(x=sprintf("`PC%s`", pcii)) +
			geom_density(size=0.3, aes(color=dataset), trim=TRUE) +
      scale_color_manual(values=c("raw"="#e08214", "adj"="#542788")) +
      ylim(0, 170) +
			theme_bw() +
			theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
						axis.text.y=element_blank(), axis.title.y=element_blank(),
            axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt"),
						panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"), legend.position="none"
			)
  } else if (pcjj < pcii) {
		ggplot(pc_gt1[dataset == "adj"]) + 
      aes_string(x=sprintf("`PC%s`", pcjj), y=sprintf("`PC%s`", pcii)) +
      rasterize(geom_hex()) +
      scale_fill_gradientn(colours=sst_pal(256), limits=c(1, 3500)) +
			theme_bw() +
			theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
						axis.text.y=element_blank(), axis.title.y=element_blank(),
            axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt"),
						panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"), legend.position="none"
      )
  } else {
		ggplot(pc_gt1[dataset == "raw"]) + 
      aes_string(x=sprintf("`PC%s`", pcjj), y=sprintf("`PC%s`", pcii)) +
      rasterize(geom_hex()) +
      scale_fill_gradientn(colours=sst_pal(256), limits=c(1, 3500)) +
			theme_bw() +
			theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
						axis.text.y=element_blank(), axis.title.y=element_blank(),
            axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt"),
						panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"), legend.position="none"
      )

  }
}
gln <- foreach(pcii = 1:npcs, .combine=c) %:% foreach(pcjj = 1:npcs, .combine=c) %do% {
  sprintf("g_%s_%s", pcii, pcjj)
}
names(gl) <- gln

g <- plot_grid(
  gl$g_1_1, gl$g_1_2, gl$g_1_3, gl$g_1_4, gl$g_1_5, gl$g_1_6, gl$g_1_7, gl$g_1_8, gl$g_1_9, gl$g_1_10, gl$g_1_11, gl$g_1_12,
  gl$g_2_1, gl$g_2_2, gl$g_2_3, gl$g_2_4, gl$g_2_5, gl$g_2_6, gl$g_2_7, gl$g_2_8, gl$g_2_9, gl$g_2_10, gl$g_2_11, gl$g_2_12,
  gl$g_3_1, gl$g_3_2, gl$g_3_3, gl$g_3_4, gl$g_3_5, gl$g_3_6, gl$g_3_7, gl$g_3_8, gl$g_3_9, gl$g_3_10, gl$g_3_11, gl$g_3_12,
  gl$g_4_1, gl$g_4_2, gl$g_4_3, gl$g_4_4, gl$g_4_5, gl$g_4_6, gl$g_4_7, gl$g_4_8, gl$g_4_9, gl$g_4_10, gl$g_4_11, gl$g_4_12,
  gl$g_5_1, gl$g_5_2, gl$g_5_3, gl$g_5_4, gl$g_5_5, gl$g_5_6, gl$g_5_7, gl$g_5_8, gl$g_5_9, gl$g_5_10, gl$g_5_11, gl$g_5_12,
  gl$g_6_1, gl$g_6_2, gl$g_6_3, gl$g_6_4, gl$g_6_5, gl$g_6_6, gl$g_6_7, gl$g_6_8, gl$g_6_9, gl$g_6_10, gl$g_6_11, gl$g_6_12,
  gl$g_7_1, gl$g_7_2, gl$g_7_3, gl$g_7_4, gl$g_7_5, gl$g_7_6, gl$g_7_7, gl$g_7_8, gl$g_7_9, gl$g_7_10, gl$g_7_11, gl$g_7_12,
  gl$g_8_1, gl$g_8_2, gl$g_8_3, gl$g_8_4, gl$g_8_5, gl$g_8_6, gl$g_8_7, gl$g_8_8, gl$g_8_9, gl$g_8_10, gl$g_8_11, gl$g_8_12,
  gl$g_9_1, gl$g_9_2, gl$g_9_3, gl$g_9_4, gl$g_9_5, gl$g_9_6, gl$g_9_7, gl$g_9_8, gl$g_9_9, gl$g_9_10, gl$g_9_11, gl$g_9_12,
  gl$g_10_1, gl$g_10_2, gl$g_10_3, gl$g_10_4, gl$g_10_5, gl$g_10_6, gl$g_10_7, gl$g_10_8, gl$g_10_9, gl$g_10_10, gl$g_10_11, gl$g_10_12,
  gl$g_11_1, gl$g_11_2, gl$g_11_3, gl$g_11_4, gl$g_11_5, gl$g_11_6, gl$g_11_7, gl$g_11_8, gl$g_11_9, gl$g_11_10, gl$g_11_11, gl$g_11_12,
  gl$g_12_1, gl$g_12_2, gl$g_12_3, gl$g_12_4, gl$g_12_5, gl$g_12_6, gl$g_12_7, gl$g_12_8, gl$g_12_9, gl$g_12_10, gl$g_12_11, gl$g_12_12,
  align="hv", nrow=12, ncol=12
)
ggsave(g, width=178, height=178, units="mm", file="paper_output/pc_pairs.pdf")

# Also plot larger version with legend so we can add legend to plot
g1 <- ggplot(pc_gt1) + aes(x=PC11, y=PC12) +
  rasterize(geom_hex()) +
  scale_fill_gradientn(colours=sst_pal(256), limits=c(1, 3500)) +
	theme_bw() +
	theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
				axis.text.y=element_blank(), axis.title.y=element_blank(),
				axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt"),
				panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
				plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"),

 	)
ggsave(g1, width=4, height=4, file="paper_output/pc_pairs_count_legend.pdf")

# Write out table of variance explained
raw_pc_summary <- dcast(raw_pc_summary, PC ~ metric, value.var="value")
setnames(raw_pc_summary, c("PC", "cumul_var_expl", "var_expl", "std_dev"))
raw_pc_summary[, var_expl := var_expl*100]
raw_pc_summary[, cumul_var_expl := cumul_var_expl*100]
raw_pc_summary <- raw_pc_summary[, .(PC, std_dev, var_expl, cumul_var_expl)]
fwrite(raw_pc_summary, sep="\t", quote=FALSE, file="paper_output/raw_pc_importance.txt")

adj_pc_summary <- dcast(adj_pc_summary, PC ~ metric, value.var="value")
setnames(adj_pc_summary, c("PC", "cumul_var_expl", "var_expl", "std_dev"))
adj_pc_summary[, var_expl := var_expl*100]
adj_pc_summary[, cumul_var_expl := cumul_var_expl*100]
adj_pc_summary <- adj_pc_summary[, .(PC, std_dev, var_expl, cumul_var_expl)]
fwrite(adj_pc_summary, sep="\t", quote=FALSE, file="paper_output/adj_pc_importance.txt")

# Load in information about participant age, sex, BMI, white british/other,
# lipid lowering medication, and fasting status to test association with PCs 
# (this seems like a reasonable starting point for where we might expect to 
# see major separation of samples in the lipid and metabolite data). 
#
# Data has been curated by the Cardiovascular Epidemiology Unit using Stata,
# variable labels have been derived and assigned names in a consistent manner
# across departmental resources.
library(readstata13)

ext_rep <- read.dta13("data/ceu_curated_phenotypes/20210302/STATA/repeats.dta")
ext_rep <- as.data.table(ext_rep)
pheno <- ext_rep[,.(eid_7439=idno, visit=repno, age=ages_rep, bmi=bmi_rep, lipid_med=lipdbin_rep, fasting_time=dummyfast)]
pheno <- pheno[visit %in% c(0,1)]
pheno[, visit := ifelse(visit == 0, "Main Phase", "Repeat Assessment")]
pheno[, lipid_med := ifelse(lipid_med == "Current", TRUE, FALSE)]

# need to load in information sheet with baseline data only to get white british designation and sex
ext_base_only <- read.dta13("data/ceu_curated_phenotypes/20210302/STATA/analysis.dta")
ext_base_only <- as.data.table(ext_base_only)
pheno[ext_base_only, on = .(eid_7439=idno), c("sex", "white_british") := .(sex, racebin)]
pheno[, sex := factor(sex, levels=c("Female", "Male"))] # More women than men 
pheno[, white_british := ifelse(white_british == "White", TRUE, FALSE)]
pheno[, visit := factor(visit, levels=c("Main Phase", "Repeat Assessment"))]

# Filter to first measurement per participant
pheno[, eid_7439 := as.integer(eid_7439)]
pheno <- pheno[first[, .(eid_7439, visit)], on = .(eid_7439, visit), nomatch=0]
pheno[, visit := factor(visit, levels=c("Main Phase", "Repeat Assessment"))]

# Join with extracted PCs
pc_gt1[, eid_7439 := as.integer(eid_7439)]
pheno <- pheno[pc_gt1, on = .(eid_7439), nomatch=0]

# Compute correlation between PCs and phenotypes
pc_cor <- foreach(dIdx = unique(pheno$dataset), .combine=rbind) %do% {
  foreach(pIdx = c("age", "sex", "bmi", "lipid_med", "fasting_time", "visit", "white_british"), .combine=rbind) %do% {
    foreach(pcIdx = sprintf("PC%s", 1:npcs), .combine=rbind) %do% {
      this_dat <- pheno[dataset == dIdx, .SD, .SDcols=c(pIdx, pcIdx)]
      if (!is.numeric(this_dat[[1]])) this_dat[[1]] <- as.integer(this_dat[[1]])
      this_cor <- cor(this_dat, method="spearman", use="pairwise.complete.obs")[1,2]
      data.table(dataset=dIdx, phenotype=pIdx, pc=pcIdx, spearman=this_cor)
    }
  }
}


# Plot correlation heatmap
pc_cor[, label := format(round(spearman*100)/100)]
pc_cor[abs(spearman) < 0.1, spearman := 0]

pal <- rev(brewer.pal(name="RdBu", n=11))
pal[6] <- "#FFFFFF"
cols <- colorRampPalette(pal)(256)

pc_cor[, pc := factor(pc, levels=sprintf("PC%s", 13:1))]
pc_cor[, phenotype := factor(phenotype,levels= c("age", "sex", "bmi", "lipid_med", "fasting_time", "visit", "white_british"))]

# Plot for post-QC data, drop PC13 because only 12 explain >1% variation in postQC data (13 in raw data)
g <- ggplot(pc_cor[dataset == "adj" & pc != "PC13"]) +
  aes(x=phenotype, y=pc, fill=spearman) +
  geom_raster() +
  geom_text(aes(label=label), size=1.8) +
  scale_x_discrete(name="", expand=c(0,0)) +
  scale_y_discrete(name="", expand=c(0,0)) +
  scale_fill_gradientn(name="Spearman correlation", colours=cols, limits=c(-1, 1)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=6), axis.text.x=element_text(size=6, angle=90, hjust=0, vjust=1),
        legend.title=element_text(size=7), legend.text=element_text(size=7), legend.key.height=unit(0.4, "cm"),
        plot.background=element_rect(size=0.5)
  ) + guides(fill = guide_colorbar(title.position="top", title.hjust=0.5))
ggsave(width=3.6, height=4.2, units="in", file="paper_output/postqc_pc_pheno_correlations.pdf")

g <- ggplot(pc_cor[dataset == "raw"]) +
  aes(x=phenotype, y=pc, fill=spearman) +
  geom_raster() +
  geom_text(aes(label=label), size=1.8) +
  scale_x_discrete(name="", expand=c(0,0)) +
  scale_y_discrete(name="", expand=c(0,0)) +
  scale_fill_gradientn(name="Spearman correlation", colours=cols, limits=c(-1, 1)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=6), axis.text.x=element_text(size=6, angle=90, hjust=0, vjust=1),
        legend.title=element_text(size=7), legend.text=element_text(size=7), legend.key.height=unit(0.4, "cm"),
        plot.background=element_rect(size=0.5)
  ) + guides(fill = guide_colorbar(title.position="top", title.hjust=0.5))
ggsave(width=3.6, height=4.2, units="in", file="paper_output/raw_pc_pheno_correlations.pdf")

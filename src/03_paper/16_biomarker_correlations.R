library(data.table)
library(ukbnmr)
library(ggplot2)
library(ggrastr)
library(WGCNA)

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

# Convert to matrices
raw <- as.matrix(dcast(dat[!is.na(raw)], eid_7439 ~ variable, value.var="raw"), rownames="eid_7439") 
adj <- as.matrix(dcast(dat[!is.na(adj5)], eid_7439 ~ variable, value.var="adj5"), rownames="eid_7439")

# Set non-finite values to missing (i.e. 0/0 or x/0 ratios)
adj[!is.na(adj) & !is.finite(adj)] <- NA

# Compute correlation coefficients.
raw <- cor(raw, use="pairwise.complete.obs", method="pearson")
adj <- cor(adj, use="pairwise.complete.obs", method="pearson")

# Get distance using TOMdist.
# This computes distance between biomarkers as a function of both their
# correlation with each other, and similarity of correlation to all other
# biomarkers
dist <- TOMdist(abs(adj)^6, TOMType="unsigned")
clust <- hclust(as.dist(dist), method="average")
biomarker_order <- colnames(adj)[clust$order]
raw_order <- intersect(biomarker_order, colnames(raw))

raw <- raw[raw_order, raw_order]
adj <- adj[biomarker_order, biomarker_order]

# Melt to long for plotting
adj <- as.data.table(adj, keep.rownames="rn")
adj <- melt(adj, id.vars="rn", variable.name="cn", value.name="pearson")
adj[,rn := factor(rn, levels=biomarker_order)]
adj[,cn := factor(cn, levels=biomarker_order)]

g <- ggplot(adj) +
  aes(x=rn, y=cn, fill=pearson) + 
  rasterize(geom_raster(), dpi=1200) +
  scale_x_discrete(name="", expand=expansion(add = 1.3)) +
  scale_y_discrete(name="", expand=expansion(add = 1.3)) +
  scale_fill_gradient2(name="Pearson correlation", low="#313695", mid="white", high="#a50026", limits=c(-1,1)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="bottom", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=1), axis.text.x=element_text(size=1, angle=90, hjust=1, vjust=0.5),
        legend.title=element_text(size=7), legend.text=element_text(size=6), legend.key.height=unit(0.4, "cm"),
        plot.background=element_rect(size=0.5)
  ) + guides(fill = guide_colorbar(title.position="top", title.hjust=0.5))
ggsave(width=3.6, height=4, units="in", file="paper_output/techqc_biomarker_correlations.pdf")

# Show raw biomarkers in same order
raw <- as.data.table(raw, keep.rownames="rn")
raw <- melt(raw, id.vars="rn", variable.name="cn", value.name="pearson")
raw[,rn := factor(rn, levels=raw_order)]
raw[,cn := factor(cn, levels=raw_order)]

g <- ggplot(raw) +
  aes(x=rn, y=cn, fill=pearson) + 
  rasterize(geom_raster(), dpi=1200) +
  scale_x_discrete(name="", expand=expansion(add = 1.3)) +
  scale_y_discrete(name="", expand=expansion(add = 1.3)) +
  scale_fill_gradient2(name="Pearson correlation", low="#313695", mid="white", high="#a50026", limits=c(-1,1)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="bottom", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=1), axis.text.x=element_text(size=1, angle=90, hjust=1, vjust=0.5),
        legend.title=element_text(size=7), legend.text=element_text(size=6), legend.key.height=unit(0.4, "cm"),
        plot.background=element_rect(size=0.5)
  ) + guides(fill = guide_colorbar(title.position="top", title.hjust=0.5))
ggsave(width=3.6, height=4, units="in", file="paper_output/raw_biomarker_correlations.pdf")


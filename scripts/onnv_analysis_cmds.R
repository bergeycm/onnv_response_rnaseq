#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Do differential expression analysis using edgeR
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

library(Rsubread)

# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)

library(RColorBrewer)
# BiocManager::install("mixOmics")
library(mixOmics)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(scales)
library(reshape2)
library(gprofiler2)

# ----------------------------------------------------------------------------------------
# --- Tally read counts using rSubread
# ----------------------------------------------------------------------------------------

bams = list.files(path = "final_bams", pattern = "*.bam$", full.names=TRUE)

# This includes all transcripts. We could run it just for known_transcripts file
# but we instead simply get rid of the XLOC genes later

gtf.file = "genomes/AgamP4.gtf"

fc = featureCounts(bams, annot.ext=gtf.file,
                   isGTFAnnotationFile=TRUE,
                   isPairedEnd=TRUE,
                   strandSpecific=2,
                   #useMetaFeatures=TRUE,
                   nthreads=20,
                   allowMultiOverlap=TRUE)

save.image("fc.Rdata")

write.table(data.frame(fc$annotation$GeneID, fc$counts),
    file="all.readCounts.txt", quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Create edgeR object
# ----------------------------------------------------------------------------------------

# load("fc.Rdata")

samp.info = data.frame(do.call(rbind,
    strsplit(
            gsub(".*bams\\.(.*)\\.final\\.bam", "\\1", colnames(fc$counts)),
        split="_")))
names(samp.info) = c("replicate", "dpi", "sex", "infection.status")
samp.groups = apply(samp.info[,c("sex", "dpi", "infection.status")], 1,
    function (x) { paste(x, collapse="_")})

d = DGEList(counts=fc$counts, genes=fc$annotation, group=samp.groups)
dim(d)
# [1] 13796    24

# Make sample names easier to read
colnames(d$counts)   = gsub(".*bams\\.(.*)\\.final.*", "\\1", colnames(d$counts))
row.names(d$samples) = gsub(".*bams\\.(.*)\\.final.*", "\\1", row.names(d$samples))

d.full = d
# head(d$counts)
# head(cpm(d)) # cpm per gene

# --- Filter for low expression genes and normalize

do.filtering.normalization = function (d) {

    keep = rowSums(cpm(d)>5) >= 2
    d = d[keep,]
    dim(d)

    d$samples$lib.size = colSums(d$counts)
    d$samples
    d = calcNormFactors(d)
    d
}

d = do.filtering.normalization(d)

# --- Do MDS plots

do.mds.plots = function(d, samp.info, suffix="") {

    # Color by sex
    pdf(paste0("reports/MDS_color_by_sex", suffix, ".pdf"))
    p1 = plotMDS(d, method="bcv", col=as.numeric(factor(samp.info$sex)))
    #print(p1)
    dev.off()

    # Color by dpi
    pdf(paste0("reports/MDS_color_by_dpi", suffix, ".pdf"))
    p2 = plotMDS(d, method="bcv", col=as.numeric(factor(samp.info$dpi)))
    #print(p2)
    dev.off()

    # Color by infection status
    pdf(paste0("reports/MDS_color_by_infection", suffix, ".pdf"))
    p3 = plotMDS(d, method="bcv", col=as.numeric(factor(samp.info$infection.status)))
    #print(p3)
    dev.off()
}

do.mds.plots(d, samp.info)

# --- Exclude outlier low-coverage sample and redo MDS plots

d = d.full

outlier = "TX2_D7_M_neg"

d$counts  = d$counts[,colnames(d$counts) != outlier]
d$samples = d$samples[rownames(d$samples) != outlier,]

outlier.idx = which(outlier == apply(samp.info, 1, function(x) {
    paste(x, sep="_", collapse="_")
}))
samp.info = samp.info[-outlier.idx,]

samp.groups = samp.groups[-outlier.idx]

d = do.filtering.normalization(d)

do.mds.plots(d, samp.info, suffix=".sans_outlier")

# ----------------------------------------------------------------------------------------

# --- Set up design and contrasts and estimate dispersion

samp.groups = factor(samp.groups)

design = model.matrix(~ 0 + samp.groups)

colnames(design) = levels(samp.groups)

# Passing it just counts works too, for future reference
# d.disp.counts = estimateDisp(d$counts, design)
d.disp = estimateDisp(d, design)

fit = glmQLFit(d.disp, design)

my.contrasts = makeContrasts(
    # Simple binary comparisons with one grouping - Sex
    sex.FvM = (F_D2_neg + F_D2_plus + F_D7_neg + F_D7_plus) / 4 -
              (M_D2_neg + M_D2_plus + M_D7_neg + M_D7_plus) / 4,
    sex.FvM.onlyD2 =
              (F_D2_neg + F_D2_plus) / 2 -
              (M_D2_neg + M_D2_plus) / 2,
    sex.FvM.onlyD7 =
              (F_D7_neg + F_D7_plus) / 2 -
              (M_D7_neg + M_D7_plus) / 2,

    # Simple binary comparisons with one grouping - Infection status
    inf.FvT =        (F_D2_neg  + F_D7_neg  + M_D2_neg  + M_D7_neg ) / 4 -
                     (F_D2_plus + F_D7_plus + M_D2_plus + M_D7_plus) / 4,
    inf.FvT.onlyD2 = (F_D2_neg  + M_D2_neg ) / 2 -
                     (F_D2_plus + M_D2_plus) / 2,
    inf.FvT.onlyD7 = (F_D7_neg  + M_D7_neg ) / 2 -
                     (F_D7_plus + M_D7_plus) / 2,
    inf.FvT.onlyF =  (F_D2_neg  + F_D7_neg ) / 2 -
                     (F_D2_plus + F_D7_plus) / 2,
    inf.FvT.onlyM =  (M_D2_neg  + M_D7_neg ) / 2 -
                     (M_D2_plus + M_D7_plus) / 2,

    # Simple binary comparisons with one grouping - DPI
    dpi.2v7 =       (F_D2_neg + F_D2_plus + M_D2_neg + M_D2_plus) / 4 -
                    (F_D7_neg + F_D7_plus + M_D7_neg + M_D7_plus) / 4,
    dpi.2v7.onlyF = (F_D2_neg + F_D2_plus) / 2 -
                    (F_D7_neg + F_D7_plus) / 2,
    dpi.2v7.onlyM = (M_D2_neg + M_D2_plus) / 2 -
                    (M_D7_neg + M_D7_plus) / 2,

    # Between-sex differences by time
    between.sex.dpi =
              ((F_D2_neg + F_D2_plus) / 2 - (F_D7_neg + F_D7_plus) / 2) / 2 -
              ((M_D2_neg + M_D2_plus) / 2 - (M_D7_neg + M_D7_plus) / 2) / 2,

    # Between-sex differences by infection status
    between.sex.inf =
              ((F_D2_neg + F_D7_neg) / 2 - (F_D2_plus + F_D7_plus) / 2) / 2 -
              ((M_D2_neg + M_D7_neg) / 2 - (M_D2_plus + M_D7_plus) / 2) / 2,
    between.sex.inf.onlyD2 =
              (F_D2_neg - F_D2_plus) / 2 -
              (M_D2_neg - M_D2_plus) / 2,
    between.sex.inf.onlyD7 =
              (F_D7_neg - F_D7_plus) / 2 -
              (M_D7_neg - M_D7_plus) / 2,
    levels = design)

#    # Levels:
#    # F_D2_neg F_D2_plus F_D7_neg F_D7_plus M_D2_neg M_D2_plus M_D7_neg M_D7_plus

contrast.list = colnames(my.contrasts)

qlf = list()

do.qlf.test = function (contrast) {

    glmQLFTest(fit, contrast=my.contrasts[,contrast])
}

qlf = lapply(contrast.list, do.qlf.test)
names(qlf) = contrast.list

# ----------------------------------------------------------------------------------------
# --- Make table of top DE genes
# ----------------------------------------------------------------------------------------

# Access top DE genes, with, e.g.:
#    topTags(qlf[["sex.FvM"]])

write.de.table = function (contrast, fdr.cutoff=0.1) {

    top.de = topTags(qlf[[contrast]], n = 1000,
        adjust.method = "fdr", sort.by = "PValue", p.value = fdr.cutoff)
    to.print = top.de$table[,c("GeneID", "logFC", "logCPM", "PValue", "FDR")]

    write.table(to.print, file=paste0("reports/top_DE_genes-", contrast, ".txt"),
        sep="\t", quote=FALSE)

    return(top.de$table)
}

fdr.cutoffs = list()

fdr.cutoffs[["sex.FvM"]]        = 1e-18
fdr.cutoffs[["sex.FvM.onlyD2"]] = 1e-18
fdr.cutoffs[["sex.FvM.onlyD7"]] = 1e-18

fdr.cutoffs[["inf.FvT"]]        = 1e-3
fdr.cutoffs[["inf.FvT.onlyD2"]] = 1e-3
fdr.cutoffs[["inf.FvT.onlyD7"]] = 0.1   # Reduced
fdr.cutoffs[["inf.FvT.onlyF"]]  = 1e-3
fdr.cutoffs[["inf.FvT.onlyM"]]  = 0.1   # Reduced

fdr.cutoffs[["dpi.2v7"]]       = 1e-3
fdr.cutoffs[["dpi.2v7.onlyF"]] = 1e-3
fdr.cutoffs[["dpi.2v7.onlyM"]] = 1e-3

fdr.cutoffs[["between.sex.dpi"]] = 1e-3

fdr.cutoffs[["between.sex.inf"]]        = 0.1  # Reduced
fdr.cutoffs[["between.sex.inf.onlyD2"]] = 1e-2
fdr.cutoffs[["between.sex.inf.onlyD7"]] = 0.5  # Reduced A LOT

all.de = lapply(contrast.list, function(x) { write.de.table(x, fdr.cutoffs[[x]]) })
names(all.de) = contrast.list

# Get count of genes in tables
lapply(all.de, dim)

# Get counts with, e.g.:
#    summary(decideTests((qlf[["sex.FvM"]])))

# ----------------------------------------------------------------------------------------
# --- Create mean-difference (MA or MD) plot
# ----------------------------------------------------------------------------------------

do.md.plot = function (contrast, fdr.cutoff=0.1, plot.title) {

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    pdf(paste0("reports/MD_plot-", contrast, ".pdf"))

    plotMD(qlf[[contrast]],
        status = adj.p < fdr.cutoff,
        main = paste0("Mean-difference - ", plot.title))
    abline(h=c(-1,1), col="black", lty=2) # Line = 2x FC
    dev.off()

    return(NA)
}

plot.titles = list()
plot.titles[["sex.FvM"]]        = "Sex (Female vs. Male)"
plot.titles[["sex.FvM.onlyD2"]] = "Sex (Female vs. Male) [2 dpi]"
plot.titles[["sex.FvM.onlyD7"]] = "Sex (Female vs. Male) [7 dpi]"

plot.titles[["inf.FvT"]]        = "Infection status (Control vs. Infected)"
plot.titles[["inf.FvT.onlyD2"]] = "Infection status (Control vs. Infected) [2 dpi]"
plot.titles[["inf.FvT.onlyD7"]] = "Infection status (Control vs. Infected) [7 dpi]"
plot.titles[["inf.FvT.onlyF"]]  = "Infection status (Control vs. Infected) [females]"
plot.titles[["inf.FvT.onlyM"]]  = "Infection status (Control vs. Infected) [males]"

plot.titles[["dpi.2v7"]] = "DPI (2 vs 7 days)"
plot.titles[["dpi.2v7.onlyF"]] = "DPI (2 vs 7 days) [females]"
plot.titles[["dpi.2v7.onlyM"]] = "DPI (2 vs 7 days) [males]"

plot.titles[["between.sex.dpi"]] = "Sex-based response to time (Female vs. Male)"

plot.titles[["between.sex.inf"]] =
                            "Sex-based response to infection (Female vs. Male)"
plot.titles[["between.sex.inf.onlyD2"]] =
                            "Sex-based response to infection (Female vs. Male) [2 dpi]"
plot.titles[["between.sex.inf.onlyD7"]] =
                            "Sex-based response to infection (Female vs. Male) [7 dpi]"

tmp = lapply(contrast.list, function(x) {
    do.md.plot(x, fdr.cutoffs[[x]], plot.titles[[x]])
})

# ----------------------------------------------------------------------------------------
# --- Create volcano plots
# ----------------------------------------------------------------------------------------

do.volcano.plot = function (contrast, fdr.cutoff=0.1, plot.title) {

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    pdf(paste0("reports/volcano_plot-", contrast, ".pdf"))

    volcano.dat = cbind(qlf[[contrast]]$table$logFC,
                        -log10(adj.p))
    colnames(volcano.dat) = c("logFC", "negLogP")

    plot(volcano.dat, pch=19, cex=0.5,
        col=c("black", "red")[factor(adj.p < fdr.cutoff)])

    dev.off()

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    do.volcano.plot(x, fdr.cutoffs[[x]], plot.titles[[x]])
})

# Consider EnhancedVolcano BioConductor package

# ----------------------------------------------------------------------------------------
# --- DE genes: Expression profile heat map
# ----------------------------------------------------------------------------------------

do.heatmap.plot = function (contrast, fdr.cutoff=0.1, plot.title) {

    write(paste0("Heatmap for contrast [", contrast,
            "]"), stderr())

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    de.genes = rownames(qlf[[contrast]]$table)[adj.p < fdr.cutoff &
                        abs(qlf[[contrast]]$table$logFC) > 1.5]

    d.cpm = cpm(d, log=TRUE, prior.count = 1)
    de.cpm = d.cpm[which(rownames(d.cpm) %in% de.genes),]

    if (length(de.genes) == 0) {
        write(paste0("No DE genes for contrast [", contrast,
            "]. Skipping heatmap."), stderr())
        return(NA)
    }

    # Stop this from being a vector instead of a matrix
    if (length(de.genes) < 2) {
        de.cpm = t(as.matrix(de.cpm))
        row.names(de.cpm) = de.genes

        write(paste0("Clustering requires >= 2 DE genes. ",
            "Skipping heatmap for contrast [", contrast, "]."), stderr())
        return(NA)
    }

    pal = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]

    pdf(paste0("reports/heatmap-", contrast, ".pdf"))
    cim(t(de.cpm), color=pal, symkey=FALSE)
    dev.off()

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    do.heatmap.plot(x, fdr.cutoffs[[x]], plot.titles[[x]])
})

# ----------------------------------------------------------------------------------------
# --- DE genes: PCA
# ----------------------------------------------------------------------------------------

do.pca.plot = function (contrast, fdr.cutoff=0.1, plot.title) {

    write(paste0("Doing PCA plot for contrast [", contrast, "]."), stderr())

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    de.genes = rownames(qlf[[contrast]]$table)[
                            adj.p < fdr.cutoff &
                            abs(qlf[[contrast]]$table$logFC) > 1.5]

    if (length(de.genes) <= 2) {
        write(paste0("Insufficient DE genes for contrast [", contrast,
            "]. Skipping PCA."), stderr())
        return(NA)
    }

    d.cpm = cpm(d, log=TRUE, prior.count = 1)
    de.cpm = d.cpm[which(rownames(d.cpm) %in% de.genes),]

    df.pca = de.cpm

    res.pca = prcomp(df.pca, scale = FALSE)

    pca.labels = rownames(res.pca$x)

    x.lims.exp = expand_range(range(res.pca$x[,1]), mul=0.2)
    y.lims.exp = expand_range(range(res.pca$x[,2]), mul=0.2)

    de.gene.info = qlf[[contrast]]$table[pca.labels,]

    # Note colored by p-value, not adjusted p-value
    p = ggplot(res.pca$x, aes(PC1, PC2, label=pca.labels,
            color=-1 * de.gene.info$logFC, size=-1 * log(de.gene.info$PValue))) +
        geom_point() +
        geom_text_repel(size=2, segment.size=0.2, color="black") +
        coord_fixed() +
        xlim(x.lims.exp) +
        ylim(y.lims.exp) +
        ggtitle(paste0(plot.title, " - DE genes (adjusted p < ", fdr.cutoff, ")")) +
        scale_colour_gradient2(name="Log(FC)") +
        scale_size_continuous(
            name = "Log(p-value)",
            range = c(0.5, 2),
        ) +
        theme_bw()

    ggsave(p, file=paste0("reports/PCA_DE_genes-", contrast, ".pdf"))

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    do.pca.plot(x, fdr.cutoffs[[x]], plot.titles[[x]])
})

# ----------------------------------------------------------------------------------------
# --- DE genes: Histograms of Fold Change
# ----------------------------------------------------------------------------------------

do.fc_histogram.plot = function (contrast, fdr.cutoff=0.1, plot.title) {

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    gene.info = qlf[[contrast]]$table
    gene.info$is.DE = FALSE

    if (sum(adj.p < fdr.cutoff) == 0) {
        write(paste0("Insufficient DE genes for contrast [", contrast,
            "]. Skipping Histogram of fold change"), stderr())
        return(NA)
    }

    # Note filtering only based on FDR for this, not log(FC) too
    gene.info[adj.p < fdr.cutoff,]$is.DE = TRUE

    p = ggplot(gene.info, aes(x=logFC, color=is.DE, fill=is.DE)) +
        geom_histogram(alpha=0.5, bins=100) +
        facet_grid(is.DE ~ ., scales="free_y") +
        ggtitle(paste0(plot.title, " - DE genes (adjusted p < ", fdr.cutoff, ")")) +
        scale_color_brewer(palette="Dark2") +
        scale_fill_brewer(palette="Dark2") +
        theme_bw()

    ggsave(p, file=paste0("reports/FC_histogram_DE_genes-", contrast, ".pdf"))

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    do.fc_histogram.plot(x, fdr.cutoffs[[x]], plot.titles[[x]])
})

# ----------------------------------------------------------------------------------------
# --- DE genes: Expression stripcharts
# ----------------------------------------------------------------------------------------

do.expr_stripchart.plot = function (contrast, fdr.cutoff=0.1, plot.title,
    first.var, second.var, third.var) {

    de.gene.info = topTags(qlf[[contrast]], n=10, adjust.method = "fdr")

    de.genes = rownames(de.gene.info)

    d.cpm = cpm(d, log=TRUE, prior.count = 1)
    de.cpm = d.cpm[which(rownames(d.cpm) %in% de.genes),]

    de.cpm.t = melt(t(de.cpm), id.vars=colnames(t(de.cpm)))
    names(de.cpm.t) = c("sample", "gene", "CPM")

    de.cpm.t = cbind(de.cpm.t, samp.info)

    p = ggplot(de.cpm.t, aes_string(x=first.var, y="CPM",
            color=second.var, fill=second.var, pch=third.var)) +
        geom_point(alpha=0.5) +
        geom_jitter(width=0.2) +
        facet_grid(gene ~ ., scales="free_y") +
        labs(title=plot.title, subtitle=paste0("DE genes (adjusted p < ", fdr.cutoff, ")")) +
        scale_color_brewer(palette="Dark2") +
        scale_fill_brewer(palette="Dark2") +
        scale_shape_manual(values=c(1, 16)) +
        theme_bw() +
        theme(strip.text = element_text(size = 6))

    ggsave(p, file=paste0("reports/stripplot_DE_genes-", contrast, ".pdf"),
        width=4, height=7)

    return(NA)
}

# --- Sex

do.expr_stripchart.plot("sex.FvM", fdr.cutoffs[["sex.FvM"]],
    plot.titles[["sex.FvM"]],
    "sex", "infection.status", "dpi")

do.expr_stripchart.plot("sex.FvM.onlyD2", fdr.cutoffs[["sex.FvM.onlyD2"]],
    plot.titles[["sex.FvM.onlyD2"]],
    "sex", "infection.status", "sex")

do.expr_stripchart.plot("sex.FvM.onlyD7", fdr.cutoffs[["sex.FvM.onlyD7"]],
    plot.titles[["sex.FvM.onlyD7"]],
    "sex", "infection.status", "sex")

# --- Infection status

do.expr_stripchart.plot("inf.FvT", fdr.cutoffs[["inf.FvT"]],
    plot.titles[["inf.FvT"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("inf.FvT.onlyD2", fdr.cutoffs[["inf.FvT.onlyD2"]],
    plot.titles[["inf.FvT.onlyD2"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("inf.FvT.onlyD7", fdr.cutoffs[["inf.FvT.onlyD7"]],
    plot.titles[["inf.FvT.onlyD7"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("inf.FvT.onlyF", fdr.cutoffs[["inf.FvT.onlyF"]],
    plot.titles[["inf.FvT.onlyF"]],
    "infection.status", "dpi", "sex")

do.expr_stripchart.plot("inf.FvT.onlyM", fdr.cutoffs[["inf.FvT.onlyM"]],
    plot.titles[["inf.FvT.onlyM"]],
    "infection.status", "dpi", "sex")

# --- DPI

do.expr_stripchart.plot("dpi.2v7", fdr.cutoffs[["dpi.2v7"]],
    plot.titles[["dpi.2v7"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("dpi.2v7.onlyF", fdr.cutoffs[["dpi.2v7.onlyF"]],
    plot.titles[["dpi.2v7.onlyF"]],
    "infection.status", "dpi", "sex")

do.expr_stripchart.plot("dpi.2v7.onlyM", fdr.cutoffs[["dpi.2v7.onlyM"]],
    plot.titles[["dpi.2v7.onlyM"]],
    "infection.status", "dpi", "sex")

# --- More complicated contrasts

do.expr_stripchart.plot("between.sex.dpi", fdr.cutoffs[["between.sex.dpi"]],
    plot.titles[["between.sex.dpi"]],
    "sex", "dpi", "infection.status")

do.expr_stripchart.plot("between.sex.inf", fdr.cutoffs[["between.sex.inf"]],
    plot.titles[["between.sex.inf"]],
    "sex", "infection.status", "dpi")

do.expr_stripchart.plot("between.sex.inf.onlyD2", fdr.cutoffs[["between.sex.inf.onlyD2"]],
    plot.titles[["between.sex.inf.onlyD2"]],
    "sex", "infection.status", "dpi")

do.expr_stripchart.plot("between.sex.inf.onlyD7", fdr.cutoffs[["between.sex.inf.onlyD7"]],
    plot.titles[["between.sex.inf.onlyD7"]],
    "sex", "infection.status", "dpi")

# ----------------------------------------------------------------------------------------
# --- Do GO enrichment analysis
# ----------------------------------------------------------------------------------------

find.overrep.go = function (contrast, fdr.cutoff=0.1) {

    write(paste0("Doing GO enrichment analysis for contrast [",
        contrast, "]."), stderr())

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    de.genes.info = qlf[[contrast]]$table[
                            adj.p < fdr.cutoff &
                            abs(qlf[[contrast]]$table$logFC) > 1.5,]
    de.genes = rownames(de.genes.info)

    de.genes.up = rownames(de.genes.info[de.genes.info$logFC > 0,])
    de.genes.dn = rownames(de.genes.info[de.genes.info$logFC < 0,])

    if (length(de.genes.up) > 0) {
        gpro.up = gost(de.genes.up,
            organism = "agambiae", ordered_query = FALSE,
            multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
            measure_underrepresentation = FALSE, evcodes = FALSE,
            user_threshold = 0.2, correction_method = "fdr",
            domain_scope = "annotated",
            custom_bg = rownames(qlf[[contrast]]$table),
            numeric_ns="NO_NAMESPACE", sources = NULL, as_short_link = FALSE)
        gpro.up.simp = gpro.up$result
        gpro.up.simp = gpro.up.simp[gpro.up.simp$significant,
            c("p_value", "term_size", "query_size", "intersection_size",
                "source", "term_id", "term_name")]

        if ((nrow(gpro.up.simp) > 0) && (is.null(gpro.up.simp) == FALSE)) {
            gpro.up.simp = cbind(contrast, gpro.up.simp)
            write.table(gpro.up.simp, file=paste0("reports/GO_overrep-", contrast, ".up.txt"),
                sep="\t", row.names=FALSE)
        } else {
            write(paste0("Insufficient significant upregulated GO terms for contrast [",
                contrast,
                "]. No output file created"), stderr())
        }
    } else {
        write(paste0("Insufficient significant upregulated GO terms for contrast [",
            contrast,
            "]. No output file created"), stderr())
    }

    if (length(de.genes.dn) > 0) {
        gpro.dn = gost(de.genes.dn,
            organism = "agambiae", ordered_query = FALSE,
            multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
            measure_underrepresentation = FALSE, evcodes = FALSE,
            user_threshold = 0.2, correction_method = "fdr",
            domain_scope = "annotated",
            custom_bg = rownames(qlf[[contrast]]$table),
            numeric_ns="NO_NAMESPACE", sources = NULL, as_short_link = FALSE)
        gpro.dn.simp = gpro.dn$result
        gpro.dn.simp = gpro.dn.simp[gpro.up.simp$significant,
            c("p_value", "term_size", "query_size", "intersection_size",
                "source", "term_id", "term_name")]

        if ((nrow(gpro.dn.simp) > 0) && (is.null(gpro.dn.simp) == FALSE)) {
            gpro.dn.simp = cbind(contrast, gpro.dn.simp)
            write.table(gpro.dn.simp,
                file=paste0("reports/GO_overrep-", contrast, ".down.txt"),
                sep="\t", row.names=FALSE)
        } else {
            write(paste0("Insufficient significant downregulated GO terms for contrast [",
                contrast,
                "]. No output file created"), stderr())
        }
    } else {
            write(paste0("Insufficient significant downregulated GO terms for contrast [",
                contrast,
                "]. No output file created"), stderr())
        }

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    find.overrep.go(x, fdr.cutoffs[[x]])
    Sys.sleep(time=3)
})

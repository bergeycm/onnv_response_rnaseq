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

    mds = plotMDS(d, method="bcv")
    mds.df = data.frame(names(mds$x), mds$x, mds$y, samp.info)

    mds.df$label = ""
    if ("TX2_D7_M_neg" %in% mds.df$names.mds.x) {
        mds.df[mds.df$names.mds.x. == "TX2_D7_M_neg",]$label = "TX2_D7_M_neg"
    }

    # Color by sex

    p = ggplot(mds.df, aes(mds.x, mds.y,
            col=sex, shape=infection.status, label=label)) +
        geom_point(cex=3) +
        geom_text_repel() +
        xlab("BCV Distance 1") +
        ylab("BCV Distance 2") +
        theme_bw() +
        scale_shape_manual(
            values = c(16, 2),
            breaks = c("plus", "neg"),
            labels = c("Infected", "Control"),
            name   = "Infection\nStatus:",
        ) +
        scale_color_manual(
            values = c("#297A45", "#923158"),
            breaks = c("F", "M"),
            labels = c("Female", "Male"),
            name   = "Sex:",
        )
    ggsave(p, file=paste0("reports/MDS_color_by_sex", suffix, ".pdf"),
        height=6, width=6)

    # Color by dpi

    p = ggplot(mds.df, aes(mds.x, mds.y,
            fill=dpi, label=label)) +
        geom_point(cex=3, shape=21) +
        geom_text_repel() +
        xlab("BCV Distance 1") +
        ylab("BCV Distance 2") +
        theme_bw() +
        scale_fill_manual(
            values = c("#D4A26A", "#552D00"),
            breaks = c("D2", "D7"),
            labels = c("2 dpi", "7 dpi"),
            name   = "Time point:",
        )
    ggsave(p, file=paste0("reports/MDS_color_by_dpi", suffix, ".pdf"),
        height=6, width=6)

    # Color by infection status

    p = ggplot(mds.df, aes(mds.x, mds.y,
            fill=infection.status, label=label)) +
        geom_point(cex=3, shape=21) +
        geom_text_repel() +
        xlab("BCV Distance 1") +
        ylab("BCV Distance 2") +
        theme_bw() +
        scale_fill_manual(
            values = c("#802A15", "#FFFFFF"),
            breaks = c("plus", "neg"),
            labels = c("Infected", "Control"),
            name   = "Infection\nStatus:",
        )
    ggsave(p, file=paste0("reports/MDS_color_by_infection", suffix, ".pdf"),
        height=6, width=6)
}

do.mds.plots(d, samp.info)

# --- Exclude outlier low-coverage sample and redo MDS plots

d = d.full

outlier = "TX2_D7_M_neg"

d$counts  = d$counts[,colnames(d$counts)  != outlier]
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
    inf.TvF =        (F_D2_plus + F_D7_plus + M_D2_plus + M_D7_plus ) / 4 -
                     (F_D2_neg  + F_D7_neg  + M_D2_neg  + M_D7_neg) / 4,
    inf.TvF.onlyD2 = (F_D2_plus + M_D2_plus ) / 2 -
                     (F_D2_neg  + M_D2_neg) / 2,
    inf.TvF.onlyD7 = (F_D7_plus + M_D7_plus ) / 2 -
                     (F_D7_neg  + M_D7_neg) / 2,
    inf.TvF.onlyF =  (F_D2_plus + F_D7_plus ) / 2 -
                     (F_D2_neg  + F_D7_neg) / 2,
    inf.TvF.onlyM =  (M_D2_plus + M_D7_plus ) / 2 -
                     (M_D2_neg  + M_D7_neg) / 2,

    # Simple binary comparisons with one grouping - DPI
    dpi.7v2 =       (F_D2_plus + F_D2_neg + M_D2_plus + M_D2_neg) / 4 -
                    (F_D7_plus + F_D7_neg + M_D7_plus + M_D7_neg) / 4,
    dpi.7v2.onlyF = (F_D2_plus + F_D2_neg) / 2 -
                    (F_D7_plus + F_D7_neg) / 2,
    dpi.7v2.onlyM = (M_D2_plus + M_D2_neg) / 2 -
                    (M_D7_plus + M_D7_neg) / 2,

    # Between-sex differences by time
    between.sex.dpi =
              ((F_D7_neg + F_D7_plus) / 2 - (F_D2_neg + F_D2_plus) / 2) / 2 -
              ((M_D7_neg + M_D7_plus) / 2 - (M_D2_neg + M_D2_plus) / 2) / 2,

    # Between-sex differences by infection status
    between.sex.inf =
              ((F_D2_plus + F_D7_plus) / 2 - (F_D2_neg + F_D7_neg) / 2) / 2 -
              ((M_D2_plus + M_D7_plus) / 2 - (M_D2_neg + M_D7_neg) / 2) / 2,
    between.sex.inf.onlyD2 =
              (F_D2_plus - F_D2_neg) / 2 -
              (M_D2_plus - M_D2_neg) / 2,
    between.sex.inf.onlyD7 =
              (F_D7_plus - F_D7_neg) / 2 -
              (M_D7_plus - M_D7_neg) / 2,
    levels = design)

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

fdr.cutoffs[["inf.TvF"]]        = 1e-3
fdr.cutoffs[["inf.TvF.onlyD2"]] = 1e-3
fdr.cutoffs[["inf.TvF.onlyD7"]] = 0.1   # Inceased
fdr.cutoffs[["inf.TvF.onlyF"]]  = 1e-3
fdr.cutoffs[["inf.TvF.onlyM"]]  = 0.1   # Inceased

fdr.cutoffs[["dpi.7v2"]]       = 1e-3
fdr.cutoffs[["dpi.7v2.onlyF"]] = 1e-3
fdr.cutoffs[["dpi.7v2.onlyM"]] = 1e-3

fdr.cutoffs[["between.sex.dpi"]] = 1e-3

fdr.cutoffs[["between.sex.inf"]]        = 0.1  # Inceased
fdr.cutoffs[["between.sex.inf.onlyD2"]] = 1e-2
fdr.cutoffs[["between.sex.inf.onlyD7"]] = 0.5  # Inceased A LOT

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

plot.titles[["inf.TvF"]]        = "Infection status (Control vs. Infected)"
plot.titles[["inf.TvF.onlyD2"]] = "Infection status (Control vs. Infected) [2 dpi]"
plot.titles[["inf.TvF.onlyD7"]] = "Infection status (Control vs. Infected) [7 dpi]"
plot.titles[["inf.TvF.onlyF"]]  = "Infection status (Control vs. Infected) [females]"
plot.titles[["inf.TvF.onlyM"]]  = "Infection status (Control vs. Infected) [males]"

plot.titles[["dpi.7v2"]] = "DPI (2 vs 7 days)"
plot.titles[["dpi.7v2.onlyF"]] = "DPI (2 vs 7 days) [females]"
plot.titles[["dpi.7v2.onlyM"]] = "DPI (2 vs 7 days) [males]"

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

    volcano.dat = data.frame(
        gene    = rownames(qlf[[contrast]]$table),
        qlf[[contrast]]$table,
        adj.p   = adj.p,
        negLogP = -log10(adj.p)
    )

    volcan.pos = volcano.dat[volcano.dat$logFC > 0,]
    volcan.neg = volcano.dat[volcano.dat$logFC < 0,]

    volcan.pos = volcan.pos[order(volcan.pos$PValue),]
    volcan.neg = volcan.neg[order(volcan.neg$PValue),]

    volcano.dat$label = ""

    to.label = volcano.dat$gene %in% c(volcan.pos$gene[1:3], volcan.neg$gene[1:3])
    volcano.dat[to.label,]$label = row.names(volcano.dat)[to.label]

    p = ggplot(volcano.dat, aes(logFC, negLogP,
            color = adj.p < fdr.cutoff, label=label)) +
        geom_point(pch=19, cex=0.5) +
        geom_text_repel(size=3, col="black",
            nudge_y = max(volcano.dat$negLogP) / 10,
            segment.alpha = 0.1) +
        xlab("log(FC)") +
        ylab("-log(p)") +
        ylim(c(0, max(volcano.dat$negLogP) * 1.1)) +
        theme_bw() +
        scale_color_manual(
            values = c("#E07B39", "black"),
            breaks = c(TRUE, FALSE),
            labels = c("DE", "Not DE"),
            name   = "Differential\nExpression:",
        )

    ggsave(p, file=paste0("reports/volcano_plot-", contrast, ".pdf"))

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    do.volcano.plot(x, fdr.cutoffs[[x]], plot.titles[[x]])
})

# ----------------------------------------------------------------------------------------
# --- DE genes: Expression profile heat map
# ----------------------------------------------------------------------------------------

do.heatmap.plot = function (contrast, fdr.cutoff=0.1, plot.title,
                            plot.de = TRUE) {

    write(paste0("Heatmap for contrast [", contrast,
            "]"), stderr())

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    if (plot.de) {
        de.genes = rownames(qlf[[contrast]]$table)[adj.p < fdr.cutoff &
                            abs(qlf[[contrast]]$table$logFC) > 1.5]
    } else {
        all.genes = qlf[[contrast]]$table

        all.pos = all.genes[all.genes$logFC > 0,]
        all.neg = all.genes[all.genes$logFC < 0,]

        all.pos = all.pos[order(all.pos$PValue),]
        all.neg = all.neg[order(all.neg$PValue),]

        de.genes = c(rownames(head(all.pos, n=50)),
                     rownames(head(all.neg, n=50)))
    }

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

    pal = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)[255:1]

    pdf(paste0("reports/heatmap-", contrast, ".pdf"),
        width=10, height=6)
    cim(t(de.cpm), color=pal, symkey=FALSE,
        margins=c(15,15), row.cex=1, col.cex=0.5)
    dev.off()

    return(NA)
}

tmp = lapply(contrast.list, function(x) {
    do.heatmap.plot(x, fdr.cutoffs[[x]], plot.titles[[x]], plot.de=FALSE)
})

# ----------------------------------------------------------------------------------------
# --- PCA
# ----------------------------------------------------------------------------------------

do.pca.plot = function (contrast, fdr.cutoff=0.1, plot.title,
                        plot.de = TRUE) {

    write(paste0("Doing PCA plot for contrast [", contrast, "]."), stderr())

    adj.p = p.adjust(qlf[[contrast]]$table$PValue, method="fdr")

    if (plot.de) {
        de.genes = rownames(qlf[[contrast]]$table)[
                                adj.p < fdr.cutoff &
                                abs(qlf[[contrast]]$table$logFC) > 1.5]
    } else {
        all.genes = qlf[[contrast]]$table

        all.pos = all.genes[all.genes$logFC > 0,]
        all.neg = all.genes[all.genes$logFC < 0,]

        all.pos = all.pos[order(all.pos$PValue),]
        all.neg = all.neg[order(all.neg$PValue),]

        de.genes = c(rownames(head(all.pos, n=5)),
                     rownames(head(all.neg, n=5)))
    }

    if (length(de.genes) <= 2) {
        write(paste0("Insufficient DE genes for contrast [", contrast,
            "]. Skipping PCA."), stderr())
        return(NA)
    }

    d.cpm = cpm(d, log=TRUE, prior.count = 1)

    df.pca = d.cpm

    res.pca = prcomp(df.pca, scale = FALSE)

    to.label = rownames(res.pca$x) %in% de.genes
    pca.labels = rep("", nrow(res.pca$x))
    pca.labels[to.label] = rownames(res.pca$x)[to.label]

    x.lims.exp = expand_range(range(res.pca$x[,1]), mul=0.2)
    y.lims.exp = expand_range(range(res.pca$x[,2]), mul=0.2)

    # Note colored by p-value, not adjusted p-value
    p = ggplot(res.pca$x, aes(PC1, PC2, label=pca.labels,
            color=-1 * qlf[[contrast]]$table$logFC,
            size=-1 * log(qlf[[contrast]]$table$PValue))) +
        geom_point() +
        geom_text_repel(size=2, segment.size=0.2, segment.alpha=0.1, color="black") +
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
    do.pca.plot(x, fdr.cutoffs[[x]], plot.titles[[x]], plot.de = FALSE)
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

    de.gene.info = topTags(qlf[[contrast]], n=1000, adjust.method = "fdr")$table
    de.gene.info.pos = de.gene.info[de.gene.info$logFC > 0,]
    de.gene.info.neg = de.gene.info[de.gene.info$logFC < 0,]

    de.gene.info.pos = de.gene.info.pos[order(de.gene.info.pos$PValue),]
    de.gene.info.neg = de.gene.info.neg[order(de.gene.info.neg$PValue),]

    de.genes = c(rownames(de.gene.info.pos)[1:12],
                 rownames(de.gene.info.neg)[1:12])

    d.cpm = cpm(d, log=TRUE, prior.count = 1)
    de.cpm = d.cpm[which(rownames(d.cpm) %in% de.genes),]

    de.cpm.t = melt(t(de.cpm), id.vars=colnames(t(de.cpm)))
    names(de.cpm.t) = c("sample", "gene", "CPM")

    de.cpm.t = cbind(de.cpm.t, samp.info)

    de.cpm.t$FC.is.pos = de.cpm.t$gene %in% rownames(de.gene.info.pos)[1:12]

    de.cpm.t = merge(de.cpm.t, de.gene.info[,c(1,7,12)], by.x="gene", by.y="GeneID")
    de.cpm.t = de.cpm.t[order(de.cpm.t$PValue),]

    de.cpm.t$gene = factor(de.cpm.t$gene, levels = de.genes)

    p = ggplot(de.cpm.t[de.cpm.t$FC.is.pos,], aes_string(x=first.var, y="CPM",
            color=second.var, fill=second.var, pch=third.var)) +
        geom_point(alpha=0.5) +
        geom_jitter(width=0.2) +
        facet_wrap(. ~ gene, scales="free_y", ncol=4) +
        labs(title=plot.title, subtitle=paste0("+ve DE genes (adj. p < ", fdr.cutoff, ")")) +
        scale_color_brewer(palette="Dark2") +
        scale_fill_brewer(palette="Dark2") +
        scale_shape_manual(values=c(1, 16)) +
        theme_bw() +
        theme(strip.text = element_text(size = 6))

    ggsave(p, file=paste0("reports/stripplot_DE_genes-", contrast, ".pos.pdf"),
        width=4, height=7)

    p = ggplot(de.cpm.t[!de.cpm.t$FC.is.pos,], aes_string(x=first.var, y="CPM",
            color=second.var, fill=second.var, pch=third.var)) +
        geom_point(alpha=0.5) +
        geom_jitter(width=0.2) +
        facet_wrap(. ~ gene, scales="free_y", ncol=4) +
        labs(title=plot.title, subtitle=paste0("-ve DE genes (adj. p < ", fdr.cutoff, ")")) +
        scale_color_brewer(palette="Dark2") +
        scale_fill_brewer(palette="Dark2") +
        scale_shape_manual(values=c(1, 16)) +
        theme_bw() +
        theme(strip.text = element_text(size = 6))

    ggsave(p, file=paste0("reports/stripplot_DE_genes-", contrast, ".neg.pdf"),
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

do.expr_stripchart.plot("inf.TvF", fdr.cutoffs[["inf.TvF"]],
    plot.titles[["inf.TvF"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("inf.TvF.onlyD2", fdr.cutoffs[["inf.TvF.onlyD2"]],
    plot.titles[["inf.TvF.onlyD2"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("inf.TvF.onlyD7", fdr.cutoffs[["inf.TvF.onlyD7"]],
    plot.titles[["inf.TvF.onlyD7"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("inf.TvF.onlyF", fdr.cutoffs[["inf.TvF.onlyF"]],
    plot.titles[["inf.TvF.onlyF"]],
    "infection.status", "dpi", "sex")

do.expr_stripchart.plot("inf.TvF.onlyM", fdr.cutoffs[["inf.TvF.onlyM"]],
    plot.titles[["inf.TvF.onlyM"]],
    "infection.status", "dpi", "sex")

# --- DPI

do.expr_stripchart.plot("dpi.7v2", fdr.cutoffs[["dpi.7v2"]],
    plot.titles[["dpi.7v2"]],
    "infection.status", "sex", "dpi")

do.expr_stripchart.plot("dpi.7v2.onlyF", fdr.cutoffs[["dpi.7v2.onlyF"]],
    plot.titles[["dpi.7v2.onlyF"]],
    "infection.status", "dpi", "sex")

do.expr_stripchart.plot("dpi.7v2.onlyM", fdr.cutoffs[["dpi.7v2.onlyM"]],
    plot.titles[["dpi.7v2.onlyM"]],
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

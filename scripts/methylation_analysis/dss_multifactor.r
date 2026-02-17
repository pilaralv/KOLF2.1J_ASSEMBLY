args <- commandArgs(trailingOnly = TRUE)

# Access individual arguments
chromosome <- as.character(args[1])
celltype <- as.character(args[2])
outputDir <- as.character(args[3])
# pval <- as.numeric(args[4])      # for float / decimal numbers
# term <- as.integer(args[5])      # for integer index

# chromosome <- "chr22"
# celltype <- "neuron"
# outputDir <- "/data/Phillippy/projects/kolf_methylation/02.processedData/ont/dmr/"
# pval <- 0.05
# term <- 2

cat("chromosome :" , chromosome, "\n")
cat("celltype :", celltype, "\n")
cat("Output file:", outputDir, "\n")
# cat("pval file:", pval, "\n")
# cat("term file:", term, "\n")

# output 

# ips
ips_alt_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/IPS/IPS.modkit.pri.alt.filt.liftover_to_chm13.bed"
ips_pri_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/IPS/IPS.modkit.pri.pri.filt.liftover_to_chm13.bed"

# mgl
mgl_alt_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/MGL/MGL.modkit.pri.alt.filt.liftover_to_chm13.bed"
mgl_pri_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/MGL/MGL.modkit.pri.pri.filt.liftover_to_chm13.bed"

# neuron
neuron_alt_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/NEURON/NEURON.modkit.pri.alt.filt.liftover_to_chm13.bed"
neuron_pri_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/NEURON/NEURON.modkit.pri.pri.filt.liftover_to_chm13.bed"

# ngn2
ngn2_alt_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/NGN2/NGN2.modkit.pri.alt.filt.liftover_to_chm13.bed"
ngn2_pri_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/NGN2/NGN2.modkit.pri.pri.filt.liftover_to_chm13.bed"


# ASTRO_NYSCF
ASTRO_NYSCF_alt_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/ASTRO_NYSCF/ASTRO_NYSCF.modkit.pri.alt.filt.liftover_to_chm13.bed"
ASTRO_NYSCF_pri_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/ASTRO_NYSCF/ASTRO_NYSCF.modkit.pri.pri.filt.liftover_to_chm13.bed"

# ASTRO_IND12
ASTRO_IND12_alt_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/ASTRO_IND12/ASTRO_IND12.modkit.pri.alt.filt.liftover_to_chm13.bed"
ASTRO_IND12_pri_file="/vf/users/Phillippy/projects/kolf_methylation/02.processedData/ont/ASTRO_IND12/ASTRO_IND12.modkit.pri.pri.filt.liftover_to_chm13.bed"


astro_alt_file <- ASTRO_NYSCF_alt_file
astro_pri_file <- ASTRO_NYSCF_pri_file

suppressPackageStartupMessages(library(DSS))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# ips
ips_alt <- fread(ips_alt_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
ips_alt <- ips_alt[, .(chr, pos, N, X)]

ips_pri <- fread(ips_pri_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
ips_pri <- ips_pri[, .(chr, pos, N, X)]

# mgl
mgl_alt <- fread(mgl_alt_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
mgl_alt <- mgl_alt[, .(chr, pos, N, X)]

mgl_pri <- fread(mgl_pri_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
mgl_pri <- mgl_pri[, .(chr, pos, N, X)]

# neuron
neuron_alt <- fread(neuron_alt_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
neuron_alt <- neuron_alt[, .(chr, pos, N, X)]

neuron_pri <- fread(neuron_pri_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
neuron_pri <- neuron_pri[, .(chr, pos, N, X)]

# neuron
ngn2_alt <- fread(ngn2_alt_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
ngn2_alt <- ngn2_alt[, .(chr, pos, N, X)]

ngn2_pri <- fread(ngn2_pri_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
ngn2_pri <- ngn2_pri[, .(chr, pos, N, X)]

# astro
astro_alt <- fread(astro_alt_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
astro_alt <- astro_alt[, .(chr, pos, N, X)]

astro_pri <- fread(astro_pri_file, header = FALSE, col.names = c("chr", "pos", "end", "N", "X"))
astro_pri <- astro_pri[, .(chr, pos, N, X)]

# chromosome="chr1"

ips_pri <- ips_pri[ips_pri$chr == chromosome, ]
mgl_pri <- mgl_pri[mgl_pri$chr == chromosome, ]
neuron_pri <- neuron_pri[neuron_pri$chr == chromosome, ]
ngn2_pri <- ngn2_pri[ngn2_pri$chr == chromosome, ]
astro_pri <- astro_pri[astro_pri$chr == chromosome, ]

ips_alt <- ips_alt[ips_alt$chr == chromosome, ]
mgl_alt <- mgl_alt[mgl_alt$chr == chromosome, ]
neuron_alt <- neuron_alt[neuron_alt$chr == chromosome, ]
ngn2_alt <- ngn2_alt[ngn2_alt$chr == chromosome, ]
astro_alt <- astro_alt[astro_alt$chr == chromosome, ]

# ips_pri <- ips_pri[!ips_pri$chr %in% c("chrX", "chrY"), ]
# mgl_pri <- mgl_pri[!mgl_pri$chr %in% c("chrX", "chrY"), ]
# neuron_pri <- neuron_pri[!neuron_pri$chr %in% c("chrX", "chrY"), ]
# ngn2_pri <- ngn2_pri[!ngn2_pri$chr %in% c("chrX", "chrY"), ]

print("Loading DSS data...")

# j=33425941
# i=37055415
# BSobj = makeBSseqData( list(ips_pri[j:i,], ips_alt[j:i,], mgl_pri[j:i,], mgl_alt[j:i,], neuron_pri[j:i,], neuron_alt[j:i,], ngn2_pri[j:i,], ngn2_alt[j:i,]), c("ips_pri","ips_alt", "mgl_pri", "mgl_alt","neuron_pri","neuron_alt","ngn2_pri","ngn2_alt") )
BSobj = makeBSseqData( list(ips_pri, ips_alt, mgl_pri, mgl_alt, neuron_pri, neuron_alt, ngn2_pri, ngn2_alt, astro_pri, astro_alt),
c("ips_pri","ips_alt", "mgl_pri", "mgl_alt","neuron_pri","neuron_alt","ngn2_pri","ngn2_alt", "astro_pri", "astro_alt") )
BSobj

celltype_lst <- c("ips", "ips", "mgl", "mgl", "neuron", "neuron", "ngn2", "ngn2", "astro","astro")
celltype_lst <- ifelse(celltype_lst == celltype, celltype, paste0("not_", celltype))
celltype_lst = factor(celltype_lst, levels = c(paste0("not_", celltype),celltype) )
# celltype_lst <- factor(celltype_lst, levels = c("ips", "mgl", "neuron", "ngn2"))

allele = c("pri","alt","pri","alt","pri","alt","pri","alt","pri","alt")

design = data.frame(celltype_lst, allele)
design$celltype_lst <- relevel(design$celltype_lst, ref = paste0("not_", celltype))
design$allele <- relevel(factor(design$allele), ref = "pri")

rm(celltype_lst, allele)

# design

# model.matrix(~celltype_lst+allele+celltype_lst:allele, design)

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~celltype_lst+allele+celltype_lst:allele)
#DMLtest.cell = DMLtest.multiFactor(DMLfit, coef="celltypenot_ips:allelepri")

# DMLtest.cell_2 = DMLtest.multiFactor(DMLfit, coef=2) # coef=2 is celltype type
# DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=3) # coef=3 is allele
# DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=4) # coef=4 is interaction term

# par(bg="white")
# options(repr.plot.width=10, repr.plot.height=5, bg = "white")  # Change numbers as needed
# png(dml_stat, width=1000, height=500)
# hist(DMLtest.cell$stat, 100, main=paste0(celltype,"_",chromosome,"_",term,"_",pval,"_dml_stat"))
# dev.off()

# par(bg="white")
# options(repr.plot.width=10, repr.plot.height=5, bg = "white")  # Change numbers as needed
# png(dml_pval, width=1000, height=500)
# hist(DMLtest.cell$pvals, 1000, main=paste0(celltype,"_",chromosome,"_",term,"_",pval,"_dml_stat_Pvalues"))
# dev.off()

# dmrs2_2 =  callDMR(DMLtest.cell_2,p.threshold = pval, minCG = 10)
# head(dmrs2)

# write.csv(dmrs2_2, paste0(outputDir,"/",celltype,"_",chromosome,"_2_",pval,".csv"), row.names = FALSE, quote = FALSE)

for (pval in c(0.00001, 0.0005, 0.05)) {
for (term in c(2, 3, 4)) {
    print(paste("Term: ", term))
    DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=term)
    dmrs2 = callDMR(DMLtest.cell, p.threshold = pval)
    write.csv(dmrs2, paste0(outputDir, "/", celltype, "_", chromosome, "_", term, "_", pval, ".csv"), row.names = FALSE, quote = FALSE)
}
}

# if (is.null(dmrs2)) {
#   cat("dmrs2 data frame is NULL\n")
# } else {
#     cat("dmrs2 data frame is not NULL\n")
#     png(dmr_area, width=1000, height=500)
#    par(bg="white")
    # options(repr.plot.width=10, repr.plot.height=5, bg = "white")  # Change numbers as needed
    # hist(dmrs2$areaStat, 1000, main=paste0(celltype,"_",chromosome,"_",term,"_",pval,"_dmr_area"))
    # dev.off()

#     i=1
#    dmrs2[i,]

    # options(repr.plot.width=15, repr.plot.height=10, bg = "white")  # Change numbers as needed
    # png(top_dmr, width=1500, height=1000)
    # par(bg="white")
    # showOneDMR(dmrs2[i,], BSobj)
    # dev.off()
# }


sessionInfo()

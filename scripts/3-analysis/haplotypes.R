library(pegas)
library(optparse)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(stringr)
library(ComplexHeatmap)
library(dendsort)

POPCOLS = c(palette.colors(palette = "Paired"), "darkmagenta", "deeppink", "darkolivegreen","darkred", "navy", "darkslateblue", "gray18", "gray")

### COMMENT THIS SECTION TO DISABLE COMMAND LINE ARGUMENTS ###

option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "Input file (noheader.strct_in)"),
  make_option(c("-g", "--gds"), type = "character", help = "GDS file"),
  make_option(c("-p", "--population"), type = "character", help = "Table matching samples to population (.tsv)")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

###############################################################

### UNCOMMENT THIS SECTION TO RUN SCRIPT IN RStudio ###
# 
# opt = list()
# opt$input = "FUT1.vcf" # Input VCF, required
# opt$population = "population.tsv" # population table
# opt$gds = 'FUT1.gds' # Input gds
# 
#######################################################

NAME = str_split_i(basename(opt$input), fixed("."), 1)
DIR = str_interp("${getwd()}/${NAME}/haplotypes")

if (!dir.exists(DIR)) {
  dir.create(DIR, recursive = T)
}

vcf = pegas::read.vcf(opt$input)
genofile = seqOpen(opt$gds, allow.duplicate = T)
samples = seqGetData(genofile, "sample.id")
variants = seqGetData(genofile, "annotation/id")
variants[variants == ""] = seqGetData(genofile, "position")[variants == ""]
alleles = seqGetData(genofile, "allele")
variants = paste(variants, alleles, sep="_")

pop = read.table(file = opt$population, sep = "\t", row.names = 1)
pop = pop[,1]
if (sum(is.na(pop))>0) {
  warning(str_interp("No population assigned to the following samples : ${samples[is.na(pop)]}.\nAssigning them to a population named 'UNK'"))
  pop[is.na(pop)] = "UNK"
}
names(pop) = samples

CSQ = sapply(strsplit(seqGetData(genofile, "annotation/info/CSQ"),fixed = T, split = "|"), `[`, 2)
names(CSQ) = variants
CSQ.colors = palette.colors(palette = "Dark2", n = length(unique(CSQ)))
names(CSQ.colors) = unique(CSQ)
CSQ.vars = CSQ.colors[CSQ]
names(CSQ.vars) = variants

setwd(DIR)

h.uncomp = haplotype(vcf, locus=1:ncol(vcf), compress = F)
nbin = as.DNAbin(t(h.uncomp))
rownames(nbin) = rep(pop, each=2)

freq.pop = haploFreq(nbin, fac = rep(pop, each = 2))
freq.ind = haploFreq(nbin, fac = rep(samples, each = 2))

write.table(freq.pop, "haplotype_number_bypop.tsv", sep = "\t", col.names = NA, row.names = T)
write.table(freq.ind, "haplotype_number.tsv", sep = "\t", col.names = NA, row.names = T)

h.full = haplotype(nbin)
h = subset(h.full, minfreq=2)

a = as.character(h)
nh = nrow(h)
df = as.data.frame(a)
colnames(df) = variants
df = df[,colSums(is.na(df)) < nrow(df)]
df[is.na(df)] = "-"
counts = df
# freqs = sapply(as.data.frame(a), table)
for (pos in 1:ncol(counts)) {
  counts[,pos] = sapply(counts[,pos], function(x) {table(counts[,pos])[x]})
}
counts = counts[,colSums(counts) < nh^2]

ref = sapply(strsplit(colnames(counts), split="_"), `[`, 2)
ref = unlist(strsplit(ref, split=","))
malleles = matrix(ref, ncol = 2, byrow = T)
malleles[nchar(malleles)>1] = '-'
malleles = tolower(malleles)

c = CSQ[colnames(counts)]
csq = unique(c)
csq = csq[order(csq)]
csq.colors = palette.colors(palette = "Dark2", n = length(csq))
names(csq.colors) = csq
csq.vars = csq.colors[c]

df = df[,colnames(counts)]

hf = freq.ind[rowSums(freq.ind)>1,]
hf = freq.ind[rownames(counts),]
write.csv(hf, "majhapind.csv")

dend = hclust(dist.hamming(h))
id.count = rowSums(hf)
id.count = sapply(names(id.count), function(x) {str_interp("${x} ($[.1f]{id.count[x]/(length(pop)*2)*100}%)")})

for (j in 1:ncol(df)) {
  for (i in 1:nrow(df)) {
    if (df[i,j] == malleles[j,1]) {
      df[i,j] = "REF"
    } else if (df[i,j] == malleles[j,2]){
      df[i,j] = "ALT"
    } else {
      df[i,j] = NA
    }
  }
}
rownames(df) = id.count

ha = HeatmapAnnotation(
  consequence = names(csq.vars),
  col = list(
    consequence = csq.colors
  ),
  annotation_legend_param = list(
    consequence = list(
      title = "Variant consequence",
      at = names(csq.colors)
    )
  )
)

row_dend = dendsort(hclust(dist.hamming(df)))
#lgd = Legend(labels = names(csq.colors), legend_gp = gpar(fill = csq.cols), ncol = 1)
#dev.off()
#draw(hb, y = unit(0.8, "npc"))
#draw(lgd, y=unit(0.2, "npc"))

png("allele_frequency_majority.png", height = 1500, width = 1500)
Heatmap(df,
        col = c("#414770","#171123"),
        #clustering_distance_rows = function(x) {dist.hamming(x)},
        cluster_rows = row_dend,
        row_dend_reorder = T,
        row_dend_width = unit(25, "mm"),
        bottom_annotation = ha,
        name = "Allele",
        show_column_names = F,
        gap = unit(0.5, "mm"),
        split = nrow(df),
        row_title = NULL
)
dev.off()

png("haplotype_distribution_majority.png", height = 1500, width = 1500)
Heatmap(hf,
        col = c("gray","coral2","sienna4"),
        cluster_rows = row_dend,
        cluster_columns = T,
        show_column_names = F,
        name = "Zygocity",
        heatmap_legend_param = list(
          at = c(1,2),
          labels = c("Heterozygous", "Homozygous")
        ),
        row_split = nrow(hf),
        gap = unit(0.5, "mm"),
        row_title = NULL,
        column_split = pop)
dev.off()

b = as.character(h.full)
nh = nrow(h.full)
df = as.data.frame(b)
colnames(df) = variants
df = df[,colSums(is.na(df)) < nrow(df)]
df[is.na(df)] = "-"
counts = df
# freqs = sapply(as.data.frame(a), table)
for (pos in 1:ncol(counts)) {
  counts[,pos] = sapply(counts[,pos], function(x) {table(counts[,pos])[x]})
}
counts = counts[,colSums(counts) < nh^2]

df = df[,colnames(counts)]

#hf = freq.ind[rowSums(freq.ind)>1,]
hf = freq.ind[rownames(counts),]

ref = sapply(strsplit(colnames(counts), split="_"), `[`, 2)
ref = unlist(strsplit(ref, split=","))
malleles = matrix(ref, ncol = 2, byrow = T)
malleles[nchar(malleles)>1] = '-'
malleles = tolower(malleles)

dend = hclust(dist.hamming(h.full))
id.count = rowSums(hf)
id.count = sapply(names(id.count), function(x) {str_interp("${x} ($[.1f]{id.count[x]/(length(pop)*2)*100}%)")})

for (j in 1:ncol(df)) {
  for (i in 1:nrow(df)) {
    if (df[i,j] == malleles[j,1]) {
      df[i,j] = "REF"
    } else if (df[i,j] == malleles[j,2]){
      df[i,j] = "ALT"
    } else {
      df[i,j] = NA
    }
  }
}
rownames(df) = id.count

c = CSQ[colnames(counts)]
csq = unique(c)
csq = csq[order(csq)]
csq.colors = palette.colors(palette = "Dark2", n = length(csq))
names(csq.colors) = csq
csq.vars = csq.colors[c]

ha = HeatmapAnnotation(
  consequence = names(csq.vars),
  col = list(
    consequence = csq.colors
  ),
  annotation_legend_param = list(
    consequence = list(
      title = "Variant consequence",
      at = names(csq.colors)
    )
  )
)

row_dend = dendsort(hclust(dist.hamming(df)))
#lgd = Legend(labels = names(csq.colors), legend_gp = gpar(fill = csq.cols), ncol = 1)
#dev.off()
#draw(hb, y = unit(0.8, "npc"))
#draw(lgd, y=unit(0.2, "npc"))

png("allele_frequency_all.png", height = 1500, width = 1500)
Heatmap(df,
        col = c("#414770","#171123"),
        #clustering_distance_rows = function(x) {dist.hamming(x)},
        cluster_rows = row_dend,
        row_dend_reorder = T,
        row_dend_width = unit(25, "mm"),
        border="black",
        bottom_annotation = ha,
        name = "Allele",
        show_column_names = F,
        row_title = NULL,
        gap = unit(1, "mm")
)
dev.off()

png("haplotype_distribution_all.png", height = 1500, width = 1500)
Heatmap(hf,
        col = c("gray","coral2","sienna4"),
        cluster_rows = row_dend,
        cluster_columns = T,
        show_column_names = F,
        name = "Zygocity",
        heatmap_legend_param = list(
          at = c(1,2),
          labels = c("Heterozygous", "Homozygous")
        ),
        column_split = pop)
dev.off()

colors = POPCOLS[1:length(unique(pop))]
names(colors) = unique(pop)

h.dist = dist.dna(h, "N")

network = rmst(h.dist)

size = summary(h)
size = size[labels(network)]

freq.pop = haploFreq(nbin, fac = rep(rownames(vcf), each = 2))
freq.pop = freq.pop[labels(network),]

filter = freq.pop[rowSums(freq.pop) > 1,]
filter = filter[,colSums(filter) > 0]
colnames(filter) = lapply(strsplit(colnames(filter), split="-"), `[`, 1)
filter = t(rowsum(t(filter), group = colnames(filter)))

png("haplotype_network.png", height = 1500, width = 1500)
plot(network, size=sqrt(size), show.mutation=1, labels = F, pie = filter, legend = c(-55,0), bg = colors)
dev.off()

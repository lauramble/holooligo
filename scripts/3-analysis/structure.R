library(LEA)
library(stringr)
library(RColorBrewer)
library(optparse)
library(ComplexHeatmap)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")

POPCOLS = c(palette.colors(palette = "Paired"), "darkmagenta", "deeppink", "darkolivegreen","darkred", "navy", "darkslateblue", "gray18", "gray")

### COMMENT THIS SECTION TO DISABLE COMMAND LINE ARGUMENTS ###
option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "Input file (noheader.strct_in)"),
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
# opt$corr_threshold = '-' # correlation threshold for SNP selection
# 
#######################################################

NAME = str_split_i(basename(opt$input), fixed('.'), 1)
DIR = str_interp("${getwd()}/${NAME}/structure")

if (!dir.exists(DIR)) {
  dir.create(DIR, recursive = T)
}

pop = read.table(file = opt$population, sep = "\t", row.names = 1)
samples = rownames(pop)
pop = pop[,1]
if (sum(is.na(pop))>0) {
  warning(str_interp("No population assigned to the following samples : ${samples[is.na(pop)]}.\nAssigning them to a population named 'UNK'"))
  pop[is.na(pop)] = "UNK"
}
names(pop) = samples

colors = POPCOLS[1:length(unique(pop))]
names(colors) = unique(pop)

Ks = 2:(length(unique(pop)))

setwd(DIR)

struct2geno(opt$input, TESS=F, diploid=T, FORMAT = 1, extra.col = 2, extra.row = 0, output = str_interp("table.geno"))
obj.snmf = snmf("table.geno", K = 2:15, ploidy = 2, alpha = 10, entropy = T, CPU = 4, project = "new", repetitions = 3, seed = 1000) 

png(str_interp("cross_entropy.png"), height = 750, width = 750)
plot.snmf = plot(obj.snmf, cex = 1.4, pch = 19, main = str_interp("Cross-entropy plot for ${NAME}"))
print(plot.snmf)
dev.off()


for (K in Ks) {
  for (run in 1:3) {
    qmatrix = Q(obj.snmf, K = K, run = run)
    print(pop)
    print(dim(qmatrix))
    rownames(qmatrix) = pop
    qmatrix = as.data.frame(qmatrix)
    #qmatrix = qmatrix[order(-qmatrix[,1], -qmatrix[,2]),]
    #qmatrix = qmatrix[do.call(order, c(decreasing=T, data.frame(qmatrix[,1:K]))),]
    
    qmat.sortq = cbind(qmatrix, max = apply(qmatrix, 1, max), match = apply(qmatrix, 1, function(x) {colnames(qmatrix)[x == max(x)]}))
    qmat.sortq = qmat.sortq[order(qmat.sortq[,"match"], -qmat.sortq[,"max"]),]
    
    x = rownames(qmat.sortq)
    
    hb = HeatmapAnnotation(
      structure = anno_barplot(qmat.sortq[,colnames(qmatrix)], 
                               gp = gpar(fill = rainbow(K)), 
                               bar_width = 1, 
                               height = unit(12, "cm"), 
                               width = unit(24, "cm")
      ), 
      breed = anno_simple(sapply(strsplit(x, split=".", fixed = T), `[`, 1),
                          width = unit(24, "cm"),
                          col = colors)
    )
    lgd = Legend(labels = unique(pop), legend_gp = gpar(fill = colors), nrow = 1)
    png(str_interp("K${K}_run${run}.png"), height = 1000, width = 1500)
    draw(hb, y = unit(0.8, "npc"))
    draw(lgd, y=unit(0.2, "npc"))
    dev.off()
  }
}


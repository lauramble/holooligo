library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(stringr)
library(optparse)
library(ComplexHeatmap)
library(circlize)


### COMMENT THIS SECTION TO DISABLE COMMAND LINE ARGUMENTS ###
if (!exists(opt)){
  option_list = list(
    make_option(c("-i", "--input"), type = "character", help = "Input file (.vcf or .gds)"),
    make_option(c("-p", "--population"), type = "character", help = "Table matching samples to population (.tsv)"),
    make_option(c("-c", "--inferchar"), type = "character", help = "Character used to split sample names for population inference (i.e. '-' if the sample names follow the format 'POP1-SAMPLE1', 'POP1-SAMPLE2', 'POP2-SAMPLE1', etc..)"),
    make_option(c("-P", "--inferpos"), type = "character", help = "Position of the population tag in sample names for population inference (i.e. 1 for 'POP-SAMPLEID', 2 for 'SAMPLEID-POP', etc...)", default = 1)
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
}
###############################################################

### UNCOMMENT THIS SECTION TO RUN SCRIPT IN RStudio ###
# 
# opt = list()
# opt$input = "FUT1.vcf" # Input VCF, required
# opt$population = "population.tsv" # population table
# opt$inferchar = '-' # character used to split sample names to infer population
# opt$inferpos = 1 # position of the population tag in sample names
# 
#######################################################

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file required")
} else {
  format = tolower(str_split_i(basename(opt$input), fixed('.'), -1))
  if (format == "vcf") {
    gds = sub(".vcf", ".gds", basename(opt$input), fixed = T)
    seqVCF2GDS(opt$input, gds)
  } else if (format == "gds") {
    gds = opt$gds
  } else {
    stop("Wrong input file format. Must be gds or vcf.")
  }
}

name = str_split_i(basename(opt$input), fixed('.'), 1)
outdir = str_interp("${getwd()}/${name}")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

genofile = seqOpen(gds, allow.duplicate = T)
samples = seqGetData(genofile, "sample.id")
variants = seqGetData(genofile, "annotation/id")
variants[variants == ""] = seqGetData(genofile, "position")[variants == ""]
alleles = seqGetData(genofile, "allele")

varheader = seqSummary(gdsfile = gds)$info$Description
varheader = varheader[grep("VEP", varheader)]
varheader = str_match(varheader, "[a-zA-Z0-9]*\\|.*[a-zA-Z0-9]")[1,1]
varheader = unlist(str_split(varheader, "\\|"))
nfields = length(varheader)

csqtab = str_split_fixed(seqGetData(genofile, "annotation/info/CSQ"), "\\|", nfields)
colnames(csqtab) = varheader
csqtab = cbind(ID = variants, chr = seqGetData(genofile, "chromosome"), position = seqGetData(genofile, "position"), REF = str_split_i(alleles, ",", 1), csqtab)
write.table(csqtab, str_interp("${outdir}/variant_table.tsv"), sep = "\t", row.names = F)

t = applyMethod(genofile, alleleFrequency)
names(t) = variants
t = t(data.frame(t))
rownames(t) = "RefFreq"
write.table(t, str_interp("${outdir}/reference_allele_freq.tsv"), sep = '\t', dec = ",", col.names = NA)

if (!is.null(opt$population)) {
  if (!is.null(opt$inferchar)) {
    warning("A population table has been provided, therefore population inference will not be performed")
  }
  pop = read.table(file = opt$population, sep = "\t", row.names = 1)
  pop = pop[samples,1]
  if (sum(is.na(pop))>0) {
    warning(str_interp("No population assigned to the following samples : ${samples[is.na(pop)]}.\nAssigning them to a population named 'UNK'"))
    pop[is.na(pop)] = "UNK"
  }
  names(pop) = samples
} else if (!is.null(opt$inferchar)) {
  message(str_interp("Using split character ${opt$inferchar} and population position ${opt$inferpos} to extract population data from sample names"))
  pop = str_split_i(samples, fixed(opt$inferchar), opt$inferpos)
} else {
  stop("No population table or split character for population inference. Skipping population statistics.")
}

### Population overview

png(str_interp("${outdir}/population-overview.png"))
pie(table(pop), main = str_interp("Population overview in ${basename(opt$input)}"))
dev.off()

pop.breakdown = data.frame(table(pop))
colnames(pop.breakdown) = c("population", "number")
write.table(pop.breakdown, str_interp("${outdir}/population_breakdown.tsv"), sep = '\t', quote = F, row.names = F, col.names = T )

### Allele frequency by population

npop = unique(pop)

af = matrix(nrow = length(npop), ncol = length(variants))
rownames(af) = npop
colnames(af) = paste(variants, alleles, sep="_")
for (p in npop) {
  samples = seqGetData(genofile, "sample.id")[pop == p]
  if (length(samples)>0){
    t = applyMethod(genofile, alleleFrequency, sample=samples)
    af[p,] = t
  }
}
write.table(af, str_interp("${outdir}/reference_allele_freq_pop.tsv"), sep = "\t", dec = ",", col.names = NA)

CSQ = sapply(strsplit(seqGetData(genofile, "annotation/info/CSQ"),fixed = T, split = "|"), `[`, 2)
names(CSQ) = colnames(af)
CSQ.colors = palette.colors(palette = "Dark2", n = length(unique(CSQ)))
names(CSQ.colors) = unique(CSQ)
CSQ.vars = CSQ.colors[CSQ]
names(CSQ.vars) = colnames(af)

ha = HeatmapAnnotation(
  consequence = CSQ,
  col = list(
    consequence = CSQ.colors
  ),
  annotation_legend_param = list(
    consequence = list(
      title = "Variant consequence",
      at = names(CSQ.colors)
    )
  )
)

png(filename = str_interp("${outdir}/reference_allele_freq_pop.png"), width = 1500, height = 1000)
Heatmap(af, cluster_columns = F, col=colorRamp2(c(0,1), colors=c("#BDB4CF","#39304B")), show_column_names = F,
        name = "Reference allele frequency by population", bottom_annotation = ha)
dev.off()

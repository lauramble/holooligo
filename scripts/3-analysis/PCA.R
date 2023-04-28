library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(optparse)
library(plotly)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(htmlwidgets)

POPCOLS = c(palette.colors(palette = "Paired"), "darkmagenta", "deeppink", "darkolivegreen","darkred", "navy", "darkslateblue", "gray18", "gray")

### COMMENT THIS SECTION TO DISABLE COMMAND LINE ARGUMENTS ###
if (!exists(opt)){
  option_list = list(
    make_option(c("-i", "--input"), type = "character", help = "Input file (.vcf or .gds)"),
    make_option(c("-p", "--population"), type = "character", help = "Table matching samples to population (.tsv)"),
    make_option(c("-c", "--corr_threshold"), type = "double", help = "Absolute correlation threshold for SNP selection", default = 0.8)
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
# opt$corr_threshold = '-' # correlation threshold for SNP selection
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

draw_pca = function(pca, samples, name, pop = NULL, outdir = '.', colors = NULL, npc = 4, type = "ggplot"){
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = T)
  }
  
  # "plotly3d", "plotly", "ggplot"
  
  pcpct = pca$varprop*100
  
  pcs = paste0("PC", 1:npc)
  pcnames = paste0(pcs, " (", format(pcpct[1:npc], digits=2), "%)", sep="")
  names(pcnames) = pcs
  eigenvect = pca$eigenvect[,1:npc]
  colnames(eigenvect) = pcs
  
  if (!is.null(pop)) {
    if (is.null(colors)) {
      colors = rainbow(length(unique(pop)))
    }
    colors = colors[1:length(unique(pop))]
    names(colors) = unique(pop)
  }
  
  if (type != "plotly3d") {
    pccomb = combn(pcs, 2)
    
    if (type  == "ggplot") {
      for (col in 1:ncol(pccomb)) {
        pc1 = pccomb[1,col]
        pc2 = pccomb[2,col]
        
        if (is.null(pop)) {
          temp = data.frame(id = samples, 
                       PC1 = eigenvect[,pc1], 
                       PC2 = eigenvect[,pc2])
          
          png(str_interp("${outdir}/${pc1}_${pc2}.png"), height = 750, width = 750)
          p = ggplot(data = temp, aes(x = PC1, y = PC2)) +
            geom_point() +
            labs(x = pcnames[pc1], y =pcnames[pc2],
                 title = str_interp("PCA for ${name}"),
                 subtitle = str_interp("Principal components ${pc1} and ${pc2}"))
          print(p)
          dev.off()
          
        } else {
          temp = data.frame(id = samples, 
                       pop = pop,
                       PC1 = eigenvect[,pc1], 
                       PC2 = eigenvect[,pc2])
          
          png(str_interp("${outdir}/${pc1}_${pc2}.png"), height = 750, width = 750)
          p = ggplot(data = temp, aes(x = PC1, y = PC2, col = pop)) +
            geom_point() +
            labs(x = pc1, y = pc2,
                 title = str_interp("PCA for ${name}"),
                 subtitle = str_interp("Principal components ${pc1} and ${pc2}")) +
            scale_color_manual(values = colors)
          print(p)
          dev.off()
        }
      }
    } else if (type == "plotly") {
      if (is.null(pop)) {
        for (col in 1:ncol(pccomb)) {
          pc1 = pccomb[1,col]
          pc2 = pccomb[2,col]
          
          p = plot_ly()
          p = p %>% add_markers(
            x = eigenvect[,pc1],
            y = eigenvect[,pc2],
            text = samples
          ) %>% layout(
            title = str_interp("PCA for ${name}"),
            xaxis = list(title = pcnames[pc1]),
            yaxis = list(title = pcnames[pc2])
          )
          saveWidget(p, str_interp("${outdir}/${pc1}_${pc2}.html"), selfcontained = F)
        }
      } else {
        for (col in 1:ncol(pccomb)) {
          pc1 = pccomb[1,col]
          pc2 = pccomb[2,col]
          
          p = plot_ly()
          p = p %>% add_markers(
            x = eigenvect[,pc1],
            y = eigenvect[,pc2],
            text = samples,
            color = pop,
            colors = colors
          ) %>% layout(
            title = str_interp("PCA for ${name}"),
            xaxis = list(title = pcnames[pc1]),
            yaxis = list(title = pcnames[pc2])
          )
          saveWidget(p, str_interp("${outdir}/${pc1}_${pc2}.html"), selfcontained = F)
        }
      }
    }
  } else  {
      pccomb = combn(pcs, 3)
      
      if (is.null(pop)) {
        for (col in 1:ncol(pccomb)) {
          pc1 = pccomb[1,col]
          pc2 = pccomb[2,col]
          pc3 = pccomb[3,col]
          
          p = plot_ly()
          p = p %>% add_markers(
            x = eigenvect[,pc1],
            y = eigenvect[,pc2],
            z = eigenvect[,pc3],
            text = samples
          ) %>% layout(
            title = str_interp("PCA for ${name}"),
            scene = list(
              xaxis = list(title = pcnames[pc1]),
              yaxis = list(title = pcnames[pc2]),
              zaxis = list(title = pcnames[pc3])
            )
          )
          saveWidget(p, str_interp("${outdir}/${pc1}_${pc2}_${pc3}.html"), selfcontained = F)
        }
      } else {
        for (col in 1:ncol(pccomb)) {
          pc1 = pccomb[1,col]
          pc2 = pccomb[2,col]
          pc3 = pccomb[3,col]
          
          p = plot_ly()
          p = p %>% add_markers(
            x = eigenvect[,pc1],
            y = eigenvect[,pc2],
            z = eigenvect[,pc3],
            text = samples,
            colors = colors,
            color = pop
          ) %>% layout(
            title = str_interp("PCA for ${name}"),
            scene = list(
              xaxis = list(title = pcnames[pc1]),
              yaxis = list(title = pcnames[pc2]),
              zaxis = list(title = pcnames[pc3])
            )
          )
          saveWidget(p, str_interp("${outdir}/${pc1}_${pc2}_${pc3}.html"), selfcontained = F)
        }
      }
    }
}

genofile = seqOpen(gds, allow.duplicate = T)
samples = seqGetData(genofile, "sample.id")
variants = seqGetData(genofile, "annotation/id")
variants[variants == ""] = seqGetData(genofile, "position")[variants == ""]
alleles = seqGetData(genofile, "allele")

pca = snpgdsPCA(genofile, num.thread=2, autosome.only = F, bayesian = T)

name = str_split_i(basename(opt$input), fixed('.'), 1)
outdir = str_interp("${getwd()}/${name}/PCA")

draw_pca(pca, samples, name, type = "ggplot", outdir = str_interp("${outdir}/ggplot"))
draw_pca(pca, samples, name, type = "plotly", outdir = str_interp("${outdir}/plotly"))
draw_pca(pca, samples, name, type = "plotly3d", outdir = str_interp("${outdir}/plotly3d"))

pcs = paste0("PC", 1:4)
pcpct = pca$varprop*100
pcnames = paste0(pcs, " (", format(pcpct[1:4], digits=2), "%)", sep="")

png(str_interp("${outdir}/pca_pairs.png"), height = 1000, width = 1000)
pairs(pca$eigenvect[,1:4], labels=pcnames, lower.panel = NULL, pch = 17, cex = 2)
dev.off()

chr = seqGetData(genofile, "chromosome")
corr = snpgdsPCACorr(pca, genofile, eig.which=1:4)

savepar = par(mfrow=c(2,2), mai=c(0.45, 0.55, 0.1, 0.25))
for (i in 1:4)
{
  png(str_interp("${outdir}/PC${i}_snp_correlation.png"), height = 1000, width = 1000)
  plot(abs(corr$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste0("PC", i),
       col=palette.colors(palette = "Dark2", n = length(unique(chr))), pch="+")
  abline(h = opt$corr_threshold)
  dev.off()
}

sel.snp = corr$snpcorr
rownames(sel.snp) = pcnames
colnames(sel.snp) = variants
sel.snp = sel.snp[,colSums(abs(sel.snp) > sel.snp) > 0]
write.table(sel.snp, str_interp("${outdir}/selected_snp.tsv"), sep = "\t", col.names = NA, row.names = T)

if (!is.null(opt$population)) {
  outdir = str_interp("${outdir}/bypop")
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = T)
  }
  
  pop = read.table(file = opt$population, sep = "\t", row.names = 1)
  pop = pop[samples,1]
  if (sum(is.na(pop))>0) {
    warning(str_interp("No population assigned to the following samples : ${samples[is.na(pop)]}.\nAssigning them to a population named 'UNK'"))
    pop[is.na(pop)] = "UNK"
  }
  names(pop) = samples
  
  colors = POPCOLS[1:length(unique(pop))]
  names(colors) = unique(pop)
  
  png(str_interp("${outdir}/pca_pairs.png"), height = 1000, width = 1000)
  pairs(pca$eigenvect[,1:4], col=colors[pop], labels=pcnames, lower.panel = NULL, pch = 17, cex = 2)
  dev.off()
  
  draw_pca(pca, samples, name, type = "ggplot", colors = POPCOLS, pop = pop, outdir = str_interp("${outdir}/ggplot"))
  draw_pca(pca, samples, name, type = "plotly", colors = POPCOLS, pop = pop, outdir = str_interp("${outdir}/plotly"))
  draw_pca(pca, samples, name, type = "plotly3d", colors = POPCOLS, pop = pop, outdir = str_interp("${outdir}/plotly3d"))
}

#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(data.table)
library(lemon)
library(bedr)
library(optparse)

# Read in inputs. They are:
# 1. [1] Markers file output from BOLT
# 2. [2] Mask name - MUST be the same as in the per-gene marker file
# 3. [3] MAF name - MUST be the same as in the per-gene marker file
# 4. [4] ENST ID of the gene you wish to plot
option_list = list(make_option(c("-m","--markers"), action="store", default=NA, type="character", help="Path to the markers.tsv.gz output from mrcepid-runassociationtesting."),
                   make_option(c("-a","--mask"), action="store", default=NA, type="character", help="Name of the MASK that you want to plot. See README.md for possible MASK values."),
                   make_option(c("-f","--maf"), action="store", default=NA, type="character", help="Name of the MAF that you want to plot. See README.md for possible MAF values."),
                   make_option(c("-e","--enst"), action="store", default=NA, type="character", help="ENST ID of the gene you want to plot. Either ENST or SYMBOL is required."),
                   make_option(c("-s","--symbol"), action="store", default=NA, type="character", help="SYMBOL of the gene you want to plot. Either ENST or SYMBOL is required."),
                   make_option(c("-p","--phenotype"), action="store", default=NA, type="character", help="Name of the phenotype. This is only used to add a phenotype name to the resulting plot."),
                   make_option(c("-o","--output"), action="store", default=NA, type="character", help="output prefix. Will name files like <output>.png & <output>.tsv"),
                   make_option(c("-v","--pvalue"), action="store", default=5.0e-8, type="double", help="P. value threshold to label Allele Counts [5.0e-8]"),
                   make_option(c("-t","--transcripts"), action="store", default='data/transcripts.tsv.gz', type="character", help="Path to the transcripts annotation file [data/transcripts.tsv.gz]"),
                   make_option(c("-x","--exons"), action="store", default='data/exon_models.txt.gz', type="character", help="Path to exon models [data/exon_models.txt.gz]"),
                   make_option(c("-y","--ymax"), action="store", default=10, type="double", help="Max -log10(p-value) (y-axis) when plotting. Values less than 1 will be -log10 transformed [10]."))
opt = parse_args(OptionParser(option_list=option_list))

if (!is.na(opt$markers)) {
  markers_file <- opt$markers
} else {
  stop("Mask name (-m, --markers) not provided. Please try again!")
}
if (!is.na(opt$mask)) {
  mask <- opt$mask
} else {
  stop("Mask name (-a, --mask) not provided. Please try again!")
}
if (!is.na(opt$maf)) {
  maf <- opt$maf
} else {
  stop("MAF name (-f, --MAF) not provided. Please try again!")
}

# Load transcripts and get info for gene of interest
if (file.exists(opt$transcripts)) {
  transcripts <- fread(opt$transcripts)
} else {
  stop("File provided to --transcripts does not exist. The default path (data/transcripts.tsv.gz) only works from the root directory of this repository.")
}
enst <- opt$enst
symbol <- opt$symbol
if (!is.na(enst)) {
  gene_info <- transcripts[ENST == enst]
} else if (!is.na(symbol)) {
  gene_info <- transcripts[SYMBOL == symbol]
  if (nrow(gene_info) > 1) {
    stop(paste0("Found more than one gene with symbol ", symbol, ". Please try again, possibly with ENST."))
  }
} else {
  stop("Neither gene ENST (-e/--enst) nor Symbol (-s/--symbol) were provided. Please try again!")
}

phenotype <- opt$phenotype
if (!is.na(opt$output)) {
  file_prefix <- opt$output
} else {
  file_prefix <- paste(paste(gene_info[,SYMBOL], mask, maf, "BOLT",sep = "_"))
}
pvalue <- opt$pvalue

if (opt$ymax < 1) {
  cat("-y / --ymax is less than 1, scaling to -log10(ymax)...")
  plot_ymax <- -log10(opt$ymax)
} else {
  plot_ymax <- opt$ymax
}
# And set y-axis limits/labels/breaks according to plot_ymax
plot_ylimits <- c((-1 * plot_ymax) - 3,plot_ymax + 3)
plot_ylabels <- seq(0, plot_ymax, by = plot_ymax / 5)
plot_ybreaks <- plot_ylabels + 1.5
plot_ybreaks <- c(-1 * rev(plot_ybreaks), plot_ybreaks)
plot_ylabels <- c(rev(plot_ylabels), plot_ylabels)

# Eugene's Default Theme:
theme <- theme(panel.background=element_rect(fill="white"),
               line=element_line(size=1,colour="black",lineend="round"),
               axis.line=element_line(size=1),
               text=element_text(size=16,face="bold",colour="black"),
               axis.text=element_text(colour="black"),
               axis.ticks=element_line(size=1,colour="black"),
               axis.ticks.length=unit(.1,"cm"),
               axis.ticks.x = element_blank(),
               strip.background=element_rect(fill="white"),
               axis.text.x = element_blank(),
               legend.position="right",
               panel.grid.major=element_blank(),
               legend.key=element_blank())


# Now load the gene's exon model:
if (file.exists(opt$exons)) {
  exon_models <- fread(opt$exons)
} else {
  stop("File provided to --exons does not exist. The default path (data/exon_models.txt.gz) only works from the root directory of this repository.")
}
setnames(exon_models, names(exon_models), c("transcript","cdsStart","cdsEnd","numExons","exonStarts","exonEnds","transcriptType"))
exon_models[,transcript:=str_split(transcript,"\\.",simplify = T)[1],by=1:nrow(exon_models)]
ensembl_annotation_gene <- exon_models[transcript == gene_info[,ENST]]

# Decide which query to use:
query <- switch (mask, 
                 "HC_PTV"='FILTER=="PASS" & PARSED_CSQ=="PTV" & LOFTEE=="HC"',
                 "PTV"='FILTER=="PASS" & PARSED_CSQ=="PTV"',
                 "MISS"='FILTER=="PASS" & PARSED_CSQ=="MISSENSE"',
                 "MISS_CADD25"='FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & CADD>=25',
                 "MISS_REVEL0_5"='PARSED_CSQ=="MISSENSE" & REVEL>=0.5',
                 "MISS_REVEL0_7"='FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & REVEL>=0.7',
                 "MISS_CADD25_REVEL0_7"='PARSED_CSQ=="MISSENSE" & CADD>=25 & REVEL>=0.7',
                 "SYN"='FILTER=="PASS" & PARSED_CSQ=="SYN"',
                 "DMG"='FILTER=="PASS" & ((PARSED_CSQ=="PTV" & LOFTEE=="HC") | (PARSED_CSQ=="MISSENSE" & CADD>=25))'
        )

if (is.null(query)) {
  stop(paste0("Mask name - ", mask, " - was not parsed properly. Please try again!"))
}


if (maf == "MAF_01") {
  query <- paste0(query, ' & AF<0.001')
} else if (maf == "AC_1") {
  query <- paste0(query, ' & AC==1')
} else {
  stop(paste0("MAF name - ", maf, " - was not parsed properly. Please try again!"))
}

# Grab Current Tabix Header and extract variants with tabix
res <- system(paste0("zcat < ", markers_file, " | head -n 1"), wait = T, intern = T)
if (length(res) == 0) {
  res <- system(paste0("zcat ", markers_file, " | head -n 1"), wait = T, intern = T)
  if (length(res) == 0) {
    stop(paste0("could not parse the header of the provided markers file..."))
  }
}
tabix.header <- str_split(res, "\t")[[1]]
variants <- data.table(tabix(region = gene_info[,coord], file.name = markers_file))
setnames(variants, names(variants), tabix.header)
# This is dumb, because the tabix parser doesn't check if things are integers or not...
variants[,AF:=as.double(AF)]
variants[,CADD:=as.double(CADD)]
variants[,REVEL:=as.double(REVEL)]
variants[,BOLT_AC:=as.double(BOLT_AC)]
variants[,POS:=as.integer(POS)]
variants[,P_BOLT_LMM_INF:=as.double(P_BOLT_LMM_INF)]
variants[,BETA:=as.double(BETA)]

# Make sure the right gene is included
variants <- variants[ENST == gene_info[,ENST]]

# Make sure the right variants are include
variants <- variants[eval(parse(text=query))]
variants <- variants[BOLT_AC > 0]

# Here we plot the actual lolliplot
## First step is to setup the appropriate gene model with "tiny" introns so that the reader can focus on coding sequence
# Create initial exon breakpoints
coding.start <- ensembl_annotation_gene[,cdsStart]
coding.end <- ensembl_annotation_gene[,cdsEnd]

exon.starts <- as.integer(unlist(str_split(ensembl_annotation_gene[,exonStarts],",")))
exon.starts <- exon.starts[1:length(exon.starts)-1]
exon.ends <-  as.integer(unlist(str_split(ensembl_annotation_gene[,exonEnds],",")))
exon.ends <- exon.ends[1:length(exon.ends)-1]

pos.map <- data.table()
gene.map <- data.table()
last.pos <- 1
for (i in 1:length(exon.starts)) {
  if (exon.starts[i] < coding.start && exon.ends[i] < coding.start) {
    # all UTR
    gene.map <- rbind(gene.map,
                      data.table(start = exon.starts[i],
                                 end = exon.ends[i],
                                 ymin = -0.25,
                                 ymax = 0.25,
                                 annotation = "utr"
                      ))
  } else if (exon.starts[i] > coding.end && exon.ends[i] > coding.end) {
    gene.map <- rbind(gene.map,
                      data.table(start = exon.starts[i],
                                 end = exon.ends[i],
                                 ymin = -0.25,
                                 ymax = 0.25,
                                 annotation = "utr"
                      ))
  } else if (exon.starts[i] < coding.start && exon.ends[i] > coding.end) {
    # Single exon gene <-  make three exons:
    gene.map <- rbind(gene.map,
                      data.table(start = c(exon.starts[i],coding.start,coding.end),
                                 end = c(coding.start,coding.end,exon.ends[i]),
                                 ymin = c(-0.25,-0.5,-0.25),
                                 ymax = c(0.25,0.5,0.25),
                                 annotation = c("utr","cds","utr")
                      ))
  } else if (exon.starts[i] < coding.start && exon.ends[i] > coding.start) {
    # CDS + UTR <- make two exons
    gene.map <- rbind(gene.map,
                      data.table(start = c(exon.starts[i],coding.start),
                                 end = c(coding.start,exon.ends[i]),
                                 ymin = c(-0.25,-0.5),
                                 ymax = c(0.25,0.5),
                                 annotation = c("utr","cds")
                      ))
  } else if (exon.starts[i] < coding.end && exon.ends[i] > coding.end) {
    # CDS + UTR <- make two exons
    # CDS + UTR <- make two exons
    gene.map <- rbind(gene.map,
                      data.table(start = c(exon.starts[i],coding.end),
                                 end = c(coding.end,exon.ends[i]),
                                 ymin = c(-0.5,-0.25),
                                 ymax = c(0.5,0.25),
                                 annotation = c("cds","utr")
                      ))
  } else {
    # All other exons <- make one exon
    gene.map <- rbind(gene.map,
                      data.table(start = exon.starts[i],
                                 end = exon.ends[i],
                                 ymin = -0.5,
                                 ymax = 0.5,
                                 annotation = "cds"
                      ))
  }
  
  len.exon = (exon.ends[i] - exon.starts[i]) + 20
  last.fake = last.pos + len.exon
  pos.map <- rbind(pos.map,
                   data.table(pos=(exon.starts[i]-10):(exon.ends[i]+10),
                              fake.pos=last.pos:last.fake))
  last.pos = last.fake + 1
}

# Then we merge that "fake" map on top of the actual gene map
gene.map <- merge(gene.map, pos.map, by.x = "start", by.y = "pos")
setnames(gene.map,"fake.pos","fake.start")
gene.map <- merge(gene.map, pos.map, by.x = "end", by.y = "pos")
setnames(gene.map,"fake.pos","fake.end")

variants <- merge(variants, pos.map, by.x="POS",by.y="pos", all.x = T)
variants <- variants[!is.na(fake.pos)]
variants[,log.p:=-log10(P_BOLT_LMM_INF)]
variants[,mod.p:=if_else(BETA < 0, T, F)]

if (max(variants[,log.p], na.rm = T) > plot_ymax) {
  warning(paste0("A variant in ", gene_info[,SYMBOL], " has a log10 p. value greater than ", plot_ymax, " (max -log10 P = ",  variants[,sprintf("%0.1f", max(log.p,na.rm = T))], "; change with -y / --ymax) and will not be plotted..."))
}

fake.coding.start <- pos.map[pos == coding.start, fake.pos]
fake.coding.end <- pos.map[pos == coding.end, fake.pos]

if (!is.na(phenotype)) {
  x_axis_name <- paste(gene_info[,SYMBOL], mask, maf, phenotype, "BOLT", sep = " - ")
} else {
  x_axis_name <- paste(gene_info[,SYMBOL], mask, maf, "BOLT", sep = " - ")
}

gene.plot <- ggplot() + 
  geom_segment(aes(x = fake.coding.start, xend = fake.coding.end, y = 0, yend = 0),size = 1) +
  geom_hline(yintercept = c(1.5 + -log10(pvalue), -1.5 - -log10(pvalue)),colour="red",linetype=2) +
  geom_segment(data = variants, aes(x = fake.pos, xend = fake.pos, y = 0, yend = if_else(mod.p == T, -1.5-log.p, 1.5+log.p))) +
  geom_segment(aes(x = fake.coding.start - 60, xend = fake.coding.start - 60, y = 2.5, yend = 8.5), arrow = arrow(length = unit(0.02, "npc")), size = 1, lineend = 'round', linejoin = 'round', colour = "darkgrey") +
  geom_segment(aes(x = fake.coding.start - 60, xend = fake.coding.start - 60, y = -2.5, yend = -8.5), arrow = arrow(length = unit(0.02, "npc")), size = 1, lineend = 'round', linejoin = 'round', colour = "darkgrey") +
  annotate(geom = "text", x = fake.coding.start - 40, y = 5.5, hjust = 0, label = "Pos. β", size = 5) +
  annotate(geom = "text", x = fake.coding.start - 40, y = -5.5, hjust = 0, label = "Neg. β", size = 5) +
  geom_rect(data = gene.map[annotation!="utr"], aes(xmin = fake.start, xmax = fake.end, ymin = ymin, ymax = ymax), fill = "lightblue", colour = "black") + 
  geom_point(data = variants, aes(fake.pos, if_else(mod.p == T, -1.5-log.p, 1.5+log.p), size = BOLT_AC)) +
  scale_x_continuous(name = x_axis_name, limits = c(fake.coding.start - 75, fake.coding.end + 50), breaks = c(fake.coding.start,fake.coding.end)) +
  scale_y_continuous(name = expression(bold(-log[10](italic(p)))), breaks = plot_ybreaks, labels = plot_ylabels, limits = plot_ylimits) +
  scale_size_continuous(guide=guide_legend(title = "Allele Count")) +
  coord_flex_cart(bottom=capped_horizontal(capped="both"),left=capped_vertical(capped="both")) + # Comes from the "lemon" package
  theme

if (nrow(variants[log.p > -log10(pvalue)]) > 0) {
  gene.plot <- gene.plot + geom_text(data = variants[log.p > -log10(pvalue)], aes(fake.pos, if_else(mod.p == T, -1.5-log.p, 1.5+log.p), label = paste0("AC = ", BOLT_AC)), nudge_x = 40, size = 5, hjust = 0)
}

ggsave(filename = paste(file_prefix,"png",sep="."), plot = gene.plot, width = 20, height = 8, dpi = 450)
fwrite(variants, file = paste(file_prefix,"tsv",sep="."), col.names = T, row.names = F, quote = F, na = "NA", sep = "\t")
  
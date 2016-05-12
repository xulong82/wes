library(ggplot2)

n.var <- sapply(1:22, function(x) { cat(x, "\n")
  load(paste0("/data/xwang/wes/geno/gMix_chr", x, ".rdt"))
  sapply(gMix, nrow)
})  # SAVED

n.var.golden <- sapply(1:22, function(x) { cat(x, "\n")
  load(paste0("/data/xwang/wes/golden/meta_chr", x, ".rdt"))
  nrow(meta)
})  # SAVED

gdt <- data.frame(num = c(n.var), type = rep(rownames(n.var), 22), chr = rep(1:22, each = 4))

gdt$num <- gdt$num * 1e-3
gdt$chr <- factor(gdt$chr, levels = 1:22)
gdt$type <- factor(gdt$type, levels = rownames(n.var))

mycol <- c("grey70", "firebrick1", "dodgerblue3", "chartreuse3")
pdf("~/Dropbox/GitHub/wes/gMix_summary.pdf", width = 8, height = 4, family = "Helvetica")

ggplot(gdt, aes(x = chr, y = num, fill = type)) + 
  geom_bar(stat = "identity", width = 0.65) +
  theme_bw() + xlab("") + ylab("Variant Number in Kilo") +
  scale_fill_manual(values = mycol)
dev.off()

n.var.maf = rbind(bi = colSums(n.var[c("bSnpInf", "bIndInf"), ]), golden = n.var.golden)

n.var.maf = rbind(pass = n.var.maf["golden", ], fail = n.var.maf[1, ] - n.var.maf[2, ])
gdt <- data.frame(num = c(n.var.maf), type = rep(c("MAF>0.01", "MAF<0.01"), 22), chr = rep(1:22, each = 2))

gdt$num <- gdt$num * 1e-3
gdt$chr <- factor(gdt$chr, levels = 1:22)
gdt$type <- factor(gdt$type, levels = c("MAF>0.01", "MAF<0.01"))

pdf("~/Dropbox/GitHub/wes/maf_summary.pdf", width = 8, height = 4, family = "Helvetica")
ggplot(gdt, aes(x = chr, y = num, fill = type)) + 
  geom_bar(stat = "identity", width = 0.65) + 
  theme_bw() + xlab("Chromosome") + ylab("Variant Number in Kilo") +
  scale_fill_manual(values = c("dodgerblue3", "firebrick1"))
dev.off()

# ---

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr_gene = c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand")
gene_ensembl = getBM(attr_gene, filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = human)
attr_exon = c(att_gene, "exon_chrom_start", "exon_chrom_end")
exon_ensembl = getBM(attr_exon, filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = human)

gene_gr = with(exon_ensembl, GRanges(chromosome_name, IRanges(start=start_position, end=end_position)))
(gene_gr_total = sum(width(disjoin(gene_gr))))
exon_gr = with(exon_ensembl, GRanges(chromosome_name, IRanges(start=exon_chrom_start, end=exon_chrom_end)))
(exon_gr_total = sum(width(disjoin(exon_gr))))
(intron_gr_total = gene_gr_total - exon_gr_total)
(genome_length = width(range(gene_gr)) %>% as.numeric %>% sum)
(intergenic_gr_total = genome_length - gene_gr_total)

exon_var = unlist(data$vep_UID_exon)
intron_var = unlist(data$vep_UID_intron)
intron_var = setdiff(intron_var, exon_var)

(exon_var_length = length(exon_var))
(intron_var_length = length(intron_var))
(intergenic_var_length = 12.6e6 - exon_var_length - intron_var_length)

(intergenic_density = intergenic_var_length / intergenic_gr_total)
(intron_density = intron_var_length / intron_gr_total)
(exon_density = exon_var_length / exon_gr_total)


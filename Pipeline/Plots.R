#Read in files, create DDS object using DESeq2

dir <- ("~/new_Kallisto_Quant/CYTO/")
samples <- read.table(file.path(dir, "meta.cyto.txt"), 
                        sep="\t", header=T, row.names="Sample")
                        
files <- file.path(dir, rownames(samples), "abundance.h5")
names(files) <- paste0(rownames(samples))

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id"), mart = mart)
tx2gene <- results[, 1:2]
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

samples$rep <- as.factor(samples$rep)
samples$condition <- as.factor(samples$condition)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ rep + condition )
dds$condition <- relevel(dds$condition, ref = "PA") 
keep <- rowSums(counts(dds)) >= 32 #remove rows with counts less than 20
dds <- dds[keep,]
dds <- DESeq(dds) 

# Log2 transformation
cts <- counts(dds, normalized=T)
log_counts <- log2(cts + 1)

# sample to sample heatmap
num_conditions <- nlevels(samples$Condition)
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(samples$condition)]
pdf(file="~/Kallisto_Quant/Control/control_sample_heatmap.pdf",width=9)
heatmap.2(cor(log_counts), RowSideColors=cond_colors, margins=c(10,10), cexCol = 0.2 + 1/log10(16),
          labCol = F, dendrogram = "row", key = F,
          trace='none', main='Control Sample correlations (log2-transformed)')
dev.off()

# Pvclust dendogram
result <- pvclust(log_counts, method.dist = "cor", method.hclust = "average", nboot = 1000)
pdf(file="~/new_Kallisto_Quant/CYTO/cyto_sample_hclust.pdf")
plot(result)
pvrect(result)
dev.off()

# Stacked PCA / Correlation plots
p1 <- biplot(p, x="PC1", y="PC2",
  colby = 'condition', colkey = c('ARIA'='royalblue', 'TYTO'='red1', 'ARIA_CD362'='forestgreen', "PA"="purple"),
  hline = 0, vline = 0,
  #legendPosition = 'right', legendLabSize = 12, legendIconSize = 8.0,
  drawConnectors = TRUE,
  title = 'Control PCA bi-plot',
  subtitle = 'PC1 versus PC2')

p2 <- biplot(p, x="PC2", y="PC3",
  colby = 'condition', colkey = c('ARIA'='royalblue', 'TYTO'='red1', 'ARIA_CD362'='forestgreen', "PA"="purple"),
  hline = 0, vline = 0,
  #legendPosition = 'right', legendLabSize = 12, legendIconSize = 8.0,
  drawConnectors = TRUE,
  title = 'Control PCA bi-plot',
  subtitle = 'PC2 versus PC3')

p3 <- biplot(p, x="PC4", y="PC3",
  colby = 'condition', colkey = c('ARIA'='royalblue', 'TYTO'='red1', 'ARIA_CD362'='forestgreen', "PA"="purple"),
  hline = 0, vline = 0,
  #legendPosition = 'right', legendLabSize = 12, legendIconSize = 8.0,
  drawConnectors = TRUE,
  subtitle = 'PC4 versus PC3')

p4 <- biplot(p, x="PC4", y="PC2",
  colby = 'condition', colkey = c('ARIA'='royalblue', 'TYTO'='red1', 'ARIA_CD362'='forestgreen', "PA"="purple"),
  hline = 0, vline = 0,
  #legendPosition = 'right', legendLabSize = 12, legendIconSize = 8.0,
  drawConnectors = TRUE,
  subtitle = 'PC4 versus PC2')
p5 <- eigencorplot(p, components = getComponents(p, 1:10), metavars=c('condition', 'rep'))
pdf(file="~/new_Kallisto_Quant/Control/PCA.pdf", width=10, height=12)
grid.arrange(arrangeGrob(p1, p2, p3, p4, ncol=2, nrow=2),
             arrangeGrob(p5 ,ncol = 1, nrow=1), heights=c(4/5, 1/5))
dev.off()

# DESeq2 RESULTS
res <- results(dds, filterFun=ihw, alpha=0.1, c("condition", "ARIA", "PA"))

# Function for Outputs
generate_output <- function(x, prefix){
  up_key <-  intersect(rownames(x)[which(x$log2FoldChange>=1)],
                       rownames(x)[which(x$pvalue<=0.1)])
  
  up_df <- as.data.frame((x)[which(rownames(x) %in% up_key),])
  up_df <- tibble::rownames_to_column(up_df, "ensembl_gene_id")
  
  info <- getBM(attributes=c("ensembl_gene_id",
                        "external_gene_name",
                        "chromosome_name",
                        "start_position",
                        "end_position",
                        "description"),
                        filters = c("ensembl_gene_id"),
                        values = up_df$ensembl_gene_id,
                        mart = mart)
  
  tmp <- merge(info, up_df, by="ensembl_gene_id")
  tmp <- tmp[order(-tmp$log2FoldChange),]
  write.table(tmp, paste0(prefix,"_up_reg.txt"),
              sep="\t", row.names = F, quote=F)
  
  down_key <- intersect(rownames(x)[which(x$log2FoldChange<=-1)],
                        rownames(x)[which(x$pvalue<=0.1)])
  
  down_df <- as.data.frame((x)[which(rownames(x) %in% down_key),])
  down_df <- tibble::rownames_to_column(down_df, "ensembl_gene_id")
  
  info <- getBM(attributes=c("ensembl_gene_id",
                        "external_gene_name",
                        "chromosome_name",
                        "start_position",
                        "end_position",
                        "description"),
                         filters = c("ensembl_gene_id"),
                         values = down_df$ensembl_gene_id,
                         mart = mart)
  
  tmp <- merge(info, down_df, by="ensembl_gene_id")
  tmp <- tmp[order(tmp$log2FoldChange),]
  write.table(tmp, paste0(prefix,"_down_reg.txt"), sep="\t", row.names = F, quote=F)
  
  sig_genes <- c(up_key, down_key)
  write.table(sig_genes, paste0(prefix,"_sig_list.txt"), 
              row.names = F, col.names = F, quote=F)
}

# call output
generate_output(res, "aria_vs_pa")

# heatmap
de_gene <- read.table("~/new_Kallisto_Quant/Control/aria-cd362_vs_pa_sig_list.txt", sep = "\t", header=F)
#counts <- read.table("~/Results/norm_counts.txt", header=T, sep="\t", row.names = "ensembl_gene_id")
key <- as.character(de_gene[,1])
mat <- as.data.frame((log_counts)[which(rownames(log_counts) %in% key),])
mat <- tibble::rownames_to_column(mat, var="gene_id")
mat <- mat  %>% dplyr::select(-c("H88-ARIA", "H88-TYTO", "H89-ARIA", "H89-TYTO", "H90-ARIA", "H90-TYTO", "H91-ARIA", "H91-TYTO"))
#mat <- tibble::rownames_to_column(mat, var="gene_id")
info <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                        filters = "ensembl_gene_id",
                        values =mat$gene_id,
                        mart = mart)

rownames(mat) <- make.names(info$external_gene_name, unique=TRUE)
mat <- subset(mat, select=-c(1))

pdf("~/new_Kallisto_Quant/Control/aria_cd362_vs_pa_heatmap.pdf",height=16)
pheatmap(mat, color=greenred(75), cluster_rows = T,
         fontsize=9, 
         show_rownames = T, scale = "row")
dev.off()

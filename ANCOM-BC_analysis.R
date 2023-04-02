library(phyloseq)
library(ANCOMBC)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(ggplotify)


#================Build a Phyloseq-Class Object from Scratch==================
########基于丰度前300的种
otu_mat <- read.delim("D:/Desktop/返稿分析/species.txt",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
otu_mat <- otu_mat[order(rowSums(otu_mat),decreasing = T),]
otu_mat <- otu_mat[2:301,]

colnames(otu_mat_group) <- c("early","middle1","middle2","late")
meta <- read.delim("D:/Desktop/返稿分析/group.txt",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
tax_mat <- read.delim("D:/Desktop/返稿分析/species-tax.txt",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
tax_mat <- as.matrix(tax_mat[rownames(otu_mat),])

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
META = sample_data(meta)
TAX = tax_table(tax_mat)
physeq = phyloseq(OTU, META, TAX)

out = ancombc(phyloseq = physeq, formula = "group",
              p_adj_method = "BH", zero_cut = 0.90, lib_cut = 1000,
              group = "group", struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = TRUE)

otures_global <- out$res_global
otures_global <- otures_global[otures_global$diff_abn == "TRUE",]
x <- out[["res"]][["beta"]]
x <- x[rownames(otures_global),]
#write.csv(x,"x.csv")

samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
log_corr_abn <- log_corr_abn[order(rowSums(log_corr_abn),decreasing = T),]

pheat <- log_corr_abn[rownames(otures_global),]
pheat <- pheat[order(rowSums(pheat),decreasing = T),]

meta$group <- factor(meta$group,levels = c("early","middle-1","middle-2","late"))
p1 <- as.ggplot(pheatmap(pheat,scale = "row",legend_breaks = c(-4,-2,0,2,4),
                         cluster_rows = F,cluster_cols = F,show_colnames = F,
                         color = colorRampPalette(colors = c("blue","white","red"))(100),
                         angle_col="45",cellwidth = 10,
                         border_color = NA,gaps_col = c(6,12,19),annotation_col = meta))

bar_ancomba <- data.frame(cbind(rownames(pheat),rowSums(pheat)))
colnames(bar_ancomba) <- c("Species","Abundance")
bar_ancomba$Abundance <- as.numeric(bar_ancomba$Abundance)
bar_ancomba$Species <- factor(bar_ancomba$Species,levels = bar_ancomba$Species)
bar_ancomba <- cbind(bar_ancomba,tax_mat[rownames(bar_ancomba),])
p2 <- ggplot(bar_ancomba,aes(Abundance,Species,fill=phylum)) + geom_bar(stat = "identity")+
  theme(panel.background = element_blank())

plot_grid(p1,p2)
setwd("D:/Desktop/")
ggsave("wwwwwwww.pdf",width = 20,height = 5)

write.csv(otures_global,"otures_global.csv")

plot_grid(p1,p2)
setwd("D:/Desktop/")
ggsave("wwwwwwww.pdf",width = 20,height = 5)

write.csv(otures_global,"otures_global.csv")


########基于keggL3
otu_mat <- read.delim("D:/Desktop/返稿分析/Pathway.l3.sample.txt",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
otu_mat_re <- sweep(otu_mat,2,colSums(otu_mat),"/")
otu_mat <- otu_mat[order(rowSums(otu_mat_re),decreasing = T),]
meta <- read.delim("D:/Desktop/返稿分析/group.txt",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
META = sample_data(meta)
physeq = phyloseq(OTU, META)

out = ancombc(phyloseq = physeq, formula = "group",
              p_adj_method = "BH", zero_cut = 0.90, lib_cut = 1000,
              group = "group", struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = TRUE)

otures_global <- out$res_global
otures_global <- otures_global[otures_global$diff_abn == "TRUE",]
x <- out[["res"]][["beta"]]
x <- x[rownames(otures_global),]
#write.csv(x,"x.csv")

samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
log_corr_abn <- log_corr_abn[order(rowSums(log_corr_abn),decreasing = T),]

pheat <- log_corr_abn[rownames(otures_global),]
pheat <- pheat[order(rowSums(pheat),decreasing = T),]

meta$group <- factor(meta$group,levels = c("early","middle-1","middle-2","late"))
p1 <- as.ggplot(pheatmap(pheat,scale = "row",legend_breaks = c(-4,-2,0,2,4),
                         cluster_rows = F,cluster_cols = F,show_colnames = F,
                         color = colorRampPalette(colors = c("blue","white","red"))(100),
                         angle_col="45",cellwidth = 10,
                         cellheight = 20,
                         border_color = NA,gaps_col = c(6,12,19),annotation_col = meta))
bar_ancomba <- data.frame(cbind(rownames(pheat),rowSums(pheat)))
colnames(bar_ancomba) <- c("Species","Abundance")
bar_ancomba$Abundance <- as.numeric(bar_ancomba$Abundance)
bar_ancomba$Species <- factor(bar_ancomba$Species,levels = bar_ancomba$Species)
p2 <- ggplot(bar_ancomba,aes(Abundance,Species)) + geom_bar(stat = "identity")+
  theme(panel.background = element_blank())
plot_grid(p1,p2)
ggsave("ancom-keggl3.pdf",width = 20,height = 8)
write.csv(otures_global,"keggl3_otures_global.csv")




library(ggplot2)
library(cowplot)
library(vegan)
library(reshape2)
library(ggpubr)
library(ggplotify)
library(ggthemes)
library(pairwiseAdonis)
library(FSA)
library(rcompanion)
library(tidyr)
compaired <- list(c("Dec","Jan"),c("Dec","Feb"),c("Dec","Mar"),c("Jan","Feb"),c("Jan","Mar"),c("Feb","Mar"))

###Taxonomic composition α-diversity and β-diversity############################################
#α-diversity
setwd("D:/Desktop/返稿分析")
data <- read.delim('species.sample.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
data <- data.frame(t(data[,-1]))
Shannon.Wiener=diversity(data,"shannon")
Simpson=diversity(data,"simpson")
da <- cbind(Shannon.Wiener,Simpson)
da <- cbind(rownames(da),da)
da[grep(pattern = "A",rownames(da)),1] <- "Early"
da[grep(pattern = "B",rownames(da)),1] <- "Middle1"
da[grep(pattern = "C",rownames(da)),1] <- "Middle2"
da[grep(pattern = "D",rownames(da)),1] <- "Late"
da <- data.frame(da)
colnames(da)[1] <- "Months"
da$Months <- factor(da$Months,levels = c("Early","Middle1","Middle2","Late"))
da_melt <- melt(da,id = "Months")
colnames(da_melt)[2] <- "class"

da_melt_st <- da_melt[25:48,]
kruskal.test(value~Months,da_melt_st)

du <- dunnTest(as.numeric(value)~Months,da_melt_st,method="bh")
#write.csv(du$res,"qqqqqqqq.csv")

p1 <- ggplot(da_melt, aes(x=Months, y=as.numeric(value),fill=Months)) + 
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=1,color="black",size=0.5)+
  scale_fill_manual(values = c("#0e72cc","#38cb7d","#fa2c7b","#05f8d6"))+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",width = 0.15,position = position_dodge( .9))+
  facet_wrap(class~.,scales = "free_y", ncol = 1)+
  stat_compare_means(method = "kruskal")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = 1),
        panel.spacing.x = unit(0.3,"cm"),
        axis.title = element_blank(),
        axis.line= element_blank(),
        axis.text.x = element_text(size = 10,angle = 45),
        axis.text.y = element_text(size = 12),
        legend.key = element_blank(),
        strip.background = element_rect(color = "black",size = 1),
        strip.text = element_text(size = 10),
        strip.placement = "inside",
        aspect.ratio = 1.3)+guides(fill="none")
p1
#β-diversity
otu <- read.delim('species.sample.relative.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
otu <- otu[,-1]
data <- as.data.frame(t(otu))
group <- data.frame(rownames(data))
group[grep(group$rownames.data.,pattern = "A",),2] <- "Early"
group[grep(group$rownames.data.,pattern = "B",),2] <- "Middle1"
group[grep(group$rownames.data.,pattern = "C",),2] <- "Middle2"
group[grep(group$rownames.data.,pattern = "D",),2] <- "Late"
adonis <- adonis2(data ~ V2, data = group, permutations = 999, method="bray")
adonis_lab <- paste0("PERMANOVA R2=",round(adonis$R2[1],2),",p=", adonis$`Pr(>F)`[1])
otu <- data.frame(t(otu))
distance <- vegdist(otu, method = 'bray') #bray  /  euclidean欧式
group <- factor(c(rep(1,6),rep(2,6),rep(3,7),rep(4,5)), labels = c("early","middle1","middle2","late"))
dispersion <- betadisper(distance, group, bias.adjust = FALSE,
                         sqrt.dist = FALSE, add = FALSE)
x <- permutest(dispersion)
betadisper_lab <- paste0("PERMDISP R2=",round(x$tab[1,2]/colSums(x$tab)[2],2),",p=",x$tab$`Pr(>F)`[1])
all_lab <- paste(betadisper_lab,"\n",adonis_lab)
pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)
plot_data <- data.frame({pcoa$point})[1:2]
plot_data$Months <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
plot_data[grep(pattern = "A",rownames(plot_data)),3] <- "Early"
plot_data[grep(pattern = "B",rownames(plot_data)),3] <- "Middle1"
plot_data[grep(pattern = "C",rownames(plot_data)),3] <- "Middle2"
plot_data[grep(pattern = "D",rownames(plot_data)),3] <- "Late"
plot_data$Months <- factor(plot_data$Months,levels = c("Early","Middle1","Middle2","Late"))
eig = pcoa$eig
p2 <- ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2,color = Months))+geom_point(alpha=1,size=3)+
  scale_color_manual(values = c("#0e72cc","#38cb7d","#fa2c7b","#05f8d6"))+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +ggtitle(all_lab)+
  stat_ellipse(level=0.95,show.legend=F,size=0.5,alpha=1)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = 1),
        panel.spacing.x = unit(0.3,"cm"),
        plot.title = element_text(hjust = 1, vjust=-15),
        legend.position = "bottom",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        legend.key = element_blank(),
        strip.background = element_rect(color = "black",size = 1),
        strip.text = element_text(size = 15),
        strip.placement = "inside",
        aspect.ratio = 1.1)+guides(fill="none")
p2
plot_grid(p1,p2,align = "hv",labels = "AUTO", nrow = 1, axis = 'l')
ggsave("组成多样性.pdf",width = 7, height = 7)

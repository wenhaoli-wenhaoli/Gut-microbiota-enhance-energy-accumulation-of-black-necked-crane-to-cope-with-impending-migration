library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(ggpubr)
library(ggplotify)
library(ggthemes)

##Ggualluvial_analysis
data <- read.delim("phylum.sample.relative.txt",row.names = 1,sep="\t",stringsAsFactors = F,check.names = F)
data$sum <- rowSums(data[,-1])
data <- data[order(data$sum,decreasing = T),]
x <- rownames(data)[c(1:2,4:11)]
data_top <- data[x,-c(1,ncol(data))]
plist <- list()
for(i in 1:dim(data_top)[2]){
  data_test <- t(data_top[i,])
  data_test <- data.frame(cbind(rownames(data_test),data_test))
  data_test[grep("A",rownames(data_test)),1] <- "Dec"
  data_test[grep("B",rownames(data_test)),1] <- "Jan"
  data_test[grep("C",rownames(data_test)),1] <- "Feb"
  data_test[grep("D",rownames(data_test)),1] <- "Mar"
  y <- colnames(data_test)[2]
  colnames(data_test) <- c("months","phylum")
  data_test$months <- factor(data_test$months,levels = c("Dec","Jan","Feb","Mar"))
  data_test$phylum <- as.numeric(data_test$phylum)*100
  plist[[i]] <- ggplot(data_test,aes(months,phylum,fill=months)) + geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=1,color="black",size=0.5)+
    stat_summary(fun.data = 'mean_sd',geom = "errorbar",colour = "black",width = 0.15,position = position_dodge( .9))+
    theme_classic()+labs(title = y,x="",y="")+stat_compare_means(method="wilcox.test",comparisons =compaired,label="p.signif",hide.ns = T)+theme(panel.background = element_rect(fill = 'white', color = 'black'),plot.title = element_text(hjust = 0.5))+guides(fill="none")
}
cowplot::plot_grid(plist[[1]],plist[[2]],plist[[3]],plist[[4]],plist[[5]],plist[[6]],plist[[7]],plist[[8]],plist[[9]],plist[[10]],nrow = 2)
#ggsave("组成top10门图.pdf",path = "D:/Desktop/8.19图/",width = 10, height = 6)
data <- data.frame(t(data[c(1,2,4:11),-c(1,ncol(data))]))
mean_dec <- c()
sd_dec <- c()
mean_jan <- c()
sd_jan <- c()
mean_feb <- c()
sd_feb <- c()
mean_mar <- c()
sd_mar <- c()

for(i in 1:10){
  dec_mean <- round(mean(data[grep("A",rownames(data)),i])*100,2)
  mean_dec <- c(mean_dec,dec_mean)
  dec_sd <- round(sd(data[grep("A",rownames(data)),i])*100,2)
  sd_dec <- c(sd_dec,dec_sd)
  jan_mean <- round(mean(data[grep("B",rownames(data)),i])*100,2)
  mean_jan <- c(mean_jan,jan_mean)
  jan_sd <- round(sd(data[grep("B",rownames(data)),i])*100,2)
  sd_jan<- c(sd_jan,jan_sd)
  feb_mean <- round(mean(data[grep("C",rownames(data)),i])*100,2)
  mean_feb <- c(mean_feb,feb_mean)
  feb_sd <- round(sd(data[grep("C",rownames(data)),i])*100,2)
  sd_feb<- c(sd_feb,feb_sd)
  mar_mean <- round(mean(data[grep("D",rownames(data)),i])*100,2)
  mean_mar <- c(mean_mar,mar_mean)
  mar_sd <- round(sd(data[grep("D",rownames(data)),i])*100,2)
  sd_mar<- c(sd_mar,mar_sd)
  print(paste(colnames(data)[i],"12月平均数±标准差:",dec_mean,"±",dec_sd,
              "1月平均数±标准差:",jan_mean,"±",jan_sd,
              "2月平均数±标准差:",feb_mean,"±",feb_sd,
              "3月平均数±标准差:",mar_mean,"±",mar_sd))
}
dec_adu <- paste(mean_dec,"±",sd_dec,sep = "")
jan_adu <- paste(mean_jan,"±",sd_jan,sep = "")
feb_adu <- paste(mean_feb,"±",sd_feb,sep = "")
mar_adu <- paste(mean_mar,"±",sd_mar,sep = "")
phylum_a <- data.frame(colnames(data),dec_adu,jan_adu,feb_adu,mar_adu)
colnames(phylum_a) <- c("Phylum","Dec","Jan","Feb","Mar")
#write.csv(phylum_a,"phylum_aundance.csv",row.names = F,quote = F)

data <- read.delim("phylum.group.relative.txt",row.names = 1,sep = "\t",stringsAsFactors = F,check.names = F)
data_top <- data[x,-1] #取topn的子集，若其中存在未分类，则n = n + 1
data_top['Others', ] <- 1 - colSums(data_top)
data_top <- data_top[,c(1,3,2,4)] #########注意未分类所在行
data_top$Phylum <- factor(rownames(data_top),levels = rev(c("Others",rownames(data_top)[-11])))
colnames(data_top) <- c("Dec","Jan","Feb","Mar","Phylum")
data_top <- melt(data_top, id = 'Phylum')
p1 <- ggplot(data_top, aes(x=variable,y=value*100,alluvium = Phylum)) +
  geom_alluvium(aes(fill=Phylum),alpha=1,size=0) +
  scale_fill_manual(values = c("#008080","#0000ff","#c8cc00","#FA3C3C","#f59311","#85c021","#ff38e0","#0e72cc","#38cb7d","#fa2c7b","#05f8d6"))+
  ylab("Relative abundance (%)")+xlab("")+ggtitle("Phylum changes among months")+theme_classic()

colnames(data)[2:5] <- c("Dec","Feb","Jan","Mar")
data <- data[x,-1]
data <- data[,c(1,3,2,4)]
ph1<- as.ggplot(pheatmap(data,fondsize=10,scale ="row",legend_breaks = c(-1, 0, 1),cluster_rows = T,cluster_cols = T,border_color=NA,
                         color = colorRampPalette(colors = c("blue","white","red"))(100),
                         clustering_method = "average",angle_col="0",cellwidth = 40))



data <- read.delim("species.sample.relative.txt",row.names = 1,sep="\t",stringsAsFactors = F,check.names = F)
data$sum <- rowSums(data[,-1])
data <- data[order(data$sum,decreasing = T),]
x <- rownames(data)[2:26]
data_top <- data[x,-c(1,ncol(data))]
plist <- list()
for(i in 1:dim(data_top)[2]){
  data_test <- t(data_top[i,])
  data_test <- data.frame(cbind(rownames(data_test),data_test))
  data_test[grep("A",rownames(data_test)),1] <- "Dec"
  data_test[grep("B",rownames(data_test)),1] <- "Jan"
  data_test[grep("C",rownames(data_test)),1] <- "Feb"
  data_test[grep("D",rownames(data_test)),1] <- "Mar"
  y <- colnames(data_test)[2]
  colnames(data_test) <- c("months","genus")
  data_test$months <- factor(data_test$months,levels = c("Dec","Jan","Feb","Mar"))
  data_test$genus <- as.numeric(data_test$genus)*100
  plist[[i]] <- ggplot(data_test,aes(months,genus,fill=months)) + geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=1,color="black",size=0.5)+
    stat_summary(fun.data = 'mean_sd',geom = "errorbar",colour = "black",width = 0.15,position = position_dodge( .9))+
    theme_classic()+labs(title = y,x="",y="")+stat_compare_means(method="wilcox.test",comparisons =compaired,label="p.signif",hide.ns = T)+theme(panel.background = element_rect(fill = 'white', color = 'black'),plot.title = element_text(hjust = 0.5))+guides(fill="none")
}
cowplot::plot_grid(plist[[1]],plist[[2]],plist[[3]],plist[[4]],plist[[5]],plist[[6]],plist[[7]],plist[[8]],plist[[9]],plist[[10]],nrow = 2)
#ggsave("组成top10种图.pdf",path = "D:/Desktop/8.19图/",width = 10, height = 6)
cowplot::plot_grid(plist[[11]],plist[[12]],plist[[13]],plist[[14]],plist[[15]],plist[[16]],plist[[17]],plist[[18]],plist[[19]],plist[[20]],nrow = 2)
#ggsave("组成top10-20种图.pdf",path = "D:/Desktop/8.19图/",width = 10, height = 6)
data <- data.frame(t(data[2:21,-c(1,ncol(data))]))
mean_dec <- c()
sd_dec <- c()
mean_jan <- c()
sd_jan <- c()
mean_feb <- c()
sd_feb <- c()
mean_mar <- c()
sd_mar <- c()

for(i in 1:20){
  dec_mean <- round(mean(data[grep("A",rownames(data)),i])*100,2)
  mean_dec <- c(mean_dec,dec_mean)
  dec_sd <- round(sd(data[grep("A",rownames(data)),i])*100,2)
  sd_dec <- c(sd_dec,dec_sd)
  jan_mean <- round(mean(data[grep("B",rownames(data)),i])*100,2)
  mean_jan <- c(mean_jan,jan_mean)
  jan_sd <- round(sd(data[grep("B",rownames(data)),i])*100,2)
  sd_jan<- c(sd_jan,jan_sd)
  feb_mean <- round(mean(data[grep("C",rownames(data)),i])*100,2)
  mean_feb <- c(mean_feb,feb_mean)
  feb_sd <- round(sd(data[grep("C",rownames(data)),i])*100,2)
  sd_feb<- c(sd_feb,feb_sd)
  mar_mean <- round(mean(data[grep("D",rownames(data)),i])*100,2)
  mean_mar <- c(mean_mar,mar_mean)
  mar_sd <- round(sd(data[grep("D",rownames(data)),i])*100,2)
  sd_mar<- c(sd_mar,mar_sd)
  print(paste(colnames(data)[i],"12月平均数±标准差:",dec_mean,"±",dec_sd,
              "1月平均数±标准差:",jan_mean,"±",jan_sd,
              "2月平均数±标准差:",feb_mean,"±",feb_sd,
              "3月平均数±标准差:",mar_mean,"±",mar_sd))
}
dec_adu <- paste(mean_dec,"±",sd_dec,sep = "")
jan_adu <- paste(mean_jan,"±",sd_jan,sep = "")
feb_adu <- paste(mean_feb,"±",sd_feb,sep = "")
mar_adu <- paste(mean_mar,"±",sd_mar,sep = "")
species_a <- data.frame(colnames(data),dec_adu,jan_adu,feb_adu,mar_adu)
colnames(species_a) <- c("Species","Dec","Jan","Feb","Mar")
#write.csv(species_a,"species_aundance.csv",row.names = F,quote = F)
data <- read.delim("species.group.relative.txt",row.names = 1,sep = "\t",stringsAsFactors = F,check.names = F)
data_top <- data[x[1:20],-1] #取topn的子集，若其中存在未分类，则n = n + 1
data_top['Others', ] <- 1 - colSums(data_top)
data_top <- data_top[,c(1,3,2,4)] #########注意未分类所在行
data_top$Species <- factor(rownames(data_top),levels = rev(c("Others",rownames(data_top)[-21])))
colnames(data_top) <- c("Dec","Jan","Feb","Mar","Species")
data_top <- melt(data_top, id = 'Species')
p3 <- ggplot(data_top, aes(x=variable,y=value*100,alluvium = Species)) +
  geom_alluvium(aes(fill=Species),alpha=1,size=0) +
  scale_fill_manual(values = c("#8EA0C9","#16afcc","#99cc00","#e30039","#fcd300","#00994e","#800080","#ff6600","#808000","#db00c2",
                               "#008080","#0000ff","#c8cc00","#FA3C3C","#f59311","#85c021","#ff38e0","#0e72cc","#38cb7d","#fa2c7b",
                               "#05f8d6"))+
  ylab("Relative abundance (%)")+xlab("")+ggtitle("Species changes among months")+theme_classic()+
  guides(fill=guide_legend(keywidth=1,keyheight=0.5,ncol=1,byrow=TRUE))+theme(legend.spacing.y = unit(0,units = "pt"))
#ggsave("种流积图.pdf",path = "D:/Desktop/8.19图/",width = 6, height = 5)

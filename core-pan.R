###core-pan cruve

library(ggplot2)

data <- read.delim2("D:/毕业/数据分析/3.gene/gene_table.taxa.txt",check.names = F,stringsAsFactors = F)
data <- data[-c(1,26)]

##pan-core gene

pan_data <- data.frame(sample="",gene="")
core_data <- data.frame(sample="",gene="")
sample_number <- 1:dim(data)[2]
for (i in 1:dim(data)[2]) {
  x <- t(combn(sample_number,i))
  num1 <- c()
  num2 <- c()
  if (dim(x)[1]>100){
    x <- x[sample(1:dim(x)[1],100),]
  }
  for (j in 1:dim(x)[1]) {
    y <- length(which(rowSums(data[x[j,]])!=0))
    num1 <- c(num1,y)
    test <- data[x[j,]]
    test[test == 0] <- NA
    z <- dim(test <- na.omit(test))[1]
    num2 <- c(num2,z)
  }
  tempdata <- data.frame(i,num1)
  colnames(tempdata) <- c("sample","gene")
  pan_data <- rbind(pan_data,tempdata)
  tempdata1 <- data.frame(i,num2)
  colnames(tempdata1) <- c("sample","gene")
  core_data <- rbind(core_data,tempdata1)
  
}

core_data <- core_data[-1,]
core_data <- data.frame(core_data,"core")
pan_data <- pan_data[-1,]
pan_data <- data.frame(pan_data,"pan")
colnames(core_data)[3] <- "type"
colnames(pan_data)[3] <- "type"
plot_data <- rbind(pan_data,core_data)
plot_data$sample <- factor(plot_data$sample,levels = 1:24)
plot_data$gene <- as.numeric(plot_data$gene)

ggplot(plot_data,aes(sample,gene,fill=type))+geom_boxplot(aes(group=sample))+
  facet_wrap(.~type,scales = "free")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = .8)
  



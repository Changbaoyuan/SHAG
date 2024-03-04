
setwd("E:\\BaiduSyncdisk\\0.extralwork\\2.RCCgxy")
load("RCCgxy.RData")
#rm(list = ls())
library("tidyverse")

tumor.path <- "E:\\BaiduSyncdisk\\0.extralwork\\2.RCCgxy";setwd(tumor.path) #create dir
data.path   <- file.path(tumor.path, "InputData")
fig1.path    <- file.path(tumor.path, "Figure1")
fig2.path    <- file.path(tumor.path, "Figure2")
fig3.path    <- file.path(tumor.path, "Figure3")
fig4.path    <- file.path(tumor.path, "Figure4")
fig5.path    <- file.path(tumor.path, "Figure5")
fig6.path    <- file.path(tumor.path, "Figure6")
fig7.path    <- file.path(tumor.path, "Figure7")
Scripts    <- file.path(tumor.path, "Scripts")

if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(data.path)) { dir.create(data.path) }
if (!file.exists(fig1.path)) { dir.create(fig1.path) }
if (!file.exists(fig2.path)) { dir.create(fig2.path) }
if (!file.exists(fig3.path)) { dir.create(fig3.path) }
if (!file.exists(fig4.path)) { dir.create(fig4.path) }
if (!file.exists(fig5.path)) { dir.create(fig5.path) }
if (!file.exists(fig6.path)) { dir.create(fig6.path) }
if (!file.exists(fig7.path)) { dir.create(fig7.path) }

# load R package
library(sva)
library(ConsensusClusterPlus)
library(pheatmap)
library(corrplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Boruta)
library(org.Hs.eg.db)
library(enrichplot)
library(ggalluvial)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(maftools)
library(GSVA)
#library(ComplexHeatmap)
library(gplots)
library(clusterProfiler)
library(tidyr)
library(ggplot2)
library(estimate)
library(ggpubr)
library('progress')
source(file.path(Scripts,"twoclasslimma.R"))
# load R package
library(NMF)
library(survival)
library(survminer)
library(sva)
library(Rtsne)
library(ComplexHeatmap)
library(gplots)
source(file.path(Scripts,"batchPCA.R"))
# set color

jco <- c("#FEB139","#18978F","#EE81B3","#001E6C","#990000")
jco2<-c("#FF6666", "#FFFF00", "#3399CC")
jco3<-c("#FFCC33","#009999","#CC0066")
jco4<-c("#146FB5","#E74811","#00A16D")
jco5<-c("#218D81","#F7BB1D","#040908")
yjp <- c("#F2B1BA","#EF8A09","#E24A19","#01764A","#0259A0")

tt<-jco4
pie(rep(1, length(tt)), col=tt)

#load R package

########## Figure 1 combat

#####GEO去批次

####combat of GEO datasets######
KIRC.expr<-read.table(file.path(data.path,"KIRC.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
KIRC.expr<-round(log2(KIRC.expr+1),2)
KIRC.expr<-KIRC.expr[!apply(KIRC.expr,1,function(x) length(x[x<1]))>(0.9*ncol(KIRC.expr)),]
KIRC.clin<-read.table(file.path(data.path,"KIRC.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
merge<-intersect(colnames(KIRC.expr),rownames(KIRC.clin))
KIRC.expr<-KIRC.expr[,merge]
KIRC.clin<-KIRC.clin[merge,]
KIRC.clin$OS.time<-KIRC.clin$OS.time/30.5
range(KIRC.expr)

EMTAB3267.expr<-read.table(file = file.path(data.path,"EMTAB3267.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
EMTAB3267.expr<-EMTAB3267.expr %>% 
  drop_na()
EMTAB3267.clin<-read.table(file = file.path(data.path,"EMTAB3267.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
EMTAB3267.expr<-EMTAB3267.expr[,rownames(EMTAB3267.clin)]
EMTAB3267.clin$OS.time<-EMTAB3267.clin$OS.time/30.5
range(EMTAB3267.expr)


GSE22541.expr<-read.table(file = file.path(data.path,"GSE22541.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE22541.clin<-read.table(file = file.path(data.path,"GSE22541.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE22541.expr<-GSE22541.expr[,rownames(GSE22541.clin)]
GSE22541.expr<-log2(GSE22541.expr+1)
GSE22541.clin$OS.time<-GSE22541.clin$OS.time/30.5
range(GSE22541.expr)

library(sva)
library(cluster)
library(oompaBase)

# 
comgene <- intersect(rownames(KIRC.expr),  rownames(GSE22541.expr))
# 
combined.expr <- cbind.data.frame(KIRC.expr[comgene,],
                                  GSE22541.expr[comgene,])
# 
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("KIRC","GSE22541"), times = c(ncol(KIRC.expr),ncol(GSE22541.expr))),
         fig.dir = fig1.path,
         PCA.fig.title = "Raw PCA for combined expression profile",
         cols = jco4[1:2],
         showID = F,
         cex = 0.7,
         showLegend = T) # 

range(combined.expr)

# 
batch <- data.frame(batch = rep(c("KIRC","GSE22541"), c(ncol(KIRC.expr),ncol(GSE22541.expr))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))


#
batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("KIRC","GSE22541"), times = c(ncol(KIRC.expr),ncol(GSE22541.expr))),
         fig.dir = fig1.path,
         PCA.fig.title = "Combat PCA for combined expression profile",
         cols = jco4[1:2],
         showID = F,
         cex = 0.7,
         showLegend = T) #


KIRC.clin$cohort<-"KIRC"
GSE22541.clin$cohort<-"GSE22541"

GSE22541.expr<-combined.expr.combat[,rownames(GSE22541.clin)]
KIRC.expr<-combined.expr.combat[,rownames(KIRC.clin)]

#########find the key pathways

library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#START GSVA

(load("hallmark.gs.RData")) 

gsym.expr <- KIRC.expr
head(gsym.expr)

# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(gsym.expr), gs)
head(gsva_es)

write.csv(gsva_es, file = file.path(data.path,"gsva_output.csv"), quote = F)

group<-KIRC.clin
group<- group[order(group$OS, decreasing = F), ]
table(group$OS)

group_list <- data.frame(sample = rownames(group), group = c(rep("a", 355), rep("b", 173)))
head(group_list)

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(a-b, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

write.csv(x, file = file.path(data.path,"gsva_limma.csv"), quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, file = file.path(data.path,"easy_input2_for39bar.csv"), quote = F, row.names = F)

#开始画图
df <- read.csv(file = file.path(data.path,"easy_input2_for39bar.csv"))
head(df)

cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('#18978F', '#BBBBBB', '#CC0066'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "inward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, Dead vs. Alive")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴

ggsave(file.path(fig1.path,"gsva.pdf"), width = 6, height = 8)

###########select hormone gene
hormone<-read.table(file = file.path(data.path,"pathways.gmt"),header=T,sep="\t",row.names=1,check.names=F)
hormone<-hormone[,-c(1:1)]

test<-t(hormone)

tt<-matrix(test,dimnames=list(t(outer(colnames(test),rownames(test),FUN=paste)),NULL))

tt[tt==""]<-NA
tt<-as.data.frame(tt)
tt<-tt %>% 
  drop_na()
select<-tt$V1[!duplicated(tt$V1)]


#################UpsetVenn plot 筛选可用的Histon基因#################
library("ggvenn")
# Default plot

KIRC.expr2<-read.table(file.path(data.path,"KIRC.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
KIRC.expr2<-KIRC.expr2[!apply(KIRC.expr2,1,function(x) length(x[x<1]))>(0.9*ncol(KIRC.expr2)),]

x <- list(
  TCGA.KIRC=rownames(KIRC.expr2),  
  GSE22541=rownames(read.table(file = file.path(data.path,"GSE22541.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)),
  hormone_gene=select
)

ggvenn(x)

library(UpSetR)
upset(fromList(x),order.by = "freq",nsets = 4,
      queries = list(list(query = intersects, params = list("TCGA.KIRC","GSE22541","hormone_gene"), active = T)))

dd2<-fromList(x)

#devtools::install_github("PhDMeiwp/basicPackages@master", force = TRUE) 安装下面的包时使用
library(basicPackages)
#basicPackages::install.yyplot()安装下面的包时使用

library(yyplot)
p2 <- yyplot::ggvenn(dd2)+
  theme_void()+
  theme(legend.position = "right")
p2

library(ggplotify) #把别的图转为ggplot2
library(ggimage) # 组合图片
p1<-upset(fromList(x),order.by = "freq",nsets = 4,
          queries = list(list(query = intersects, params = list("TCGA.KIRC","GSE22541","hormone_gene"), active = T)))
g1<-as.ggplot(p1) # 转换为ggplot2图片

library(yyplot)

g5<-g1 + geom_subview(subview = p2 + theme_void(), x=.7, y=.7, w=.5, h=.5)
g5

ggsave(file.path(fig1.path,"upsetR Venn.pdf"),width = 10,height = 5)

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are KIRCepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

# load meta signature
hormone.sig <- gmt2list(file.path(data.path,"pathways.gmt"))
hormone.sig <- sapply(hormone.sig, function(x) setdiff(x,""))
hormone.class <- NULL
for (i in names(hormone.sig)) {
  tmp <- hormone.sig[i]
  for (j in tmp) {
    hormone.class <- rbind.data.frame(hormone.class,
                                      data.frame(gene = j,
                                                 path = i,
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
  }
}


# calculate GSVA enrichment score计算通路分数
hormone.score <- gsva(expr = as.matrix(KIRC.expr),
                      gset.idx.list = hormone.sig,
                      method = "gsva")
###保存计算出来的分数，备用，投稿的时候可能会要原始数据
write.table(hormone.score,file.path(fig1.path,"gsva enrichment score of hormone signature in tcga cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
range(hormone.score)

library(ClassDiscovery) #用于外部聚类

annCol <-KIRC.clin[,c("OS","Age","Gender","Race","laterality","Stage")]
#annCol$PSA<-ifelse(annCol$PSA=="unknow","unknow",
 #                  ifelse(annCol$PSA<10,"<=10",">10"))

head(annCol)
#annCol[is.na(annCol) | annCol == ""] <- "N/A"
##设定每个因子的颜色
annCol$OS<-ifelse(annCol$OS=="0","Alive","Dead")
table(annCol$Race)
annColors <- list(OS  = c("Alive" = "grey",
                              "Dead"   = "black"),
                  Age    =  c("white", "#990000"),
                  Gender  = c("Male" = "#009999",
                           "Female"   = "#CC0066"),
                  Race  = c("unknow" = "grey",
                           "Asian"   = "#34A12E",
                           "White"   = "#FF7E02",
                           "BoAA"   = "#53207C"),
                  laterality = c(
                    "Left"    = "#F6BF02",
                    "Right"    = "#00A16D",
                    "Bilateral"    = "#E24A19"),
                  Stage = c(
                    "Stage I"    = yjp[1],
                    "Stage II"    = yjp[2],
                    "Stage III"    = yjp[3],
                    "Stage IV"    = yjp[4],
                    "[Discrepancy]"    = "grey"))
annColors
#画图看一下基本情况
#pdf(file.path(fig1.path, "heatmap in tcga.pdf"), width = 10,height = 12)
pheatmap(hormone.score,
         scale = "row",
         color = c(
           "#2421BA","#63BFA5","#F8FCB4","#EC6146","#D41714"),
         annotation_col = annCol,
         #cutree_cols = 2,
         #kmeans_k  = 2,
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F,
         filename = "raw_heatmap.pdf")
#dev.off()

##聚类分组，调整，选取差别最明显的分组。可以尝试不同的k值
mat<-hormone.score

hcs <- hclust(distanceMatrix(as.matrix(mat), "manhattan"), "ward.D") # 请阅读distanceMatrix()以及hclust()，了解更多distance测度和linkage方法
#the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
# one of "ward.D", "ward.D2", "single", "complete", "average" 
hcg <- hclust(distanceMatrix(t(as.matrix(mat)), "maximum"), "ward.D") # 注意距离函数是针对列的，所以对行聚类要转置
group <- cutree(hcs,k=3)

# 增加一行annotation，及其配色
annCol$Clust2 <- paste0("C",group[rownames(annCol)])
annColors[["Clust2"]] <- c("C1"=jco3[1],"C2"=jco3[2],"C3"=jco3[3],"C4"="blue")
#head(annCol)
pdf(file.path(fig1.path, "heatmap in tcga.pdf"), width = 10,height = 6)
pheatmap(mat,
         scale = "row",
         color = c("#4BFDF8","black", "#F40000"),
         cutree_cols = 3,
         cluster_rows = hcg,
         cluster_cols = hcs,
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T,show_colnames = F,
         filename = "heatmap_with_outside_Cluster.pdf")
dev.off()

####给分组加一个前缀，123变成C1C2C3
KIRC.clin$Clust<-paste0("C",group[rownames(annCol)])
table(KIRC.clin$Clust,KIRC.clin$OS)

library(survival)
library("survminer")

#KIRC.clin$OS.time<-KIRC.clin$OS.time/30.5
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ Clust, data = KIRC.clin)
ggsurvplot(fit)


pdf(file=file.path(fig1.path,"KIRC three hormone subgroup survival.pdf"), width=5,height=4.2,onefile = FALSE)
ggsurvplot(fit, data = KIRC.clin,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_classic2(), 
           conf.int = F, 
           conf.int.style = "step",
           censor = T, 
           palette = jco3, #
           ylim = c(0,1),
           
           xlab = 'Time in months',
           
           risk.table.y.text.col = T, 
           risk.table.y.text = T, 
           
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           #fun = "event",
           
           pval = T,
           pval.coord = c(0, 0.15)
)

dev.off()


############twoclass limma 寻找两组间差异基因，即组蛋白修饰最强和最弱组间差异基因
source(file.path(Scripts,"twoclasslimma.R"))
condition<-as.matrix(KIRC.clin$Clust)
colnames(condition)<-c("condition")

subt <- data.frame(condition = condition,
                   row.names = rownames(KIRC.clin))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = KIRC.expr[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "C1", # name of treatment group
              ctrlVar  = "C3", # name of control group
              prefix   = "KIRC", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = fig1.path) # path for result

# extract group specific pathways
tmp1 <- read.table(file.path(fig1.path,"KIRC_limma_test_result.C1_vs_C3.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

###########筛选差异基因，筛选与之可调整

deg1up <- rownames(tmp1[which(tmp1$log2fc > 0.6 & tmp1$padj < 0.05),])
deg1dw <- rownames(tmp1[which(tmp1$log2fc < -0.6& tmp1$padj < 0.05),])

genes<-c(deg1up,deg1dw)
genes<-intersect(genes,comgene)

##保存
write.table(genes,file=file.path(fig1.path,"intersectgene.txt"),sep="\t",row.names=T,quote=F)

##画火山图

tmp1$change = ifelse(tmp1$pvalue < 0.01& abs(tmp1$log2fc) >= 0.6, 
                     ifelse(tmp1$log2fc> 0.6 ,'Up','Down'),
                     'Stable')
tmp1$symbol<-rownames(tmp1)
p <- ggplot(data = tmp1, 
            aes(x = tmp1$log2fc, 
                y = -log10(tmp1$pvalue), 
                colour=change,
                label = tmp1$symbol)) +
  geom_point(alpha=0.6, size=3) +
  scale_color_manual(values=jco4)+
  #xlim(c(-2.5, 2.5)) +
  geom_vline(xintercept=c(-0.6,0.6),lty=4,col='black',lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col='black',lwd=0.8) +
  labs(x='log2(fold change)',
       y='-log10 (p-value)',
       title='hormone modification associated genes')  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position=c(0.1,0.8), 
        legend.title = element_blank())
p
ggsave(file.path(fig1.path,"hormone volcano.pdf"),width = 5,height = 5)

#--------新型火山图---------
diff_express2<-tmp1
diff_express2$name<-rownames(diff_express2)
diff_express2$log2FoldChange<-diff_express2$log2fc
# Default plot
ggmaplot(diff_express2, main = expression("C1" %->% "C3"),
         fdr = 0.05, fc = 1.5157, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(diff_express2$name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())
ggsave(file.path(fig1.path,"hormone volcano new.pdf"),width = 6,height = 5)


#####差异基因通路富集
library("clusterProfiler")
options(connectionObserver = NULL)
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

#
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gene<-entrezIDs <- as.character(entrezIDs)

genelist_input<-tmp1
genelist_input$enzid<-mget(genelist_input$symbol, org.Hs.egSYMBOL2EG, ifnotfound=NA)

genelist_input<-subset(genelist_input,genelist_input$enzid!="NA")

geneList = genelist_input[,2]
names(geneList) = as.character(genelist_input[,7])
geneList = sort(geneList, decreasing = TRUE)
edo2 <- gseDO(geneList)
ridgeplot(edo2)

#---Go通路富集
Go<- enrichGO(gene = gene,
              keyType = "ENTREZID",
              OrgDb = org.Hs.eg.db, 
              pvalueCutoff =0.05, 
              qvalueCutoff = 0.5,
              ont="all",
              readable =T)
write.table(Go,file.path(fig1.path,file="GO-up.txt"),sep="\t",quote=F,row.names = F)

upsetplot(Go)
library(DOSE)

###保存点图
pdf(file=file.path(fig1.path,"Go-up.pdf"),width = 10,height = 6)
dotplot(Go,showCategory = 10,label_format = 70,split="ONTOLOGY") + ##label_format 改前面名称显示的长度
  facet_grid(ONTOLOGY~., scale='free') + #是否根据BP，CC，MF分类
  scale_color_continuous(low='#009966', high='#FF0034')+ #设置颜色
  aes(shape=GeneRatio > 0.04)#设置三角形圆形
dev.off()

###保存调控网络图
pdf(file=file.path(fig1.path,"Go-tree.pdf"),width = 10,height = 6)
Go2<- pairwise_termsim(Go)
treeplot(Go2)
dev.off()


###保存通路联系、关键基因图
pdf(file=file.path(fig1.path,"Go-circ.pdf"),width = 8,height = 5)
cnetplot(Go2, showCategory = 5,categorySize="count", colorEdge = TRUE)
dev.off()


#---KEGG分析
Kegg <- enrichKEGG(gene = gene, organism = "hsa",
                   pvalueCutoff =0.5, qvalueCutoff =0.5)
write.table(Kegg,file.path(fig1.path,"KEGG-up.txt"),sep="\t",quote=F,row.names = F)


#设定你只想展示的通路
intest<-c("Nucleocytoplasmic transport","Cell cycle")
#          "Signaling pathways regulating pluripotency of stem cells",
#          "Axon guidance","Cellular senescence",
#          "mTOR signaling pathway","Proteoglycans in cancer",
#          "Ubiquitin mediated proteolysis","Phagosome")

pdf(file=file.path(fig1.path,"KEGG-up.pdf"),width = 8,height = 3)
dotplot(Kegg, showCategory = 8)+ scale_color_continuous(low='#009999', high='#53207C')+ aes(shape=GeneRatio > 0.04)
dev.off()

Kegg2<- pairwise_termsim(Kegg)
Kegg2<-setReadable(Kegg2, 'org.Hs.eg.db', 'ENTREZID')

pdf(file=file.path(fig1.path,"KEGG-net.pdf"),width = 8,height = 4)
#test<-c("ECM-receptor interaction","Focal adhesion","PI3K-Akt signaling pathway")
cnetplot(Kegg2, showCategory = 5,categorySize="count",colorEdge = TRUE)
dev.off()


pdf(file="Kegg-net.pdf",width = 10,height = 8)
emapplot(Kegg2, showCategory = 30,color = "p.adjust")+ scale_color_continuous(low='yellow', high='red')
dev.off()

####保存富集通路树图
pdf(file.path(fig1.path,"Kegg-tree.pdf"),width = 10,height = 5)
treeplot(Kegg2,showCategory = 30)
dev.off()

upsetplot(Kegg2)

heatplot(Kegg2, showCategory=5)

#---HALLMARK分析
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gene<-entrezIDs <- as.character(entrezIDs)
library(msigdbr)
msigdbr_species() #支持的物种
#Hu_msigdbr <- msigdbr(species="Homo sapiens")
#head(Hu_msigdbr, 2) %>% as.data.frame

HALLMARK<- msigdbr(species="Homo sapiens",category="H") %>% 
  dplyr::select(gs_name, entrez_gene, gene_symbol)
head(HALLMARK)

HALLMARKresult<- enricher(gene,TERM2GENE=HALLMARK[,c(1,2)])

write.table(HALLMARKresult,file=file.path(fig1.path,"HALLMARKresult.txt"),sep="\t",quote=F,row.names = F)

pdf(file=file.path(fig1.path,"HALLMARKresult.pdf"),width = 6,height = 4)
dotplot(HALLMARKresult, showCategory = 10)+ scale_color_continuous(low='#009966', high='#FF0034')+ aes(shape=GeneRatio > 0.1)
dev.off()

##################################
#################
###计算每个队列中，2620个组蛋白差异基因的预后预测能力
library(survival)
rt=cbind(KIRC.clin[,c("OS.time","OS")],t(KIRC.expr[genes,]))
#rt=cbind(KIRC.clin[,c("OS.time","OS")],t(KIRC.expr))
KIRCoutTab=data.frame()
pb <- progress_bar$new(total=ncol(rt)-2)
for(i in colnames(rt[,3:ncol(rt)])){
  rt[,i]<- factor(ifelse(rt[,i]>median(rt[,i]),"high"," low"))
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  KIRCoutTab=rbind(KIRCoutTab,
                   cbind(id=i,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  pb$tick()
  Sys.sleep(0.05)
}
rownames(KIRCoutTab)<-KIRCoutTab$id
write.table(KIRCoutTab,file.path(fig1.path,"KIRCuniCox.txt"),sep="\t",row.names=F,quote=F)


rt=cbind(GSE22541.clin[,c("OS.time","OS")],t(GSE22541.expr[genes,]))
#rt=cbind(GSE22541.clin[,c("OS.time","OS")],t(GSE22541.expr))
GSE22541outTab=data.frame()
pb <- progress_bar$new(total=ncol(rt)-2)
for(i in colnames(rt[,3:ncol(rt)])){
  rt[,i]<- factor(ifelse(rt[,i]>median(rt[,i]),"high"," low"))
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  GSE22541outTab=rbind(GSE22541outTab,
                       cbind(id=i,
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR.95L=coxSummary$conf.int[,"lower .95"],
                             HR.95H=coxSummary$conf.int[,"upper .95"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  pb$tick()
  Sys.sleep(0.05)
}

rownames(GSE22541outTab)<-GSE22541outTab$id
write.table(GSE22541outTab,file.path(fig1.path,"GSE22541uniCox.txt"),sep="\t",row.names=F,quote=F)

rt=cbind(EMTAB3267.clin[,c("OS.time","OS")],t(EMTAB3267.expr[genes,]))
#rt=cbind(EMTAB3267.clin[,c("OS.time","OS")],t(EMTAB3267.expr))
EMTAB3267outTab=data.frame()
pb <- progress_bar$new(total=ncol(rt)-2)
for(i in colnames(rt[,3:ncol(rt)])){
  rt[,i]<- factor(ifelse(rt[,i]>median(rt[,i]),"high"," low"))
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  EMTAB3267outTab=rbind(EMTAB3267outTab,
                       cbind(id=i,
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR.95L=coxSummary$conf.int[,"lower .95"],
                             HR.95H=coxSummary$conf.int[,"upper .95"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  pb$tick()
  Sys.sleep(0.05)
}
rownames(EMTAB3267outTab)<-EMTAB3267outTab$id
write.table(EMTAB3267outTab,file.path(fig1.path,"EMTAB3267uniCox.txt"),sep="\t",row.names=F,quote=F)

#这里是上面把两个GEO分开去筛选预后基因，但是结果太少。所以最后把两个GEO合并计算
rtt1=cbind(EMTAB3267.clin[,c("OS.time","OS")],t(EMTAB3267.expr[genes,]))
rtt2=cbind(GSE22541.clin[,c("OS.time","OS")],t(GSE22541.expr[genes,]))
rt<-rbind(rtt1,rtt2)
#rt=cbind(GEO.clin[,c("OS.time","OS")],t(GEO.expr))
pb <- progress_bar$new(total=ncol(rt)-2)
GEOoutTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  rt[,i]<- factor(ifelse(rt[,i]>median(rt[,i]),"high"," low"))
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  GEOoutTab=rbind(GEOoutTab,
                  cbind(id=i,
                        HR=coxSummary$conf.int[,"exp(coef)"],
                        HR.95L=coxSummary$conf.int[,"lower .95"],
                        HR.95H=coxSummary$conf.int[,"upper .95"],
                        pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  pb$tick()
  Sys.sleep(0.05)
}
rownames(GEOoutTab)<-GEOoutTab$id
write.table(GEOoutTab,file.path(fig1.path,"GEOuniCox.txt"),sep="\t",row.names=F,quote=F)


##筛选预后基因
Risk1<-rownames(KIRCoutTab[which(KIRCoutTab$pvalue<0.05&KIRCoutTab$HR>1),])
Prot1<-rownames(KIRCoutTab[which(KIRCoutTab$pvalue<0.05&KIRCoutTab$HR<1),])

#survgene<-c(intersect(Risk2,Risk1),intersect(Prot2,Prot1))

#Risk2<-rownames(GEOoutTab[which(GEOoutTab$pvalue<0.01&GEOoutTab$HR>1),])
#Prot2<-rownames(GEOoutTab[which(GEOoutTab$pvalue<0.01&GEOoutTab$HR<1),])

#Risk3<-rownames(EMTAB3267outTab[which(EMTAB3267outTab$pvalue<0.05&EMTAB3267outTab$HR>1),])
#Prot3<-rownames(EMTAB3267outTab[which(EMTAB3267outTab$pvalue<0.05&EMTAB3267outTab$HR<1),])

Risk4<-rownames(GSE22541outTab[which(GSE22541outTab$pvalue<0.05&GSE22541outTab$HR>1),])
Prot4<-rownames(GSE22541outTab[which(GSE22541outTab$pvalue<0.05&GSE22541outTab$HR<1),])

#survgene<-c(intersect(intersect(Risk3,Risk1),Risk4),intersect(intersect(Prot3,Prot1),Prot4))
survgene<-c(intersect(Risk4,Risk1),intersect(Prot4,Prot1))
library("ggVennDiagram")
# Default plot
venn1<-ggVennDiagram(list(Risk1,Risk4),
                     label_alpha = 0,
                     category.names = c("TCGA-KIRC cohort","GSE22541 cohort"))+
  ggplot2::scale_fill_gradient(low="#32AEEC",high = "#EB262F")
venn1
ggsave(file.path(fig1.path,"Risky gene Venn.pdf"),width = 5,height = 3)

venn2<-ggVennDiagram(list(Prot1,Prot4),
                     label_alpha = 0,
                     category.names = c("TCGA-KIRC cohort","GSE22541 cohort"))+
  ggplot2::scale_fill_gradient(low="#299D92",high = "#732A7C")
venn2
ggsave(file.path(fig1.path,"Protective gene Venn.pdf"),width = 5,height = 3)

##############################
library(yyplot)
library(ggplot2)

KIRCoutTab$symbol<-rownames(KIRCoutTab)
KIRCoutTab$pvalue<-as.numeric(KIRCoutTab$pvalue)
KIRCoutTab$HR<-as.numeric(KIRCoutTab$HR)
KIRCoutTab$change = ifelse(KIRCoutTab$pvalue < 0.05, 
                           ifelse(KIRCoutTab$HR> 1 ,'Riskey','Protective'),
                           'Stable')
KIRCvolcano <- ggplot(data = KIRCoutTab, 
                      aes(x = KIRCoutTab$HR, 
                          y = -log10(KIRCoutTab$pvalue), 
                          colour=change,
                          label = KIRCoutTab$symbol)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c('#32AEEC', '#EB262F','#299D92'))+
  xlim(c(0, 3)) +
  geom_vline(xintercept=c(1),lty=4,col='black',lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col='black',lwd=0.8) +
  labs(x='Hazard Ratio',
       y='-log10 (p-value)',
       title='TCGA-KIRC cohort')  +
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5), 
    legend.position=c(0.8,0.3), 
    legend.title = element_blank())
KIRCvolcano

GSE22541outTab$symbol<-rownames(GSE22541outTab)
GSE22541outTab$pvalue<-as.numeric(GSE22541outTab$pvalue)
GSE22541outTab$HR<-as.numeric(GSE22541outTab$HR)
GSE22541outTab$change = ifelse(GSE22541outTab$pvalue < 0.05, 
                               ifelse(GSE22541outTab$HR> 1 ,'Riskey','Protective'),
                               'Stable')
GSE22541volcano <- ggplot(data = GSE22541outTab, 
                          aes(x = GSE22541outTab$HR, 
                              y = -log10(GSE22541outTab$pvalue), 
                              colour=change,
                              label = GSE22541outTab$symbol)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c('#E8863D','#732A7C', '#299D92'))+
  xlim(c(0, 4)) +
  geom_vline(xintercept=c(1),lty=4,col='black',lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col='black',lwd=0.8) +
  labs(x='Hazard Ratio',
       y='-log10 (p-value)',
       title='GSE22541 cohort')  +
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5), 
    legend.position=c(0.8,0.3), 
    legend.title = element_blank())
GSE22541volcano


####拼图
ggarrange(KIRCvolcano, GSE22541volcano,venn1,venn2,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

ggsave(file.path(fig2.path,"cohort HR volcano.pdf"),width = 8,height = 8)

###########
rt=cbind(annCol,t(KIRC.expr[survgene,rownames(annCol)]))
#----------热图---------
library(dplyr)
library(pheatmap)                   #引用包
head(rt)
rt<-rt[order(rt$OS,decreasing = F),]

rt2 <- rt %>% 
  dplyr::select(-c(1:7)) 

rt2 <- scale(t(rt2)) #scale标化

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(rt2), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

rt3 <- standarize.fun(rt2,halfwidth =1)

head(rt3[1:3,1:3])


cluster<- rt %>% 
  dplyr::select(c(1:6)) 


#绘制热图
pdf(file.path(fig2.path,"KIRC_heatmap_192genes.pdf"),height=8,width=10)

pheatmap(rt3,
         scale = "row",
         color = c("#05D5B3","black", "#F00A6D"),
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = cluster,
         annotation_colors = annColors,
         show_rownames = F,show_colnames = F)
dev.off()

###########GSE22541
rt=cbind(GSE22541.clin,t(GSE22541.expr[survgene,rownames(GSE22541.clin)]))
#----------热图---------
library(dplyr)
library(pheatmap)                   #引用包
head(rt)
rt<-rt[order(rt$OS,decreasing = F),]

rt2 <- rt %>% 
  dplyr::select(-c(1:5)) 

rt2 <- scale(t(rt2)) #scale标化

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(rt2), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

rt3 <- standarize.fun(rt2,halfwidth =1)

head(rt3[1:3,1:3])


cluster<- rt %>% 
  dplyr::select(c(2:4)) 
table(cluster$Type)
annColors2 <- list(OS  = c("0" = "grey",
                          "1"   = "black"),
                  Gender  = c("male" = "#009999",
                              "female"   = "#CC0066"),
                  Type = c(
                    "metastasis"    = "#F6BF02",
                    "primary"    = "#00A16D"))
#绘制热图
pdf(file.path(fig2.path,"GSE22541_heatmap_192genes.pdf"),height=8,width=5)

pheatmap(rt3,
         scale = "row",
         color = c("#32AEEC","black", "#FDE805"),
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = cluster,
         annotation_colors = annColors2,
         show_rownames = F,show_colnames = F)
dev.off()





########################
###在TCGA中训练预后预测模型

LASSO.input<-cbind(KIRC.clin[,1:2],t(KIRC.expr[survgene,]))

#-------------LASSO----------

df<-LASSO.input
dim(df)
head(df)
mydata<-df[,3:ncol(df)]
mydata<-as.matrix(mydata)

#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
#####set.seed是为了保证每次重复结果都一样。如果一次结果不理想，可以调整seed

df<-LASSO.input
dim(df)
head(df)
mydata<-df[,3:ncol(df)]
mydata<-as.matrix(mydata)

#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
#####set.seed是为了保证每次重复结果都一样。如果一次结果不理想，可以调整seed

#n<-round(runif(1, min = 0, max = 1000000),0)

set.seed(510181)
#做10倍交叉验证，算出lambda值
cvfit = cv.glmnet(mydata, Surv(df$OS.time,df$OS), 
                  family = "cox",
                  nfold=20) #10倍交叉验证
cvfit$lambda.min
plot(cvfit)


coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- colnames(mydata)[active.min]
geneids
index.min = coef.min[active.min]
index.min

fit <- glmnet(mydata, Surv(df$OS.time,df$OS), 
              family = "cox",
              nfold=100) #10倍交叉验证


####fit plot 美化图片

library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)

x <- coef(fit)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1] # extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm  


ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Lambda (log scale)") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  #annotate('text',x = -3.3,y=0.26,label='Optimal Lambda = 0.012',color='black')+
  guides(col=guide_legend(ncol = 2))

ggsave(filename = file.path(fig2.path,"LASSO.fit.pdf"), width = 8,height = 5)

## 准备数据LASSO.cvfit.
xx <- data.frame(lambda=cvfit[["lambda"]],cvm=cvfit[["cvm"]],cvsd=cvfit[["cvsd"]],
                 cvup=cvfit[["cvup"]],cvlo=cvfit[["cvlo"]],nozezo=cvfit[["nzero"]])
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo,' vars')

ggplot(xx,aes(ll,cvm,color=NZERO))+
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),width=0.05,size=0.5)+
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_point(size=2)+
  xlab("Log Lambda")+ylab('Partial Likelihood Deviance')+
  theme_bw(base_rect_size = 1.5)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12)))+
  scale_x_continuous(expand = c(0.02,0.02))+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+ #记得修改下方的最佳lambda值
  annotate('text',x = -6.8,y=13.2,label=paste0('Optimal Lambda = ',round(cvfit$lambda.min,3)),color='black')+
  guides(col=guide_legend(ncol = 2))
ggsave(filename = file.path(fig2.path,"LASSO.cvfit.pdf"), width = 7.5,height = 5)

cvfit$lambda.min #最佳lambda值

cvfit$lambda.1se #一倍SE内的更简洁的模型

# 输出基因顺序
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- colnames(mydata)[active.min]
geneids
index.min = coef.min[active.min]

combine<-cbind(geneids, index.min)
write.csv(combine,file = file.path(fig3.path,"gene_index.csv"))

# 将纳入signature的变量拟合成一个变量，作为nomogram的输入

geneids
index.min

signature <- as.matrix(df[, geneids]) %*% as.matrix(index.min) 
write.csv(signature,file = file.path(fig3.path,"KIRC-signature.csv"))

total_signature<-signature1<-signature

KIRC.os<-LASSO.input[,1:2]

KIRC_KM_input<-cbind(KIRC.os,signature)
KIRC_KM_input1<-cbind(KIRC.os,signature1)

#------------K-M

library(survival)
library("survminer")

KIRC_KM_input$signature<- factor(ifelse(KIRC_KM_input$signature>median(KIRC_KM_input$signature),"high"," low"))
#KIRC_KM_input$OS.time<-KIRC_KM_input$OS.time/30.5
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ signature, data = KIRC_KM_input)

cox <- coxph(Surv(OS.time, OS) ~ signature, data = KIRC_KM_input)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab=rbind(outTab,
             cbind(
               HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
outTab

HR <- paste("Hazard Ratio = ", round(outTab$HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(outTab$HR.95L,3), round(outTab$HR.95H,3), sep = " - "), sep = "")


pdf(file=file.path(fig3.path,"KIRC survival.pdf"), width=4,height=4.2,onefile = FALSE)
ggsurvplot(fit, data = KIRC_KM_input,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_classic2(), 
           conf.int = F, 
           conf.int.style = "step",
           censor = T, 
           palette = jco4, #
           ylim = c(0,1),
           
           xlab = 'Time in months',
           legend.title='riskscore', 
           legend.labs=c('Low','High'), 
           
           risk.table.y.text.col = T, 
           risk.table.y.text = T, 
           
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           #fun = "event",
           
           pval = paste(pval = ifelse(outTab$pvalue < 0.001, "p < 0.001", 
                                      paste("P = ",round(outTab$pvalue,3), sep = "")),
                        HR, CI, sep = "\n"),
           pval.coord = c(0, 0.15)
)

dev.off()

outTab


#----------ROC————————————

library(timeROC)
library(survival)

pdf(file=file.path(fig3.path,"TCGA HS ROC 1 3 5.pdf"), width=4,height=4.5)

ROC.DSST<-timeROC(T=KIRC_KM_input1$OS.time,#结局时间
                  delta=KIRC_KM_input1$OS,#生存结局
                  marker=KIRC_KM_input1$signature1,#预测变量
                  cause=1,#阳性结局赋值，比如死亡，复发的赋值
                  weighting="marginal",# 权重计算方法，marginal是默认值，采用km计算删失分布
                  times=c(12,36,30),# 时间点，选取10年和20年生存率
                  ROC = TRUE,
                  iid = TRUE
)
plot(ROC.DSST,time=12,col=jco2[1],title=FALSE,lwd=2)
plot(ROC.DSST,time=36,col=jco2[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=30,col=jco2[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC.DSST$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC.DSST$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC.DSST$AUC[3],3))),
       col=jco2,lwd=2,bty = 'n')
dev.off()

#####################Risk table
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

data <- cbind(LASSO.input,signature)
data[1:2, 1:4]

bestvars <-geneids
bestvars

# risk score，用于画顶部散点图
rs <- data$signature
names(rs) <- rownames(data)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)


# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=data[names(sort(rs)),'OS.time']/12*12,
                        s=data[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]
exp_data[1:2,1:4]

plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=1.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#EB5353","#187498")) + 
  
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = min(rs_data$rs), 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  geom_segment(aes(x=0,y=median(rs_data$rs),
                   xend=nrow(rs_data),
                   yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
  theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A


plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=1.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=1,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#187498","#EB5353"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B

tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=1) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#187498", 
                                    high="#EB5353", mid="#F7E0A8") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

KIRCriskmap<-plot_grid(plot.A, plot.B, plot.C,
          labels = c("", "",""), # 或者按顺序标注ABC
          rel_heights = c(1,1,2), # 3个图的比例
          #label_x=0,
          #label_y=1,
          align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)
KIRCriskmap

ggsave(file.path(fig3.path,paste0("KIRC riskmap figures.pdf")), width = 6, height = 10)



#-----------KIRC-临床信息森林图-----------
  library("tidyverse")
#KIRC.sig<-signature
Risk<-KIRC_KM_input1[rownames(KIRC.clin),"signature1"]
KIRCmulti<-cbind(KIRC.clin,Risk)
KIRCmulti$Risk<- factor(ifelse(KIRCmulti$signature>median(KIRCmulti$signature),"high-risk"," low-risk"))
KIRCmulti<-KIRCmulti[,c("OS.time","OS","Age","Gender","Race","Stage","Grade","Risk","Tumor_status")]
KIRCmulti$Grade<-ifelse(KIRCmulti$Grade=="G1","G1+G2",
                        ifelse(KIRCmulti$Grade=="G2","G1+G2",KIRCmulti$Grade))

library(survival)
library(survminer)
KIRCmulti[KIRCmulti=="unknow"]<-NA
KIRCmulti[KIRCmulti==""]<-NA
KIRCmulti[KIRCmulti=="[Discrepancy]"]<-NA

KIRCmulti<-KIRCmulti %>% 
  drop_na()
head(KIRCmulti)



#构建模型
model <- coxph( Surv(OS.time, OS) ~Age+Gender+Race+Stage+Grade+Tumor_status+Risk, data =KIRCmulti,na.action=na.exclude )
coxSummary = summary(model)
coxSummary


outTab=data.frame()
outTab=cbind(
  HR=coxSummary$conf.int[,"exp(coef)"],
  HR.95L=coxSummary$conf.int[,"lower .95"],
  HR.95H=coxSummary$conf.int[,"upper .95"],
  pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file=file.path(fig3.path,"TCGAmultiCox2.xls"),sep="\t",row.names=F,quote=F)

pdf(file=file.path(fig3.path,"KIRC TOTAL clinical forest.pdf"),onefile = FALSE,
    width = 7,             #图片的宽度
    height = 11,            #图片的高度
)

ggforest(model,
         data=KIRCmulti,
         main = "KIRC Cohort",
         cpositions = c(0.01,0.14,0.36), 
         fontsize = 1, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()

model

####
library(ggpubr)
library(rstatix)
library(survminer)

boxtt<-KIRCmulti[,c("signature","Stage")]
boxtt[boxtt=="unknow"]<-NA
boxtt[boxtt==""]<-NA
boxtt[boxtt=="[Discrepancy]"]<-NA

boxtt<-boxtt %>% 
  drop_na()

boxtt<-boxtt[order(boxtt$Stage, decreasing = F), ]


compare_means(signature ~ Stage,  data = boxtt,
              ref.group = ".all.", method = "t.test")

ggboxplot(boxtt, x = "Stage", y = "signature", 
          fill = "Stage", 
          #color = "Stage",
          legend = "none") +
  geom_jitter(width=0.15, alpha=0.6)+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(boxtt$signature[ which(!is.na(boxtt$Stage))]), linetype = 2)+
  stat_compare_means(method = "anova", label.y = 0.8)+
  scale_fill_brewer(palette="OrRd")+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")

ggsave(file.path(fig4.path,paste0("Stage signature subgroup.pdf")), width = 6, height = 4)

#
library(ggpubr)
library(rstatix)
library(survminer)

boxtt<-KIRCmulti[,c("signature","Grade")]
boxtt[boxtt=="GX"]<-NA
boxtt[boxtt==""]<-NA
boxtt[boxtt=="[Discrepancy]"]<-NA

boxtt<-boxtt %>% 
  drop_na()

boxtt<-boxtt[order(boxtt$Grade, decreasing = F), ]


compare_means(signature ~ Grade,  data = boxtt,
              ref.group = ".all.", method = "t.test")

ggboxplot(boxtt, x = "Grade", y = "signature", 
          fill = "Grade", 
          #color = "Grade",
          legend = "none") +
  geom_jitter(width=0.15, alpha=0.6)+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(boxtt$signature[ which(!is.na(boxtt$Grade))]), linetype = 2)+
  stat_compare_means(method = "anova", label.y = -2.3)+
  scale_fill_brewer(palette="YlGnBu")+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")

ggsave(file.path(fig4.path,paste0("Grade signature subgroup.pdf")), width = 6, height = 4)


#------森力图后multiROC
library("survival")
library("survminer")
library("riskRegression")
library(dplyr)

write.table(KIRCmulti,file=file.path(fig6.path,"KIRCmulti.txt"),sep="\t",row.names=T,quote=F)

KIRCmulti<-cbind(KIRC.clin,signature)
all_expdata<-cbind(KIRC.clin,signature)%>% 
  select(c(signature))

all_expdata<-t(all_expdata)

mark_expdata<-all_expdata

clindata<-dplyr::select(KIRCmulti,c("Stage","Gender","Grade","Age","Tumor_status"))
traitData<-dplyr::select(KIRCmulti,c("OS.time","OS"))

expdata = all_expdata[rownames(mark_expdata),]
rownames(clindata) = gsub("-", ".", rownames(clindata))
rownames(traitData) = gsub("-", ".", rownames(traitData))
colnames(traitData) = c("time", "status")
clindata2 = clindata
traitData2 = traitData
exp_sig = expdata

#Nomogram
predict_mat_all = na.omit(cbind(traitData2, clindata2, exp_sig))
cox_form = as.formula(paste("Surv(time, status)", paste("`", paste(colnames(predict_mat_all)[3:ncol(predict_mat_all)], collapse = "` + `"), "`", sep = ""), sep = " ~ "))
res.cox.all = coxph(cox_form, data = predict_mat_all, x = TRUE)
res.cox.lst = list(res.cox.all)

#Classifier
predict_mat_classifier = cbind(traitData2, exp_sig)[rownames(predict_mat_all),]
cox_form = as.formula(paste("Surv(time, status)", paste("`", paste(colnames(predict_mat_classifier)[3:ncol(predict_mat_classifier)], collapse = "` + `"), "`", sep = ""), sep = " ~ "))
res.cox.classifier = coxph(cox_form, data = predict_mat_classifier, x = TRUE)
res.cox.lst = c(res.cox.lst, list(res.cox.classifier))

# Single variate Cox regression ROC comparation
predict_mat_tmp = cbind(traitData2, clindata2)[rownames(predict_mat_all),]
predict_mat_tmp$Nomogram = predict(res.cox.all,predict_mat_all)
predict_mat_tmp$Classifier = predict(res.cox.classifier,predict_mat_classifier)
predict_mat_for_single = predict_mat_tmp[,c(1:2,ncol(predict_mat_tmp)-1,ncol(predict_mat_tmp),3:(ncol(predict_mat_tmp)-2))]
res.cox.lst = NULL
for (i in 3:(ncol(predict_mat_for_single))){
  test_mat_single = predict_mat_for_single[,c(1,2,i)]
  cox_form = as.formula(paste("Surv(time, status)", paste("`", colnames(predict_mat_for_single)[i], "`", sep = ""), sep = " ~ "))
  cox.res.test = coxph(cox_form, data = predict_mat_for_single, x = TRUE)
  if(is.null(res.cox.lst)){
    res.cox.lst = list(cox.res.test)
  }else{
    res.cox.lst = c(res.cox.lst, list(cox.res.test))
  }
}
names(res.cox.lst) = c(colnames(predict_mat_for_single)[3:ncol(predict_mat_for_single)])
xs <- Score(res.cox.lst, Hist(time,status)~1,data=predict_mat_for_single, plots="roc",metrics="auc")
pdf(file=file.path(fig3.path,"/total KIRC ROC combined.pdf", sep = ""), height = 6.5, width = 6)
plotROC(xs, xlab = "False negative rate", ylab = "Ture negative rate", legend = TRUE, auc.in.legend = TRUE)
dev.off()

#####################################
##################GSE22541 validation
#####################################

GSE22541.signature <- as.matrix(t(GSE22541.expr)[,geneids]) %*% as.matrix(index.min) 
write.csv(signature,file=file.path(fig3.path,"GSE22541 raw signature.csv"))

signature1<-signature<-GSE22541.signature

GSE22541.os<-GSE22541.clin[,c("OS.time","OS")]

GSE22541_KM_input<-cbind(GSE22541.os,signature)
GSE22541_KM_input1<-cbind(GSE22541.os,signature1)

#------------K-M

library(survival)
library("survminer")

GSE22541_KM_input$signature<- factor(ifelse(GSE22541_KM_input$signature>median(GSE22541_KM_input$signature),"high"," low"))
#GSE22541_KM_input$OS.time<-GSE22541_KM_input$OS.time/30.5
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ signature, data = GSE22541_KM_input)

cox <- coxph(Surv(OS.time, OS) ~ signature, data = GSE22541_KM_input)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab=rbind(outTab,
             cbind(
               HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
outTab

HR <- paste("Hazard Ratio = ", round(outTab$HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(outTab$HR.95L,3), round(outTab$HR.95H,3), sep = " - "), sep = "")


pdf(file=file.path(fig3.path,"GSE22541 survival.pdf"), width=4,height=4.2,onefile = FALSE)
ggsurvplot(fit, data = GSE22541_KM_input,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_classic2(), 
           conf.int = F, 
           conf.int.style = "step",
           censor = T, 
           palette = jco4, #
           ylim = c(0,1),
           
           xlab = 'Time in months',
           legend.title='riskscore', 
           legend.labs=c('Low','High'), 
           
           risk.table.y.text.col = T, 
           risk.table.y.text = T, 
           
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           #fun = "event",
           
           pval = paste(pval = ifelse(outTab$pvalue < 0.001, "p < 0.001", 
                                      paste("P = ",round(outTab$pvalue,3), sep = "")),
                        HR, CI, sep = "\n"),
           pval.coord = c(0, 0.2)
)

dev.off()

outTab


#----------ROC————————————

library(timeROC)
library(survival)

pdf(file=file.path(fig3.path,"GSE22541 ROC 1 3 5.pdf"), width=4,height=4.5)

ROC.DSST<-timeROC(T=GSE22541_KM_input1$OS.time,#结局时间
                  delta=GSE22541_KM_input1$OS,#生存结局
                  marker=GSE22541_KM_input1$signature1,#预测变量
                  cause=1,#阳性结局赋值，比如死亡，复发的赋值
                  weighting="marginal",# 权重计算方法，marginal是默认值，采用km计算删失分布
                  times=c(36,60,120),# 时间点，选取10年和20年生存率
                  ROC = TRUE,
                  iid = TRUE
)
plot(ROC.DSST,time=36,col=jco2[1],title=FALSE,lwd=2)
plot(ROC.DSST,time=60,col=jco2[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=120,col=jco2[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC.DSST$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC.DSST$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC.DSST$AUC[3],3))),
       col=jco2,lwd=2,bty = 'n')
dev.off()

#####################Risk table
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

data <- cbind(GSE22541.os,signature,t(GSE22541.expr[geneids,]))
data[1:2, 1:4]

bestvars <-geneids
bestvars

# risk score，用于画顶部散点图
rs <- data$signature
names(rs) <- rownames(data)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)


# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=data[names(sort(rs)),'OS.time']/12*12,
                        s=data[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]
exp_data[1:2,1:4]
plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=1.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#EB5353","#187498")) + 
  
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = min(rs_data$rs), 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  geom_segment(aes(x=0,y=median(rs_data$rs),
                   xend=nrow(rs_data),
                   yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
  theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A

plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=1.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=1,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#187498","#EB5353"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B

tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=1) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#187498", 
                                    high="#EB5353", mid="#F7E0A8") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

GSE22541risk<-plot_grid(plot.A, plot.B, plot.C,
                        labels = c("", "",""), # 或者按顺序标注ABC
                        rel_heights = c(1,1,2), # 3个图的比例
                        #label_x=0,
                        #label_y=1,
                        align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)


ggsave(file.path(fig3.path,paste0("GSE22541 riskmap figures.pdf")), width = 6, height = 10)

#-----------GSE22541-临床信息森林图-----------
  library("tidyverse")
#GSE22541.sig<-signature

GSE22541multi<-cbind(GSE22541.clin,signature)
GSE22541multi$Risk<- factor(ifelse(GSE22541multi$signature>median(GSE22541multi$signature),"high-risk"," low-risk"))

GSE22541multi[GSE22541multi==""]<-NA
GSE22541multi<-GSE22541multi %>% 
  drop_na()


library(survival)
library(survminer)
head(GSE22541multi)

#构建模型
model <- coxph( Surv(OS.time, OS) ~Gender+Risk, data =GSE22541multi,na.action=na.exclude )
coxSummary = summary(model)
coxSummary 
pdf(file=file.path(fig3.path,"GSE22541 TOTAL clinical forest.pdf"),onefile = FALSE,
    width = 5,             #图片的宽度
    height = 3,            #图片的高度
)

ggforest(model,
         data=GSE22541multi,
         main = "GSE22541 Cohort",
         cpositions = c(0.01,0.14,0.36), 
         fontsize = 1, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()

model

#------森力图后multiROC
library("survival")
library("survminer")
library("riskRegression")

write.table(GSE22541multi,file=file.path(fig3.path,"GSE22541multiROCinput.txt"),sep="\t",row.names=T,quote=F)

all_expdata<-dplyr::select(GSE22541multi, signature)
all_expdata<-t(all_expdata)

mark_expdata<-all_expdata

clindata<-dplyr::select(GSE22541multi,c("Gender","Type"))
traitData<-dplyr::select(GSE22541multi,c("OS.time","OS"))

expdata = all_expdata[rownames(mark_expdata),]
rownames(clindata) = gsub("-", ".", rownames(clindata))
rownames(traitData) = gsub("-", ".", rownames(traitData))
colnames(traitData) = c("time", "status")
clindata2 = clindata
traitData2 = traitData
exp_sig = expdata

#Nomogram
predict_mat_all = na.omit(cbind(traitData2, clindata2, exp_sig))
cox_form = as.formula(paste("Surv(time, status)", paste("`", paste(colnames(predict_mat_all)[3:ncol(predict_mat_all)], collapse = "` + `"), "`", sep = ""), sep = " ~ "))
res.cox.all = coxph(cox_form, data = predict_mat_all, x = TRUE)
res.cox.lst = list(res.cox.all)

#Classifier
predict_mat_classifier = cbind(traitData2, exp_sig)[rownames(predict_mat_all),]
cox_form = as.formula(paste("Surv(time, status)", paste("`", paste(colnames(predict_mat_classifier)[3:ncol(predict_mat_classifier)], collapse = "` + `"), "`", sep = ""), sep = " ~ "))
res.cox.classifier = coxph(cox_form, data = predict_mat_classifier, x = TRUE)
res.cox.lst = c(res.cox.lst, list(res.cox.classifier))

# Single variate Cox regression ROC comparation
predict_mat_tmp = cbind(traitData2, clindata2)[rownames(predict_mat_all),]
predict_mat_tmp$Nomogram = predict(res.cox.all,predict_mat_all)
predict_mat_tmp$Classifier = predict(res.cox.classifier,predict_mat_classifier)
predict_mat_for_single = predict_mat_tmp[,c(1:2,ncol(predict_mat_tmp)-1,ncol(predict_mat_tmp),3:(ncol(predict_mat_tmp)-2))]
res.cox.lst = NULL
for (i in 3:(ncol(predict_mat_for_single))){
  test_mat_single = predict_mat_for_single[,c(1,2,i)]
  cox_form = as.formula(paste("Surv(time, status)", paste("`", colnames(predict_mat_for_single)[i], "`", sep = ""), sep = " ~ "))
  cox.res.test = coxph(cox_form, data = predict_mat_for_single, x = TRUE)
  if(is.null(res.cox.lst)){
    res.cox.lst = list(cox.res.test)
  }else{
    res.cox.lst = c(res.cox.lst, list(cox.res.test))
  }
}
names(res.cox.lst) = c(colnames(predict_mat_for_single)[3:ncol(predict_mat_for_single)])
xs <- Score(res.cox.lst, Hist(time,status)~1,data=predict_mat_for_single, plots="roc",metrics="auc")
pdf(file=file.path(fig3.path,"/total GSE22541 ROC combined.pdf", sep = ""), height = 5.5, width = 5)
plotROC(xs, xlab = "False negative rate", ylab = "Ture negative rate", legend = TRUE, auc.in.legend = TRUE)
dev.off()
#-------------- KIRC single gene bestsep 每个基因预后

svdata <- cbind(KIRC.os,t(KIRC.expr[geneids,]))

dim(svdata) #

head(svdata[1:3,1:3])
#svdata$OS.time<-svdata$OS.time/30.5
library(survival)
library(survminer)

res.cut <- svdata[,1:2]
  
  for(i in colnames(svdata[,3:ncol(svdata)])){
    svdata[,i]<- factor(ifelse(svdata[,i]>median(svdata[,i]),"high"," low"))
    res.cut=cbind(res.cut,svdata[,i])
  }
  

res.cat <- res.cut
colnames(res.cat)<-colnames(svdata)
##统计作图
my.surv <- Surv(res.cat$OS.time, res.cat$OS)

pl<-list()
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  #i="SAA1"
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  ##计算HR以及95%CI
  ##修改分组参照
  group <- factor(group, levels = c(" low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  
  #只画出p value<=0.05的基因，如果不想筛选，就删掉下面这行
  #if (p.val>0.05) next
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  #按照基因表达量从低到高排序，便于取出分界表达量
  svsort <- svdata[order(svdata[,i]),]
  
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      
                      ggtheme = theme_classic2(), #想要网格就运行这行
                      conf.int = T, #不画置信区间，想画置信区间就把F改成T
                      conf.int.style = "step",#置信区间的类型，还可改为ribbon
                      censor = T, #不显示观察值所在的位置
                      palette = c(jco[2],jco[1]), #线的颜色对应高、低
                      ylim = c(0,1),
                      font.legend = 12,
                      font.main = c(14, "bold", "darkblue"),
                      font.x = c(14, "bold", "black"),
                      font.y = c(14, "bold", "black"),
                      font.tickslab = c(12, "plain", "black"),
                      xlab = 'Time (months)',
                      #legend.title='CTSZ', # 自定义图例的标题
                      legend.labs=c('Low','High'), # 自定义分组变量的名字
                      legend.title = i,#基因名写在图例题目的位置
                      
                      #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                      
                      #在图例上标出高低分界点的表达量，和组内sample数量
                      #legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                      #               paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      
                      #在左下角标出pvalue、HR、95% CI
                      #太小的p value标为p < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"),
                      pval.coord = c(0, 0.2)
  )
  
  #如果想要一个图保存为一个pdf文件，就把下面这行前面的“#"删掉
  #ggsave(paste0(i,".pdf"),width = 4,height = 4)
}

length(pl)

res <- arrange_ggsurvplots(pl, 
                           print = T,
                           ncol = 3, nrow = 4)#每页纸画几列几行

ggsave(file=file.path(fig3.path,"KIRC 33 genes single sur.pdf"),res,width = 12,height = 16)


######GSE22541 单基因预后
svdata <- cbind(GSE22541.os,t(GSE22541.expr[geneids,]))

dim(svdata) #

head(svdata[1:3,1:3])
#svdata$OS.time<-svdata$OS.time/30.5
library(survival)
library(survminer)

res.cut <- svdata[,1:2]

for(i in colnames(svdata[,3:ncol(svdata)])){
  svdata[,i]<- factor(ifelse(svdata[,i]>median(svdata[,i]),"high"," low"))
  res.cut=cbind(res.cut,svdata[,i])
}


res.cat <- res.cut
colnames(res.cat)<-colnames(svdata)
##统计作图
my.surv <- Surv(res.cat$OS.time, res.cat$OS)

pl<-list()
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  #i="SAA1"
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  ##计算HR以及95%CI
  ##修改分组参照
  group <- factor(group, levels = c(" low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  
  #只画出p value<=0.05的基因，如果不想筛选，就删掉下面这行
  #if (p.val>0.05) next
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  #按照基因表达量从低到高排序，便于取出分界表达量
  svsort <- svdata[order(svdata[,i]),]
  
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      
                      ggtheme = theme_classic2(), #想要网格就运行这行
                      conf.int = F, #不画置信区间，想画置信区间就把F改成T
                      conf.int.style = "step",#置信区间的类型，还可改为ribbon
                      censor = T, #不显示观察值所在的位置
                      palette =   c(jco[2],jco[1]), #线的颜色对应高、低
                      ylim = c(0,1),
                      font.legend = 12,
                      font.main = c(14, "bold", "darkblue"),
                      font.x = c(14, "bold", "black"),
                      font.y = c(14, "bold", "black"),
                      font.tickslab = c(12, "plain", "black"),
                      xlab = 'Time (months)',
                      #legend.title='CTSZ', # 自定义图例的标题
                      legend.labs=c('Low','High'), # 自定义分组变量的名字
                      legend.title = i,#基因名写在图例题目的位置
                      
                      #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                      
                      #在图例上标出高低分界点的表达量，和组内sample数量
                      #legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                      #               paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      
                      #在左下角标出pvalue、HR、95% CI
                      #太小的p value标为p < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"),
                      pval.coord = c(0, 0.2)
  )
  
  #如果想要一个图保存为一个pdf文件，就把下面这行前面的“#"删掉
  #ggsave(paste0(i,".pdf"),width = 4,height = 4)
}

length(pl)

res <- arrange_ggsurvplots(pl, 
                           print = T,
                           ncol = 3, nrow = 4)#每页纸画几列几行

ggsave(file=file.path(fig3.path,"GSE22541 33 genes single sur.pdf"),res,width = 12,height = 16)



##################################
############Figure 5 Nomogram#####
##################################

#########nomogram

library(survival)

pbc<-KIRCmulti[,c("OS.time","OS","Age","Stage","Tumor_status","signature")]
pbc[pbc=="[Discrepancy]"]<-NA
pbc[pbc==""]<-NA
pbc<-pbc %>% 
  drop_na()

pbccox <- coxph(formula = Surv(OS.time,OS) ~ Age+Stage+Tumor_status+signature, data = pbc)
pbccox

library(regplot)

regplot(pbccox,
        #对观测2的六个指标在列线图上进行计分展示
        observation=pbc[c("TCGA-VS-A9UJ-01"),], #也可以不展示
        points = T, #If FALSE the regression scores of each βx contribution are shown. Otherwise contributions are represented by a 0-100 "points" scale.
        plots = c("bean", #可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"
                  "bars"), #可选"no plot" "boxes" "bars" "spikes"
        subticks = TRUE, dencol="#F6BF02",boxcol="#18978F",obscol="red",spkcol="brown",
        #预测3年和5年的死亡风险，此处单位是day
        center = T,
        failtime = c(60,36,12),
        prfail = TRUE, #cox回归中需要TRUE
        rank="sd", #rank="range" is by the range of the βx's, and rank="sd" is by the standard deviation of the βx's. 
        showP = T, #是否展示统计学差异
        droplines = T,#观测2示例计分是否画线
        #colors = mycol, #用前面自己定义的颜色
        #rank=NULL, #根据统计学差异的显著性进行变量的排序
        #interval="confidence"
) #展示观测的可信区间

dev.copy2pdf(file=file.path(fig5.path, "nomogram_new.pdf"), width = 9,height = 7)

#------------K-M

library(survival)
library("survminer")
points<-pbccox$linear.predictors
KM_input<-cbind(pbc,points)

KM_input$points<- factor(ifelse(KM_input$points>median(KM_input$points),"high"," low"))
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ points, data = KM_input)

cox <- coxph(Surv(OS.time, OS) ~ points, data = KM_input)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab=rbind(outTab,
             cbind(
               HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
outTab

HR <- paste("Hazard Ratio = ", round(outTab$HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(outTab$HR.95L,3), round(outTab$HR.95H,3), sep = " - "), sep = "")


pdf(file=file.path(fig5.path,"nomogram survival.pdf"), width=4,height=4.2,onefile = FALSE)
ggsurvplot(fit, data = KM_input,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_linedraw(), 
           #fun = "event",
           conf.int = T, 
           #conf.int.style = "",
           censor = T, 
           palette = c( "#18978F","#CC0066"), #
           ylim = c(0,1),
           
           xlab = 'Time in months',
           legend.title='riskscore', 
           legend.labs=c('Low points','High points'), 
           
           risk.table.y.text.col = T, 
           risk.table.y.text = T, 
           
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           
           
           
           
           pval = paste(pval = ifelse(outTab$pvalue < 0.001, "p < 0.001", 
                                      paste("P = ",round(outTab$pvalue,3), sep = "")),
                        HR, CI, sep = "\n"),
           pval.coord = c(20, 0.2)
)

dev.off()

outTab

#----------ROC————————————

library(timeROC)
library(survival)
KM_input1<-cbind(pbc,points)
head(KM_input1)

pdf(file.path(fig5.path,"nomograme ROC2.pdf"), 5, 5)
ROC.DSST<-timeROC(T=KM_input1$OS.time,#结局时间
                  delta=KM_input1$OS,#生存结局
                  marker=KM_input1$points,#预测变量
                  cause=1,#阳性结局赋值，比如死亡，复发的赋值
                  weighting="marginal",# 权重计算方法，marginal是默认值，采用km计算删失分布
                  times=c(12,36,30),# 时间点，选取10年和20年生存率
                  ROC = TRUE,
                  iid = TRUE
)
plot(ROC.DSST,time=12,col=jco3[1],title=FALSE,lwd=2)
plot(ROC.DSST,time=36,col=jco3[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=30,col=jco3[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC.DSST$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC.DSST$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC.DSST$AUC[3],3))),
       col=jco3,lwd=2,bty = 'n')
dev.off()

#------------------Calibration---------------
#------------------Calibration---------------

library(ResourceSelection)
pbccox <- coxph(formula = Surv(OS.time,OS) ~Age+Stage+Tumor_status+signature, data = pbc)
library(survival)
library(regplot)
library(rms)

f1<-cph(formula = Surv(OS.time,OS) ~Age+Stage+Tumor_status+signature, data = pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 12) 
#参数m=50表示每组50个样本进行重复计算
cal1<-calibrate(f1, cmethod="KM", method="boot",u=12,m=100,B=1000) 

pbc$OS1<-ifelse(pbc$OS.time>12,"0",pbc$OS)
pbc$OS1<-as.numeric(pbc$OS1)

dat1 <- as.data.frame(cbind(pbc$OS1,f1$linear.predictors))
colnames(dat1)<-c("OS","model")

fullmodel_glm1 <- glm(OS ~ model, 
                      data = dat1, 
                      family = "binomial", 
                      control = list(maxit = 50))

p.hoslem1 <- hoslem.test(fullmodel_glm1$y, fitted(fullmodel_glm1), g=10)$p.value


f3<-cph(formula = Surv(OS.time,OS) ~Age+Stage+Tumor_status+signature, data = pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 36) 
#参数m=50表示每组50个样本进行重复计算
cal3<-calibrate(f3, cmethod="KM", method="boot",u=36,m=100,B=1000) 

pbc$OS3<-ifelse(pbc$OS.time>36,"0",pbc$OS)
pbc$OS3<-as.numeric(pbc$OS3)

dat3 <- as.data.frame(cbind(pbc$OS3,f3$linear.predictors))
colnames(dat3)<-c("OS","model")

fullmodel_glm3 <- glm(OS ~ model, 
                      data = dat3, 
                      family = "binomial", 
                      control = list(maxit = 50))

p.hoslem3 <- hoslem.test(fullmodel_glm3$y, fitted(fullmodel_glm3), g=10)$p.value

f5<-cph(formula = Surv(OS.time,OS) ~Age+Stage+Tumor_status+signature, data = pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 60) 
#参数m=50表示每组50个样本进行重复计算
cal5<-calibrate(f5, cmethod="KM", method="boot",u=60,m=100,B=1000) 

pbc$OS5<-ifelse(pbc$OS.time>60,"0",pbc$OS)
pbc$OS5<-as.numeric(pbc$OS5)

dat5 <- as.data.frame(cbind(pbc$OS5,f5$linear.predictors))
colnames(dat5)<-c("OS","model")

fullmodel_glm5 <- glm(OS ~ model, 
                      data = dat5, 
                      family = "binomial", 
                      control = list(maxit = 50))

p.hoslem5 <- hoslem.test(fullmodel_glm5$y, fitted(fullmodel_glm5), g=10)$p.value


pdf("calibration_compare.pdf",width = 6,height = 6)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#14B96A"),
     bty = "o", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#14B96A"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#14B96A"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FFA900"),
     xlim = c(0,1),ylim= c(0,1),col = c("#FFA900"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#FFA900"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#5C0A98"),
     xlim = c(0,1),ylim= c(0,1),col = c("#5C0A98"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#5C0A98"), pch = 16)

abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#14B96A","#FFA900","#5C0A98"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框

text(0.25,0.205,bquote("Hosmer-Lemeshow 1-year"~italic(P)~" = "~.(round(p.hoslem1,3))),adj = 0)
text(0.25,0.135,bquote("Hosmer-Lemeshow 3-year"~italic(P)~" = "~.(round(p.hoslem3,3))),adj = 0)
text(0.25,0.065,bquote("Hosmer-Lemeshow 5-year"~italic(P)~" = "~.(round(p.hoslem5,3))),adj = 0)

dev.off()


pdf(file=file.path(fig5.path,"calibration_compare 1-year.pdf"),width = 5,height = 5)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#14B96A"),
     bty = "o", #只画左边和下边框
     xlim = c(0.6,1),ylim= c(0.6,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#14B96A"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#14B96A"), pch = 16)
mtext("")


abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

text(0.60,0.63,bquote("Hosmer-Lemeshow 1-year"~italic(P)~" = "~.(round(p.hoslem1,3))),adj = 0)

dev.off()


pdf(file=file.path(fig5.path,"calibration_compare 3-year.pdf"),width = 5,height = 5)
plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FFA900"),
     bty = "o", #只画左边和下边框
     xlim = c(0.35,1),ylim= c(0.35,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#FFA900"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#FFA900"), pch = 16)
mtext("")


abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

text(0.45,0.40,bquote("Hosmer-Lemeshow 3-year"~italic(P)~" = "~.(round(p.hoslem3,3))),adj = 0)

dev.off()

pdf(file=file.path(fig5.path,"calibration_compare 5-year.pdf"),width = 5,height = 5)
plot(cal5,lwd = 2,lty = 0,errbar.col = c("#5C0A98"),
     bty = "o", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#5C0A98"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#5C0A98"), pch = 16)
mtext("")


abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

text(0.15,0.12,bquote("Hosmer-Lemeshow 5-year"~italic(P)~" = "~.(round(p.hoslem5,3))),adj = 0)

dev.off()

##-------------------

for (i in 3:5) {
  fit <- coxph(Surv(lenfol, fstat)~rcspline.eval(whas500$bmi,nk=i,inclx = T)+gender,data=whas500, x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} #nk保存最优的knots数目，且样条数目为nk-2。具体见参考文献“理论1.pdf"第2页2.3节公式1
}


#-------------DCA-------------------


library(rms)
library(rmda)

modul<- decision_curve(OS ~ Age+Stage+Tumor_status+signature,
                       data= pbc,
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95)
modul1<- decision_curve(OS ~ Age,
                        data= pbc,
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals = 0.95)
modul2<-decision_curve(OS ~ Stage,
                       data= pbc,
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95)
modul3<- decision_curve(OS ~ Tumor_status,
                        data= pbc,
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals = 0.95)

pbccox <- glm(OS ~ ., data = pbc)

list<-list(modul,modul1,modul2,modul3)
jcot <- c("#C72E24","#97B959","#F16A2F","#0088A2","#80A7DE","#411050")

pdf(file=file.path(fig5.path,"DCA.pdf"),width = 4,height = 4.5)
plot_decision_curve(list,
                    curve.names=c("Nomogram","Age","Stage","Tumor_status"),
                    xlab="Threshold probability",
                    cost.benefit.axis =FALSE,col= jcot,
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

modul


#-------------------clinical impact curve nomogram--------------

pdf(file=file.path(fig5.path,"Clinical_impact_curve.pdf"),width = 4,height = 4.5)
plot_clinical_impact(modul,xlim = c(0, 1),legend.position = "topright",
                     col = c("#EA4440", "#1C819E"))
dev.off()
#--------------------points correlation

KM_input1$points2<-(71+KM_input1$signature*31)+
  ifelse(KM_input1$Tumor_status=="Tumor Free",35,
         ifelse(KM_input1$Tumor_status=="With Tumor",66,71))+
  ifelse(KM_input1$Stage=="Stage I",35,
         ifelse(KM_input1$Stage=="Stage II",32,
                ifelse(KM_input1$Stage=="Stage III",47,65)))+
  (KM_input1$Age-26)

library(ggstatsplot)
ggscatterstats(
  data = KM_input1,
  x = points2,
  y = OS.time,
  bf.message = FALSE
)
ggsave(file=file.path(fig5.path,"Correlation points OS time.pdf"),width = 5,height = 5)

KM_input1$OS2<-as.factor(KM_input1$OS)
ggscatter(KM_input1, x = "points2", y = "OS.time",
          add = "reg.line",               # Add regression line
          conf.int = TRUE,                # Add confidence interval
          color = "OS2", palette = "Dark2", # Color by groups "cyl"
          #shape = "OS"                   # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = OS2), label.x = 180)       # Add correlation coefficient
ggsave(file=file.path(fig5.path,"Correlation points OS time subgroup.pdf"),width = 5,height = 5)


p <- ggboxplot(KM_input1, x = "OS", y = "points2",
               color = "OS", palette = "Dark2",
               add = "jitter")

p + stat_compare_means(method = "t.test")

ggsave(file=file.path(fig5.path,"points subgroup.pdf"),width = 5,height = 5)

##计算预测C-index 并绘图

library(dynpred)
cindex.sig<-cindex(Surv(OS.time,OS) ~Age+Stage+Tumor_status+signature, data = pbc)
cindex.age<-cindex(Surv(OS.time,OS) ~Age, data = pbc)
cindex.stage<-cindex(Surv(OS.time,OS) ~Stage, data = pbc)
cindex.grade<-cindex(Surv(OS.time,OS) ~Tumor_status, data = pbc)

df <- data.frame(dose=c(" Nomogram", "Age", "Stage","Tumor_status"),
                 len=c(round(cindex.sig$cindex,3), round(cindex.age$cindex,3),
                       round(cindex.stage$cindex,3),round(cindex.grade$cindex,3)))

# Change barplot fill colors by groups
p<-ggplot(data=df, aes(x=dose, y=len,fill=dose)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_bw()
p2<-p+scale_fill_brewer(palette="Dark2")+ theme(legend.position="none")

p2
ggsave(file=file.path(fig5.path,"C-index nomogram others.pdf"),width = 4,height = 4)

#######################Figure 6  

##################KIRC

library(MOVICS)
KIRCmulti$Risk<- factor(ifelse(KIRCmulti$signature>median(KIRCmulti$signature),"high"," low"))

KIRCmulti$group<-ifelse(KIRCmulti$Risk==" low","1","2")
head(KIRCmulti)
#write.table(KIRCmulti,file.path(fig6.path,"KIRCmulti.group.txt"),sep = "\t",col.names=NA,quote = F)

pseudo.moic.KIRC <- list("clust.res" = KIRCmulti[,c("Risk","group")],
                         "mo.method" = "test")
pseudo.moic.KIRC$clust.res$samID <- rownames(pseudo.moic.KIRC$clust.res)
pseudo.moic.KIRC$clust.res$clust <- pseudo.moic.KIRC$clust.res$group

library(MOVICS)

clin.KIRC <- compClinvar(moic.res      = pseudo.moic.KIRC,
                         var2comp      = KIRCmulti, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("OS","PFI"), # features that are considered categorical variables
                         nonnormalVars = "longest_dimension", # feature(s) that are considered using nonparametric test
                         #exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         res.path = fig6.path,
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")

print(clin.KIRC$compTab)


mut <- read.delim(file.path(data.path,"KIRC_mc3_gene_level.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
mut <- as.data.frame(na.omit(mut))
rownames(mut) <- mut$sample; mut <- mut[,-1]
mut <- mut[rowSums(mut) > 0,]
mut <- mut[,which(substr(colnames(mut), 14, 15) == "01")]
write.table(mut, "TCGA_KIRC_mut1.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

mut<-mut[,intersect(colnames(mut),pseudo.moic.KIRC$clust.res$samID)]
rowSums(mut)
intmut<-intersect(colnames(mut),pseudo.moic.KIRC$clust.res$samID)

pseudo.moic.KIRC2<-list("clust.res" = KIRCmulti[intmut,c("Risk","group")],
                        "mo.method" = "test")
pseudo.moic.KIRC2$clust.res$samID <- rownames(pseudo.moic.KIRC2$clust.res)
pseudo.moic.KIRC2$clust.res$clust <- pseudo.moic.KIRC2$clust.res$group


mut.KIRC <- compMut(moic.res     = pseudo.moic.KIRC2,
                    mut.matrix   = mut, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.04, # keep those genes that mutated in at least 5% of samples
                    p.cutoff = 0.1,
                    clust.col = jco,
                    p.adj.cutoff = 0.2, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    #innerclust   = T, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    #annColors    = annColors, # same annotation color for heatmap
                    fig.path      = fig6.path,
                    width        = 11, 
                    height       = 4,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

print(mut.KIRC)

maf <- read_tsv(file.path(data.path,"data_mutations_extended.txt"), comment = "#")

label <- c("Tumor_Sample_Barcode",
           "Hugo_Symbol",
           "Chromosome",
           "Start_Position",
           "End_Position",
           "Variant_Classification",
           "Variant_Type",
           "Reference_Allele",
           "Tumor_Seq_Allele1",
           "Tumor_Seq_Allele2")
maf <- maf[,label]

head(maf)

tmb.KIRC <- compTMB(moic.res     = pseudo.moic.KIRC,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

head(tmb.KIRC$TMB)


segment <- read.table(file.path(data.path,"KIRC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
colnames(segment) <- c("sample","chromosome","start","end","num.mark","seg.mean")
segment$sample <- substr(segment$sample, start = 1,stop = 15)
segment <- segment[which(substr(segment$sample, 14, 15) == "01"),]

segment$value<-segment$seg.mean
segment$chrom<-segment$chromosome
head(segment)

fga.KIRC <- compFGA(moic.res     = pseudo.moic.KIRC,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA",
                    width = 7,
                    height = 4)
head(fga.KIRC$summary)


GSET.FILE <-
  system.file("extdata", "28immunemarker.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.KIRC,
          norm.expr     = KIRC.expr,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          clust.col = jco,
          annColors     = annColors,
          color         =c("white","black", "#D0482E"),
          fig.path      = fig6.path,
          fig.name      = "t28immune GENE SETS OF INTEREST HEATMAP",
          height        = 6,
          width         = 12)

GSET.FILE <- 
  system.file("extdata", "h.all.v7.3.symbols.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.KIRC,
          norm.expr     = KIRC.expr[,pseudo.moic.KIRC$clust.res$samID],
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          clust.col = jco,
          annCol        = annCol,
          annColors     = annColors,
          color         =c("#00CBFD","black","#FDFD00" ),
          fig.path      = fig6.path,
          fig.name      = "HALLMARK GENE SETS OF INTEREST HEATMAP",
          height        = 8,
          width         = 12)

GSET.FILE <- 
  system.file("extdata", "Renal genesets.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = pseudo.moic.KIRC,
          norm.expr     = KIRC.expr[,pseudo.moic.KIRC$clust.res$samID],
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          clust.col = jco,
          annCol        = annCol,
          annColors     = annColors,
          color         =c("#25ABC0","#461254","#D34487" ),
          fig.path      = fig6.path,
          fig.name      = "Renal genesets GENE SETS OF INTEREST HEATMAP",
          height        = 8,
          width         = 12)


drug.KIRC.KIRC<- compDrugsen(moic.res    = pseudo.moic.KIRC,
                        norm.expr   = KIRC.expr[,pseudo.moic.KIRC$clust.res$samID], # double guarantee sample order
                        drugs       = c("Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine",
                                        "AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine",
                                        "Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84",
                                        "BI-2536","BMS-536924","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X",
                                        "MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal",
                                        "Lapatinib","GSK269962A","Doxorubicin","Etoposide","Gemcitabine","Mitomycin C","Vinorelbine",
                                        "NSC-87877","Bicalutamide","QS11","CP466722","Midostaurin","CHIR-99021","AP-24534","AZD6482",
                                        "JNK-9L","PF-562271","HG-6-64-1","JQ1","JQ12","DMOG","FTI-277","OSU-03012","Shikonin",
                                        "AKT inhibitor VIII","Embelin","FH535","PAC-1","IPA-3","GSK-650394","BAY 61-3606",
                                        "5-Fluorouracil","Thapsigargin","Obatoclax Mesylate","BMS-754807","Lisitinib",
                                        "Bexarotene","Bleomycin","LFM-A13","GW-2580","AUY922","Phenformin","Bryostatin 1"), # a vector of names of drug in GDSC
                        tissueType  = "all", # choose specific tissue type to construct ridge regression model
                        clust.col   = c(jco[2],jco[1]),
                        test.method = "parametric", # statistical testing method
                        seed = 123456,
                        fig.path = fig7.path,
                        prefix      = "KIRC IC50") 
#1,2,3,5,8,9,10,12,13,14,15,18,20,22,23,25-29,31-33,43,44,46,49,50,52-55,57,60,63-69,72-75,77-79,82
  
library(listr)

drugKIRC<-list_select(drug.KIRC.KIRC,1,2,3,5,8,9,10,12,13,14,15,18,20,22,23,25:29,31:33,43,44,46,49,50,52:55,57,60,63:69,72:75,77:79,82)
drugKIRC<-as.data.frame(drugKIRC)
drugKIRC$Risk<-drugKIRC$Erlotinib.Subtype
del <- seq(2, nrow(drugKIRC), by = 2)

KIRCdrug<-drugKIRC[,-del]

KIRCdrug<- KIRCdrug[order(KIRCdrug$Risk, decreasing = F), ]

KIRCdrugs<-data.frame(colMeans(KIRCdrug[which(KIRCdrug$Risk=="CS1"),1:49])-
  colMeans(KIRCdrug[which(KIRCdrug$Risk=="CS2"),1:49]))
colnames(KIRCdrugs)<-"x"
KIRCdrugs$ID<-rownames(KIRCdrugs)

rownames(KIRCdrugs[which(KIRCdrugs$x<0),])

library(ggpubr)
library(rstatix)
KIRCdrug1<-KIRCdrug[,c("Risk",rownames(KIRCdrugs[which(KIRCdrugs$x<0),]))]
#KIRCdrug1$clust<-paste0("C",KIRCdrug1$clust)
KIRCdrug1<-KIRCdrug1  %>% 
  as_tibble() %>%
  gather(Drugs, IC50,  -Risk)%>% 
  mutate(Drugs=factor(Drugs, levels=colnames(KIRCdrug1)[-1]))

KIRCdrug1$Drugs<-gsub(".Est.IC50","",KIRCdrug1$Drugs)
KIRCdrug1$Risk<-ifelse(KIRCdrug1$Risk=="CS1"," low-Risk","High-Risk")
KIRCdrug1<-as.data.frame(KIRCdrug1)
head(KIRCdrug1)

stat.test <- KIRCdrug1 %>%
  group_by(Drugs) %>%
  t_test(IC50 ~ Risk) %>%
  add_significance()
#stat.test$y.position<-stat.test$y.position+1
stat.test <- stat.test %>% add_xy_position(x = "Risk")

p <- ggboxplot(KIRCdrug1, x = "Risk", y = "IC50",
               color = "Risk",
               fill = "Risk",
               lwd=0.5, alpha = 0.3,
               palette = c("#169FB0","#F64C32","darkblue"),
               add.params = list(size = 0.5, jitter = 0.2), 
               add = "jitter",
               facet.by = "Drugs", 
               nrow = 6,
               ncol = 7,
               scales = "free",
               short.panel.labs = FALSE
)
p
p+ stat_pvalue_manual(
  stat.test, 
  bracket.nudge.y =-1, 
  label = "{p}{p.signif}"
) 

ggsave(file.path(fig6.path,"KIRCdrug low.pdf"),width = 18,height = 15)

library(ggpubr)
library(rstatix)
KIRCdrug2<-KIRCdrug[,c("Risk",rownames(KIRCdrugs[which(KIRCdrugs$x>0),]))]
#KIRCdrug2$clust<-paste0("C",KIRCdrug2$clust)
KIRCdrug2<-KIRCdrug2  %>% 
  as_tibble() %>%
  gather(Drugs, IC50,  -Risk)%>% 
  mutate(Drugs=factor(Drugs, levels=colnames(KIRCdrug2)[-1]))

KIRCdrug2$Drugs<-gsub(".Est.IC50","",KIRCdrug2$Drugs)
KIRCdrug2$Risk<-ifelse(KIRCdrug2$Risk=="CS1"," low-Risk","High-Risk")
KIRCdrug2<-as.data.frame(KIRCdrug2)
head(KIRCdrug2)

stat.test <- KIRCdrug2 %>%
  group_by(Drugs) %>%
  t_test(IC50 ~ Risk) %>%
  add_significance()
#stat.test$y.position<-stat.test$y.position+1
stat.test <- stat.test %>% add_xy_position(x = "Risk")

p <- ggboxplot(KIRCdrug2, x = "Risk", y = "IC50",
               color = "Risk",
               fill = "Risk",
               lwd=0.5, alpha = 0.3,
               palette = c("#5C0A98","#00CBFD","darkblue"),
               add.params = list(size = 0.5, jitter = 0.2), 
               add = "jitter",
               facet.by = "Drugs", 
               nrow = 2,
               ncol = 7,
               scales = "free",
               short.panel.labs = FALSE
)
p
p+ stat_pvalue_manual(
  stat.test, 
  bracket.nudge.y =-1, 
  label = "{p}{p.signif}"
) 

ggsave(file.path(fig6.path,"KIRCdrug High.pdf"),width = 18,height = 6)

#










###预测高低风险组差异激活的基因

subt <- data.frame(condition =KIRCmulti$Risk,
                   row.names = rownames(KIRCmulti))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  =KIRC.expr[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "high-risk", # name of treatment group
              ctrlVar  = " low-risk", # name of control group
              prefix   = "KIRCmulti", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = fig6.path ) # path for result

deg.KIRC<- read.table(file.path(fig6.path,"KIRCmulti_limma_test_result.high-risk_vs_ low-risk.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

library(fgsea)
deg.KIRC$log2fc <- as.numeric(as.character(deg.KIRC$log2fc))
deg.KIRC <- deg.KIRC[order(deg.KIRC$log2fc, decreasing = T), ]
si.id <- deg.KIRC$log2fc; names(si.id) <- rownames(deg.KIRC)
head(si.id)
gmtfile <- "h.all.v7.0.symbols.gmt"
hallmark <- read.gmt(file.path(data.path,gmtfile))
hallmark$term<-gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

fgseaRes <- fgsea(pathways = hallmark.list, 
                  stats = si.id,
                  minSize=5,
                  maxSize=500,
                  nperm=10000)
sig<-fgseaRes[fgseaRes$padj<0.05,]
sig<-sig[order(sig$NES,decreasing = T)]

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
length(topPathways)

pdf(file.path(fig6.path,"KIRCmulti fgsea.pdf"), width = 8, height = 6)

plotGseaTable(hallmark.list[topPathways], si.id, fgseaRes, 
              gseaParam = 0.5)
dev.off()


################## KIRC SubMap 免疫治疗预测结果。根据结果调整，不好可以不放


dat <- KIRC.expr[,pseudo.moic.KIRC$clust.res$samID]
dat[1:3, 1:3]

ann <- pseudo.moic.KIRC$clust.res[,c(21:22)]
rownames(ann)<-ann$samID

colnames(ann)<-c("ID","ImmClust")
ann$ImmClust<-paste0("C",ann$ImmClust)
head(ann)

table(ann$ImmClust)

TIDE <- log2(dat+1)
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
write.table(TIDE,file.path(fig6.path,"TIDE_input.self_subtract"),sep = "\t",row.names = T,col.names = NA,quote = F)

#------------------------------------#
# 
#------------------------------------#

# 
TIDE.res <- read.csv(file.path(data.path,"TIDE_output.csv"),header = T,row.names = 1,check.names = F,stringsAsFactors = F)
ann$TIDE <- TIDE.res[rownames(ann),"Responder"]
print(table(ann$TIDE,ann$ImmClust))

print(fisher.test(table(ann$TIDE,ann$ImmClust))) 

# 
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 
skcm.immunotherapy.logNC <- read.table(file.path(data.path,"skcm.immunotherapy.47samples.log2CountsNorm.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原???峁???log2转???谋?准??count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # ??????写????为??使?玫??????前鸦?????????写??
skcm.immunotherapy.info <- read.table(file.path(data.path,"skcm.immunotherapy.47sampleInfo.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# 
tmp <- dat
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

#
samples.C1 <- rownames(ann[which(ann$ImmClust == "C1"),])
samples.C2 <- rownames(ann[which(ann$ImmClust == "C2"),])
samples.C3 <- rownames(ann[which(ann$ImmClust == "C3"),])

sam_info <- data.frame("ImmClust"=c(samples.C1,samples.C2,samples.C3),row.names = c(samples.C1,samples.C2,samples.C3))
sam_info$rank <- rep(c(1,2,3),times=c(length(samples.C1),length(samples.C2),length(samples.C3))) #1: C1,??HPV16-IMM 2: C2,??HPV16-KRT

# 
gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- tmp[GENELIST,rownames(sam_info)] # ??????示?????????频???式??log2转???谋?准??count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")


#-------------------------------

heatmap.YlGnPe <- c("#9A1424","#DE518D","#FF9F1C","#B4D834","#2EC4B6")
cherry    <- "#700353"
lightgrey <- "#dcddde"

tmp <- matrix(c(  0.1818182, 0.4155844, 0.5824176, 0.92007992,
                  0.9360639, 0.4255744, 0.8961039, 0.01198801,
                  1,  1,  1, 1.0000000,
                  1,  1,  1, 0.0959041
), # Bonferroni
nrow = 4,byrow = T,dimnames = list(c("Low-Risk","High-Risk","Low-Risk-b","High-Risk-b"),c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))
library(pheatmap)
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[5:1],
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "KIRC heatmap_submap.pdf")


#####table

library(table1) 
library(boot)
library(flextable)
library(magrittr)

KIRCtable<-KIRC.clin[,c("OS.time","OS","Age","Stage","Grade")]
#colnames(KIRCtable)<-c("OS.time","OS","Gender","Age","Metastatic","Tumor site")
KIRCtable$subgroup<-"TCGA-KIRC"

EMTAB3267.clin$Age<-"NA"
EMTAB3267.clin$Grade<-""
EMTAB3267table<-EMTAB3267.clin[,c("OS.time","OS","Age","Stage","Grade")]
EMTAB3267table$subgroup<-"EMTAB3267"

GSE22541.clin$Grade<-""
GSE22541table<-GSE22541.clin[,c("OS.time","OS","Age","Stage","Grade")]
GSE22541table$subgroup<-"GSE22541"


data <- rbind(KIRCtable,EMTAB3267table,GSE22541table)
colnames(data)<-c("Survival time","Status","Age","Stage","Grade","subgroup")
#data<-data[,c(3:8,16:17,25:25)]
data$Status<-ifelse(data$Status=="0","No","Yes")
head(data)
data$Age<-as.numeric(data$Age)
data[data==""]<-"unknow"

pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    p <- t.test(y ~ g)$p.value #数值型数据用t-test(两组比较)
  } else {
    p <- fisher.test(table(y, g))$p.value #因子型数据用卡方
  }
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

tbl1<-table1(~  `Survival time`+Status+Age+Stage+Grade| subgroup, 
             data=data)
tbl1

# Convert to flextable
t1flex(tbl1) %>% 
  save_as_docx(path="Table 1.docx")


geneids


#---------------TCGA compared other signatures
#--------------------------#
# HMAG et.al # # complete #
#--------------------------#

tcga.gene <- rownames(KIRC.expr)

HMAG.sig <- c("TNFRSF12A","MYO10","PDK1","FASN","MYLIP","FGFR2","HLF","SOX21","ZDHHC11")
HMAG.sig <- intersect(HMAG.sig,tcga.gene)
KIRC.expr.sinfo <- KIRC.clin[colnames(KIRC.expr),]
tmp <- as.data.frame(t(KIRC.expr[HMAG.sig,]))
tmp$OS <- KIRC.expr.sinfo[rownames(tmp),"OS"]
tmp$OS.time <- KIRC.expr.sinfo[rownames(tmp),"OS.time"]

HMAG.coef <- c(0.115713151, 0.032286676, 0.300717617, 0.131853347, -0.015900246, -0.033923872, -0.099329649, -0.034241118, -0.042765486)

#HMAG.coef <- as.numeric(cox.HMAG$coefficients[,1])
HMAG.risk <- apply(KIRC.expr[HMAG.sig,],2,function(x){x %*% HMAG.coef})

tcgacomp<-tmp[,10:11]
tcgacomp$HMAG<-HMAG.risk


KIRC.expr<-read.table(file.path(data.path,"KIRC.tumor.txt"),header=T,sep="\t",row.names=1,check.names=F)
merge<-intersect(colnames(KIRC.expr),rownames(KIRC.clin))
KIRC.expr<-KIRC.expr[,merge]
KIRC.clin<-KIRC.clin[merge,]
KIRC.expr<-log2(KIRC.expr+1)
range(KIRC.expr)
#------------------------#
#Yu-relatedet.al # # complete #
#------------------------#


Yu.sig <- c("CHIT1","GNG8","GTSF1L","PLA2G2D")
Yu.sig <- intersect(Yu.sig,rownames(KIRC.expr))

tmp <- as.data.frame(t(KIRC.expr[Yu.sig,]))
tmp$OS <- KIRC.expr.sinfo[rownames(tmp),"OS"]
tmp$OS.time <- KIRC.expr.sinfo[rownames(tmp),"OS.time"]

Yu.coef <- c(0.7071341,0.6469606,0.5702634,0.7482603)
Yu.risk <- apply(log(KIRC.expr[Yu.sig,]+1),2,function(x){x %*% Yu.coef})

tcgacomp$Yu<-Yu.risk


#------------------------#
# Zhou et.al # # complete #
#------------------------#


Zhou.sig <- c("ACTR2","TEX12","UBE2V1","HSF1","FBXO6")
Zhou.sig <- intersect(Zhou.sig,rownames(KIRC.expr))

Zhou.coef <- c(0.418,-1.995,0.147,0.543 ,-0.217)
Zhou.risk <- apply(log(KIRC.expr[Zhou.sig,]+1),2,function(x){x %*% Zhou.coef})

tcgacomp$Zhou<-Zhou.risk

#------------------------#
# Qi et.al # # complete #
#------------------------#
Qi.sig <- c("TFRC","ACACA","SQLE","PHKG2")
Qi.sig <- intersect(Qi.sig,rownames(KIRC.expr))

Qi.coef <- c(0.195,0.104,0.097,-0.512)
Qi.risk <- apply(log(KIRC.expr[Qi.sig,]+1),2,function(x){x %*% Qi.coef})

tcgacomp$Qi<-Qi.risk

#------------------------#
# Chen et.al # # complete #
#------------------------#
Chen.sig <- c("ATG3","BCL2","CD46","IFNG","NAMPT","TM9SF1")
Chen.sig <- intersect(Chen.sig,rownames(KIRC.expr))

Chen.coef <- c(-0.63,-0.42,0.85,-0.38,0.23,0.82)
Chen.risk <- apply(log(KIRC.expr[Chen.sig,]+1),2,function(x){x %*% Chen.coef})

tcgacomp$Chen<-Chen.risk



#--------------TCGA time ROC

library("pROC")
yjp <- c("#5C0A98","#EA4440","#9A8266","#FFA900","#14B96A")
Yu.roc <- roc(tcgacomp$OS, tcgacomp$Yu, smooth = F);Yu.roc$auc
HMAG.roc <- roc(tcgacomp$OS, tcgacomp$HMAG, smooth = F);HMAG.roc$auc
Zhou.roc <- roc(tcgacomp$OS, tcgacomp$Zhou, smooth = F);Zhou.roc$auc
Qi.roc <- roc(tcgacomp$OS, tcgacomp$Qi, smooth = F);Qi.roc$auc
Chen.roc <- roc(tcgacomp$OS, tcgacomp$Chen, smooth = F);Chen.roc$auc

pdf("TCGA-ROC.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25,las = 1)
# plot the first roc
plot(1-Yu.roc$specificities, Yu.roc$sensitivities,type="l", xlim=c(0,1), ylim=c(0,1),col=yjp[1],
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "")
# add diagonal line
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.22,0.24,"Yu AUC: 0.632",cex = 1,col = yjp[1],adj = 0)
# add other rocs
lines(1-HMAG.roc$specificities, HMAG.roc$sensitivities, col=yjp[2],lwd = 2)
text(0.22,0.18,"HMAG AUC: 0.711",cex = 1,col = yjp[2],adj = 0)
# add other rocs
lines(1-Zhou.roc$specificities, Zhou.roc$sensitivities, col=yjp[3],lwd = 2)
text(0.22,0.12,"Zhou AUC: 0.572",cex = 1,col = yjp[3],adj = 0)
#
lines(1-Qi.roc$specificities, Qi.roc$sensitivities, col=yjp[4],lwd = 2)
text(0.22,0.06,"Qi AUC: 0.652",cex = 1,col = yjp[4],adj = 0)
#
lines(1-Chen.roc$specificities, Chen.roc$sensitivities, col=yjp[5],lwd = 2)
text(0.22,0.00,"Chen AUC: 0.661",cex = 1,col = yjp[5],adj = 0)

invisible(dev.off())

###################比较signaure与其他的预测C-index
library(dynpred)


cindex.sig<-cindex(Surv(OS.time,OS) ~HMAG, data = tcgacomp)
cindex.sig1<-cindex(Surv(OS.time,OS) ~Yu, data = tcgacomp)
cindex.sig2<-cindex(Surv(OS.time,OS) ~Zhou, data = tcgacomp)
cindex.sig3<-cindex(Surv(OS.time,OS) ~Qi, data = tcgacomp)
cindex.sig4<-cindex(Surv(OS.time,OS) ~Chen, data = tcgacomp)

df <- data.frame(dose=c("AHMAG","Yu","Zhou","Qi","Chen"),
                 len=c(round(cindex.sig$cindex,3), 
                       round(cindex.sig1$cindex,3), 
                       round(cindex.sig2$cindex,3), 
                       round(cindex.sig3$cindex,3), 
                       round(cindex.sig4$cindex,3)))
df

# Change barplot fill colors by groups
p<-ggplot(data=df, aes(x=dose, y=len,fill=dose)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_bw()
p2<-p+scale_fill_brewer(palette="Paired")+ theme(legend.position="none")

p2
ggsave("C-index nomogram  TCGA.pdf",width = 4,height = 4)





ggplot(data=aosi_data, aes(x=GROUP, y=V12.aosi.total_score_1_18, color=Gender))+
  geom_boxplot()

wnt<-intersect(c("WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","WNT10A","WNT10B","WNT11","WNT16"),
               rownames(KIRC.expr))
intmut

WNTS<-data.frame(cbind(t(mut["PBRM1",intmut]),t(KIRC.expr[wnt,intmut])))


table(WNTS$PBRM1)
WNTS<-WNTS[order(WNTS$PBRM1,decreasing = F),] 

conNum=219                           
treatNum=142 

rt=WNTS[,c(2:2,5:5,7:7,9:9,12:12,14:14)]

type=c(rep(" Wild",conNum),rep("Mut",treatNum))

data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=rt[,i],gene=i,type))
}
write.table(data,file.path(data.path,"data.txt"),sep="\t",row.names=F,quote=F)

data=read.table(file.path(data.path,"data.txt"),sep="\t",header=T,check.names=F)       #??取????图?????募?
p=ggboxplot(data, x="gene", y="expression", notch = T,
            lwd=0.5, alpha = 0.6,color = "type", fill = 'type',
            #width=1,
            ylab="Gene expression, log2(TPM+1)",
            xlab="",
            palette = c( "grey","#B41947","#5A58A6","#ED2024") )

p=p+rotate_x_text(45)
pdf(file=file.path(fig6.path,"boxplot for Wnts PBRM1 mut expr.pdf"),width=5,height=4)                          #????图片?募?  #t.test
p+stat_compare_means(aes(group=type),method = "anova", symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),label = "p.signif")
dev.off()


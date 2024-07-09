#### 0. What is the GWAS signal density plot?

**GWAS signal density plot is a heatmap plot, summarizing the number of signals (p < 0.05) per 10-Mbp window on the genome.** 

After performing GWAS analysis, a GWAS signal density plot can help illustrate the genetic architecture of the targeted phenotype and conduct comparisons between traits. For instance, in one recent study investigating distinct biological ages of organs, GWAS signal density plots were used to demonstrate the genetic distinctness between different organ ages (Nie et al., 2022) [^1]. 

![image-20240709145905926](https://github.com/EnricMartin/GWAS-signal-density-plot/blob/main/image-20240709145905926.png)

In other studies, a similar idea is also expressed by the Brisbane plot [^2].

![image-20240709150740445](https://github.com/EnricMartin/GWAS-signal-density-plot/blob/main/image-20240709150740445.png)

#### 1. Why you need this tutorial

The reason I wrote this tutorial is very simple: there isn't an open-access code for drawing such a plot, though it is actually very simple. There indeed is an R package "CMplot" for painting the Manhattan plot and the SNP-density plot [^3]. However, using "CMplot" cannot put all 22 chromosomes into one row, which would take much place and lack clarity. 

#### 2. Preparations

**2.1 R packages preparation**

```R
library(pheatmap)
library(vegan)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grDevices)
library(data.table)
library(dplyr)
```

**2.2 Data preparation**

GWAS summary data contains 4 columns: CHR, BP, P-value.

#### 3. Code

**3.1 The data we used**

SNP_data: CHR, BP, P-value.

**3.2 Set the Heatmap color**

```R
col=c("darkgreen", "yellow", "orange","red")
col=colorRampPalette(col)(9)
col_fun=colorRamp2(c(0,100,200,300,400,500,600,700,800), col) ##The heatmap color depends on the unit of SNP density [the number of signals (p < 0.05) per 10-Mbp window, if it exceeds 1000, than the code above should update accordingly].  
```

**3.3 Preprocess the data for Heatmap**

The core idea is to create a dataframe contains each interval of each chromosome and the number of SNPs within it.

```R
signif<-SNP_data[SNP_data$P_value<0.05,] ##Restrict to SNP with P<0.05.
max.chr <- max(SNP_data$CHR)
chr.num <- unique(SNP_data$CHR) ##Number of CHR. 
chorm.maxlen <- max(SNP_data$BP) ##The max basepair of this data.
bin=1e+7 ##Define the block size.
pos.x <- list() ##Create a list to store all the basepair of Chromosomes. 
chr.pos.max.v <- NULL ##Create a vector to store all the max BP of each Chromosome.
maxbin.num <- NULL
windinfo <- list()
signal<-NULL
for(i in 1 : length(chr.num)){
  pos.x[[i]] <- SNP_data$BP[SNP_data$CHR == chr.num[i]]
  maxposindx <- which.max(pos.x[[i]])
  max.pos <- pos.x[[i]][maxposindx]
  chr.pos.max.v <- c(chr.pos.max.v, max.pos)
  cut.breaks <- seq(0, max.pos, bin) ## calculate the break value for each chromosome
  cut.len <- length(cut.breaks)
  if(cut.breaks[length(cut.breaks)] < max.pos){  
    cut.breaks <- c(cut.breaks, cut.breaks[length(cut.breaks)] + bin)
  } ##If the basepair/binsize is not a integar, the number of blocks should be added one.
  cut.r <- cut(signif$BP[signif$CHR == chr.num[i]], cut.breaks, labels=FALSE)
  eachbin.num <- table(cut.r)
  eachbin.num<-as.data.frame(eachbin.num)
  match_bin<-data.frame(cut.r=c(1:length(cut.breaks)))
  match_bin<-merge(match_bin,eachbin.num,all.x = T)
  match_bin$Freq[is.na(match_bin$Freq)]<-0
  
  match_bin$chr<-rep(paste0(i),nrow(match_bin))
  signal<-rbind(signal,match_bin)
  maxbin.num <- c(maxbin.num, max(match_bin$Freq))
}
my_data<-t(signal[,2])
```

**3.4 Prepare annotations for Heatmap plot**

```R
##column split
list_or<-as.data.frame(signal[,3])
colnames(list_or)<-"Chromosome"
list_or$Chromosome<-as.numeric(list_or$Chromosome)

##bar plot
at_1<-round(max(my_data)/100)*100
ha <- HeatmapAnnotation(
  `SNP per 10 Mbp` = anno_barplot(
    signal$Freq,
    bar_width = 1, height = unit(0.5,"cm"),
    gp=gpar(fill="#9BB0C1",col="#9BB0C1"),
    axis_param = list(at=c(0,at_1*0.5,at_1))
  )
)
```

**3.5 Draw the Heatmap**

```R
a<-Heatmap(my_data,col = col_fun,
        cluster_columns = F,cluster_rows = F,
        width = unit(25,"cm"),height = unit(0.55,"cm"),
        show_heatmap_legend = F,
        column_split = list_or,column_gap = unit(2, "mm"),
        top_annotation = ha,
        show_column_names = FALSE,column_title_side = "top",
        show_row_names = T,column_title = " ",
        name = " "
)
```

#### 4. Outcome

![image-20240709175149367](https://github.com/EnricMartin/GWAS-signal-density-plot/blob/main/image-20240709175149367.png)

[^1]: Nie C, Li Y, Li R, et al. Distinct biological ages of organs and systems identified from a multi-omics study. *Cell Rep*. 2022;38(10):110459. doi:10.1016/j.celrep.2022.110459
[^2]: Yengo, L., Vedantam, S., Marouli, E. *et al.* A saturated map of common genetic variants associated with human height. *Nature* **610**, 704–712 (2022). https://doi.org/10.1038/s41586-022-05275-y
[^3]: [GitHub - YinLiLin/CMplot: 📊 Circular and Rectangular Manhattan Plot](https://github.com/YinLiLin/CMplot)



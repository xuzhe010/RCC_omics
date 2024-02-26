library(GenomicRanges)
library(tidyverse)
## read miRNA from miRNAbaseV22
allmirna<-read_tsv("KidneyCancer/CircleSeq/results.count/miRNA_primary_transcript.bed",col_names = F)
names(allmirna)<-c("seqnames","start","end","strand","id","name")
mirna<-makeGRangesFromDataFrame(allmirna,keep.extra.columns = T)
rm(allmirna)
find_intact<-function(x){
  print(x)
  ## readin ecc
  df<-read_tsv(stringr::str_glue("KidneyCancer/CircleSeq/raw.count/miRna.bed/cKca_{x}_miRNA.bed"),col_names = F)
  names(df)<-c("seqnames","start","end","eccDNA")
  gr<-makeGRangesFromDataFrame(df,keep.extra.columns = T)
  rm(df)
  ## find overlap with intact mirna
  hits <- findOverlaps(gr,mirna,ignore.strand=T)
  rst<-gr[queryHits(hits)]
  rst$mirna<-mirna[subjectHits(hits)]$id
  ints<-pintersect(gr[queryHits(hits)],mirna[subjectHits(hits)])
  rst$intW=width(ints)
  rst$geneC <- mirna[subjectHits(hits)]
  intact<-rst[width(rst$geneC)==rst$intW]
  if(length(intact)>0){
    return(tibble(sample=x,mirna=unique(intact$mirna)))
  }else{
    return(tibble(sample="NA",mirna="NA"))
  }
}

infos<-openxlsx::read.xlsx("KidneyCancer/CircleSeq/results.count/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2
cases<-paste0(need_samples,"T")
normals<-paste0(need_samples,"N")

fin<-do.call('rbind',lapply(c(cases,normals),find_intact))
fin<-fin%>%
  filter(sample!="NA")
fin<-fin%>%group_by(sample)%>%distinct(mirna,.keep_all = T)
fin$num<-1
## expand into matrix

miRNA_mat<-fin%>%tidyr::pivot_wider(names_from ="sample" ,values_from ="num")
miRNA_mat<-miRNA_mat%>%
  tibble::column_to_rownames(var="mirna")%>%
  as.matrix()
miRNA_mat[is.na(miRNA_mat)]<-0
left<-miRNA_mat[,1:68]
right<-miRNA_mat[,69:243]
middle<-matrix(rep(0,nrow(left)),ncol = 1)
colnames(middle)<-"82T"
rownames(middle)<-rownames(left)
fin<-cbind(left,middle,right)

### calculate numbers
case_mat<-fin[,1:122]
normal_mat<-fin[,123:244]

case_numbers<-apply(case_mat,1,function(x){sum(x==1)})
control_numbers<-apply(normal_mat,1,function(x){sum(x==1)})

stat_table<-tibble(miRNA=rownames(fin),Casein=case_numbers,Casenon=122-case_numbers,Normalin=control_numbers,Normalnon=122-control_numbers)

### fisher test

ha<-stat_table%>%
  rowwise()%>%
  mutate(pvalues=fisher.test(matrix(c(Casein,Casenon,Normalin,Normalnon),
                                    nrow = 2,
                                    byrow = T),
                             alternative = "two.sided")$p)
ha<-ha%>%arrange(pvalues,desc(Casein))

topGenes<-ha[1:30,]$miRNA

plotData<-fin[topGenes,]
lieshu<-apply(plotData,2,function(x) sum(x,na.rm=T))

ht<-Heatmap(plotData,
        show_column_names = F,
        cluster_columns = T,
        cluster_rows = T,
        show_column_dend = F,
        show_row_dend = F,
        show_heatmap_legend = F,
        rect_gp = gpar(type="none"),
        column_split = c(c(rep("RCC",122),rep("NAT",122))),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "white", fill = "#CCCCCC",lwd=1.4))
          if(plotData[i,j]==1){
            grid.rect(x = x, y = y, width = width, height = height/2, 
                      gp = gpar(fill = "#2E318C",lwd=0,col="#CCCCCC"))
          }
        },
        top_annotation = HeatmapAnnotation(genes=anno_barplot(lieshu,
                                                              border=F,
                                                              axis=T,
                                                              gp=gpar(fill="#2E318C",
                                                                      col="transparent"))),
        bottom_annotation = HeatmapAnnotation(Grade=infos$Grade,
                                              Stage = infos$stage,
                                              Histology=infos$Histology,
                                              col = list(Histology=c("ccRCC"="#e41a1c","CDC"="#377eb8","chRCC"="#4daf4a","pRCC"="#984ea3","sRCC"="#ff7f00"),
                                                         Grade=c("0"="#1b9e77","1"="#d95f02","2"="#7570b3","3"="#e7298a","4"="#a6761d"),
                                                         Stage=c("1"="#386cb0","2"="#beaed4","3"="#fdc086","4"="#f0027f"))),
        row_names_gp = gpar(fontface="italic",cex=0.7)
        
)

ht
pdf("ecc-miRNA-heatmap.pdf",width = 7.42,height = 5.24)
draw(ht)
dev.off()

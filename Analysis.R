####====Packages loading====####
update.packages()
list.of.packages <- c("locfit","hwriter","coda","jsonlite","askpass","pracma","stringi",
                      "plotrix","pkgconfig","memoise","cli","ashr","numDeriv","urltools","ggnewscale",
                      "gplots", "plyr","scales","zeallot","RMariaDB","ggplot2","nVennR","ggprism",
                      "RColorBrewer","data.table", "RSpectra","MatrixCorrelation","resample","grid",
                      "pheatmap","xml2", "devtools","resample", "foreach","wordcloud","stopwords",
                      "iterators","tcltk","parallel","doParallel", "stringr","doSNOW", "pals", "rvest",
                      "dplyr","readr","factoextra","cluster","pvclust","openxlsx","Hmisc","combinat","gmp")
list.of.packages.Bio <- c("IRanges","DOSE","clusterProfiler","geneplotter","genefilter","SummarizedExperiment","Rsamtools",
                          "Biostrings","AnnotationDbi","GenomicRanges","GenomeInfoDb","S4Vectors","BiocGenerics",
                          "doSNOW","mygene","apeglm","RUVSeq","EnhancedVolcano","GO.db", "goseq","simplifyEnrichment",
                          "biomaRt", "tximport","rtracklayer","geneLenDataBase", "ensembldb","gage",
                          "gageData","pathview","org.Mm.eg.db","DESeq2", "edgeR", "limma","qvalue")
#Packages that does not install yet
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
new.packages.Bio <- list.of.packages.Bio[!(list.of.packages.Bio %in% installed.packages()[, "Package"])]
#install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
if(length(new.packages.Bio)) {
  BiocManager::install(new.packages.Bio)
  devtools::install_github("grimbough/biomaRt", force = TRUE)
}
if(length(new.packages)) {install.packages(new.packages)}
#Loading all packages
lapply(c(list.of.packages,list.of.packages.Bio), library, character.only = TRUE)
####====Functions loading====####
pseudoLog10 <- function(x) { asinh(x/2)/log(10) }
OCcutOff<-function(geneSets, X, value=0.25){
  temp2<-c()
  for (i in 1:(nrow(temp)-1)) {
    A<-length(geneSets[[X$ID[i]]])
    temp3<-0
    temp4<-0
    for (k in (i+1):nrow(X)) {
      B<-length(geneSets[[X$ID[k]]])
      OC<-length(intersect(geneSets[[X$ID[i]]],geneSets[[X$ID[k]]]))/#overlap coefficient
        min(A,B)
      if(OC>=value){
        if(!X$Count[i]==X$Count[k]){
          if(X$Count[i]>X$Count[k]){
            next()
          }else{
            temp4<-1
          }
        }else if(!X$p.adjust[i]==X$p.adjust[k]){
          if(X$p.adjust[i]>X$p.adjust[k]){
            next()
          }else{
            temp4<-1
          }
        }else{
          next()
        }
      }
      if(temp4==1){next()}
    }
    if(temp3==0&temp4==0){
      if(sum(temp2%in%X$ID[i])==0){
        temp2<-c(temp2,X$ID[i])
      }
    }
    if(temp4==1){
      if(sum(temp2%in%X$ID[k])==0){
        temp2<-c(temp2,X$ID[k])
      }
    }
  }
  X<-X[X$ID%in%temp2,]
  return(X)
}
GOMostDescend<-function(x){
  Par_GOCC <- as.list(GOCCANCESTOR)
  Par_GOBP <- as.list(GOBPANCESTOR)
  Par_GOMF <- as.list(GOMFANCESTOR)
  GOlist<-c(Par_GOCC, Par_GOBP, Par_GOMF)
  for (i in 1:length(x)) {
    x<-x[!(x %in% GOlist[[x[i]]])]
  }
  return(x)
}
GOMostAscend<-function(x){
  Par_GOCC <- as.list(GOCCOFFSPRING)
  Par_GOBP <- as.list(GOBPOFFSPRING)
  Par_GOMF <- as.list(GOMFOFFSPRING)
  GOlist<-c(Par_GOCC, Par_GOBP, Par_GOMF)
  for (i in 1:length(x)) {
    x<-x[!(x %in% GOlist[[x[i]]])]
  }
  return(x)
}
####====|01.Get Runx2 associated protein list from STRING Database|====####
download.file(
  read_html('https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus')%>%
    html_nodes(".download_file")%>%html_attr("href")%>%.[1],
  destfile = "protein.links.txt.gz"
)
protein.links<-read.table(gzfile("protein.links.txt.gz"),header = TRUE)
protein.links<-protein.links[protein.links$combined_score>=500,]#combined_score>=500
Runx2<-"10090.ENSMUSP00000109202"
protein.links.protein2_Runx2<-protein.links[protein.links$protein2==Runx2,]
#Get annotation of the retrieved proteins from UniPort
temp<-data.frame()
for (i in 1:length(protein.links.protein2_Runx2$protein1)) {
  temp<-rbind(
    temp,
    cbind(
      protein.links.protein2_Runx2$protein1[i],
      fread(
        paste0(
          "https://www.uniprot.org/uniprot/?query=",
          protein.links.protein2_Runx2$protein1[i],
          "+organism:10090&format=tab&limit=1&columns=id,genes,go,database(KEGG),keywords,mass"
        ),verbose=FALSE,showProgress=FALSE
      )
    )
  )
  cat(i, paste0("of ",length(protein.links.protein2_Runx2$protein1),"\r"))
  flush.console()
}
UniportAnno<-temp
UniportAnno$Mass<-as.numeric(gsub(",", "", UniportAnno$Mass))
UniportAnno<-rbind(
  UniportAnno[UniportAnno$Mass>=105000&UniportAnno$Mass<=186000,],
  UniportAnno[UniportAnno$Mass>=66000&UniportAnno$Mass<=86000,],
  UniportAnno[UniportAnno$Mass>=21000&UniportAnno$Mass<=45000,]
)
UniportAnno<-UniportAnno[order(UniportAnno$Mass, decreasing=TRUE),]
UniportAnno$`Gene names`<-unlist(lapply(strsplit(UniportAnno$`Gene names`," "), `[[`, 1))
colnames(UniportAnno)[3]<-"GeneSymbol"
UniportAnno$V1<-unlist(lapply(strsplit(UniportAnno$V1,"\\."), `[[`, 2))
colnames(UniportAnno)[1]<-"ensembl_peptide_id"
#|Set up Biomart parameters-
#View(listMarts()) #Check database list
Ann_0mart <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
#View(listDatasets(mart=Ann_0mart)) #Check speices list
Ann_0mart <- useDataset("mmusculus_gene_ensembl", Ann_0mart)
#|Get ENSEMBL Gene ID
ensembl_gene_id <- getBM(
  mart=Ann_0mart,
  attributes=c("ensembl_gene_id","ensembl_peptide_id"),
  filter="ensembl_peptide_id",
  values=as.vector(UniportAnno$ensembl_peptide_id),
  uniqueRows = TRUE)
UniportAnno<-merge(ensembl_gene_id,UniportAnno,'ensembl_peptide_id')
write.csv(UniportAnno,file = "Runx2Associated.csv", row.names = FALSE)
####|literature-mining|####
N<-paste0(#To search for the total number of PubMed citations, enter all[sb] in the search box.(20210506)
  "https://pubmed.ncbi.nlm.nih.gov/?term=",
  URLencode('all[sb]'),
  "&sort=date"
)%>%read_html()%>%html_nodes(".value")%>%html_text()%>%.[1]%>%gsub(",","",.)%>%as.numeric()
K<-paste0(#the amounts of literature related to osteoblast differentiation #20210506# PubMed Query: 
  "https://pubmed.ncbi.nlm.nih.gov/?term=",
  URLencode('("Osteoblasts"[Mesh]) AND "Cell Differentiation"[Mesh]'),
  "&sort=date"
)%>%read_html()%>%html_nodes(".value")%>%html_text()%>%.[1]%>%gsub(",","",.)%>%as.numeric()
Thresh<-K/N*10
#Loading Mesh term for the genes with recorded O-GlcNAcylation according to O-GlcNAcAtlas (https://oglcnac.org/)
Mesh_term4Genes<-read.csv("Mesh_term4Genes.csv")
UniportAnno<-merge(UniportAnno, Mesh_term4Genes, by="GeneSymbol")
#Calculate results and construct the literature mining table
literature.mining<-data.frame(matrix(nrow = 0, ncol=14))
colnames(literature.mining)<-c("GeneSymbol","ensembl_peptide_id","ensembl_gene_id",
                               "Entry", "n", "k", "R", "10x Background of R","P.value","Gene ontology (GO)",
                               "Cross-reference (KEGG)","Keywords","Mass","Mesh Term")
for (j in 1:nrow(Mesh_term4Genes)) {
  n<-paste0(#the amounts of papers for each gene #20210506# PubMed Query:
    "https://pubmed.ncbi.nlm.nih.gov/?term=",
    URLencode(UniportAnno$Mesh_Term[j]),
    "&sort=date"
  )%>%read_html()%>%html_nodes(".value")%>%html_text()%>%.[1]%>%gsub(",","",.)%>%as.numeric()
  temp<-paste0(#the amounts of papers of each gene on osteoblast differentiation #20210506# PubMed Query:
    "https://pubmed.ncbi.nlm.nih.gov/?term=",
    URLencode(
      paste0(
        "(",UniportAnno$Mesh_Term[j],")",
        " AND ",
        '("Osteoblasts"[Mesh]) AND "Cell Differentiation"[Mesh]'
      )
    ),
    "&sort=date"
  )%>%read_html()
  if(identical(
    html_nodes(temp,".single-result-redirect-message")%>%html_text(),
    character(0)
  )){
    k<-html_nodes(temp,".value")%>%html_text()%>%.[1]%>%gsub(",","",.)%>%as.numeric()
  }else{k<-1}
  if(is.na(k)){k<-0}
  R<-k/n
  P<-0
  for (i in 0:(k-1)) {
    P<-P+((chooseZ(K,i)*chooseZ(N-K,n-i))/chooseZ(N,n))
  }
  P<-asNumeric(1-P)
  if(P==0){P<-1*10^-300}
  temp<-as.data.frame(
    c(UniportAnno[j,1:4],n,k,R,Thresh,P,
      UniportAnno[j,5:9])
  )
  colnames(temp)<-c("GeneSymbol","ensembl_peptide_id","ensembl_gene_id",
                    "Entry", "n", "k", "R", "10x Background of R","P.value","Gene ontology (GO)",
                    "Cross-reference (KEGG)","Keywords","Mass","Mesh Term")
  literature.mining<-rbind(literature.mining,temp)
}
#Correcting P-value for multiple testing
literature.mining<-cbind(literature.mining[1:9], 
                         p.adj=p.adjust(as.numeric(literature.mining$P.value),method = "BH"),
                         literature.mining[10:14])
#Plotting Results
literature.mining$n<-as.numeric(literature.mining$n)
literature.mining$k<-as.numeric(literature.mining$k)
literature.mining$R<-as.numeric(literature.mining$R)
df<-literature.mining  
df$R<-df$R*1000
df$n<-log10(df$n)
df$R<-pseudoLog10(df$R)
df$p.adj<-pseudoLog10(log(1/df$p.adj,2))
gplot<-ggplot(df,aes(x=p.adj, y=R, colour=n, size=k )) +
  geom_point() +
  expand_limits(x=0) +
  scale_y_continuous(
    limits = c(0,3),
    breaks = 0:3,
    guide = "prism_offset_minor",
    minor_breaks = pseudoLog10(rep(1:9, 4)*(10^rep(0:3, each = 9))),
    labels =  function(lab){
      do.call(
        expression,
        lapply(paste(lab), function(x) {
          if(x==0){bquote(bold("0"))}else{bquote(bold("10"^.(x)))}
        })
      )
    }
  )+
  scale_x_continuous(
    limits = c(0,3.2),
    breaks = 0:3.2,
    guide = "prism_offset_minor",
    minor_breaks = pseudoLog10(rep(1:9, 4)*(10^rep(0:3.2, each = 9))),
    labels = function(lab){
      do.call(
        expression,
        lapply(paste(lab), function(x) {
          if(x==0){bquote(bold("0"))}else{bquote(bold("10"^.(x)))}
        })
      )
    }
  )+
  scale_size(range = c(0,20))+
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    axis.ticks.length = unit(5,"pt"),
    axis.ticks = element_line(size = 0.7),
  )+
  geom_text_repel(
    data = df[df$R>=pseudoLog10(Thresh*1000) & df$p.adj>=pseudoLog10(log(1/0.05,2)),],
    aes(label=GeneSymbol),hjust=0, vjust=0, size=8, color="darkred"
  )+
  coord_cartesian(clip = "off")+
  geom_hline(yintercept=pseudoLog10(Thresh*1000), linetype="dashed", 
             color = "red", size=0.4)+
  geom_vline(xintercept=pseudoLog10(log(1/0.05,2)), linetype="dashed", 
             color = "red", size=0.4)+
  labs(x="log2(1/FDR)", 
       y="R (‰)", 
       title="",
       colour="log10(Gene Related Articles No.)", size="Osteoblast Differentiation\nRelated Articles No.")+
  annotation_custom(
    grob = grid::linesGrob(gp=gpar(lwd=2)),
    xmin = -Inf,
    xmax = -Inf,
    ymin = 0,
    ymax = 3
  )+
  annotation_custom(
    grob = grid::linesGrob(gp=gpar(lwd=2)),
    xmin = 0,
    xmax = 3.3,
    ymin = -Inf,
    ymax = -Inf
  )
svg(filename = "literature.mining.plot.svg",width = 10, height = 7)
print(gplot)
dev.off()
#Results Output
write.csv(literature.mining,"Literature_Mining.csv",row.names = FALSE)
####|01.Make Wordcloud from GO terms|####
#read protein annotation list from string-db website
wordcloud<-as.data.frame(literature.mining[literature.mining$R>=Thresh&literature.mining$p.adj<=0.05,])
#combine all text
temp<-c()
for (i in c(11,13)) {
  temp<-paste0(temp,paste0(wordcloud[,i],collapse = " "), collapse = " ")
}
#text clean & measure frequency
clean<-paste(c("\\;","\\[","\\]","\\)","\\(","/","-", "~",",",":","[[:digit:]]+"),collapse = "|")
temp<-gsub(clean," ",temp)
temp<-unlist(strsplit(temp," |,"))
temp<-tolower(temp)
word.freq<-table(temp)
word.cloud<-data.frame(word=names(word.freq),counts=as.integer(word.freq))
word.cloud<-word.cloud[order(as.integer(word.freq),decreasing = TRUE),]
word.cloud<-word.cloud[-1,]
clean<-c(
  "into","restricted","activated","activity", "alpha", "beta", "by", "cell", "cellular", "factor", "gamma", "i", "ii", "iii", "involved", "mediated", "negative", "pathway", "positive", "process", "region", "regulation", "regulatory", "response", "signaling", "stimulus", "through", "vi", "via", "vii", "viii", "ix","go","kegg","to","in","on","a","b","c","d","e","f","g","h","i","g","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","an","of","alpha","beta","gamma", "and", "side", "up", "down", "protein", "dna", "rna","gene","specific","like", "specification","from","component","space","activation","core","signal","type", "binding","expression", "system", "templated", "organic", "canonical","complete","genome","proteome",stopwords()
)
word.cloud<-word.cloud[!match(word.cloud[,1],clean,nomatch = 0)>0,]
#generate wordcloud
svg(filename = "Word.Cloud.plot.svg", width = 7, height = 7)
par(mar=c(0,0,0,0))
set.seed(1)
wordcloud(words = word.cloud$word, freq = word.cloud$counts, min.freq=2, max.words = Inf,
          rot.per=0.1, random.order=FALSE,
          colors=c(
            rep("#A6CEE3",2),
            rep("#33A02C",3),
            rep("#6A3D9A",12),
            rep("blue",40),
            rep("red",40)
          ))
dev.off()
set.seed(Sys.time())
####|02.Gene Ontology (GO) enrichment (Over Representation Test)|####
FEA<-data.frame()
for (j in 1:3) {
  FEA_GO_OR <-enrichGO(
    gene          = literature.mining[literature.mining$R>=Thresh&literature.mining$p.adj<=0.05,]$ensembl_gene_id,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENSEMBL",
    ont           = c("BP","MF","CC")[j],
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.02,
    readable      = TRUE
  )
  if(is.null(FEA_GO_OR)){next()}
  temp<-as.data.frame(FEA_GO_OR@result)
  temp<-temp[temp$p.adjust<=0.05,]
  temp<-temp[as.numeric(unlist(lapply(strsplit(temp$GeneRatio,"/"), `[[`, 1)))>=2,]
  temp<-temp[temp$ID%in%GOMostDescend(temp$ID),]
  temp<-OCcutOff(geneSets=FEA_GO_OR@geneSets, X=temp, value=0.25)
  geneSets<-list()
  for (i in 1:nrow(temp)) {
    geneSets[[i]]<-unlist(strsplit(temp$geneID[i],"/"))
    names(geneSets)[i]<-temp$ID[i]
  }
  temp<-OCcutOff(geneSets=geneSets, X=temp, value=0.6)
  temp<-cbind(Category=rep(c("BP","MF","CC")[j],nrow(temp)),temp)
  FEA<-rbind(FEA,temp)
}
KEGGgenels <- getBM(
  mart=Ann_0mart,
  attributes=c("ensembl_gene_id","entrezgene_id",
               "external_gene_name"),
  filter="ensembl_gene_id",
  values=literature.mining[literature.mining$R>=Thresh&literature.mining$p.adj<=0.05,]$ensembl_gene_id,
  uniqueRows = TRUE)
FEA_KEGG_OR <- enrichKEGG(gene         = KEGGgenels$entrezgene_id,
                          organism     = 'mmu',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.02,
                          pAdjustMethod = "BH")
temp<-as.data.frame(FEA_KEGG_OR@result)
temp<-temp[temp$p.adjust<=0.05,]
temp<-temp[!temp$ID%in%c("mmu05205","mmu05146","mmu05144","mmu05217","mmu05224"),]
temp<-cbind(Category=rep("KEGG",nrow(temp)),temp)
for (i in 1:nrow(KEGGgenels)) {
  temp$geneID<-str_replace(temp$geneID, as.character(KEGGgenels$entrezgene_id[i]), KEGGgenels$external_gene_name[i])
}
FEA<-rbind(FEA[],temp)
FEA<-cbind(
  FEA,
  List.Hit=(as.numeric(unlist(lapply(strsplit(FEA$GeneRatio,"/"), `[[`, 1)))/
    as.numeric(unlist(lapply(strsplit(FEA$BgRatio,"/"), `[[`, 1))))*100,
  Term=paste0(FEA$ID,"_",FEA$Description)
)
write.csv(FEA,file = "FEA.csv",row.names = FALSE)
#Make file for ring chart in Cytoscape
cyto.ring<-literature.mining[literature.mining$R>=Thresh&literature.mining$p.adj<=0.05,][,c(1,4,14)]
cyto.ring<-cbind(cyto.ring, MW_Group=rep("",nrow(cyto.ring)))
cyto.ring[cyto.ring$Mass>133000&cyto.ring$Mass<186000,"MW_Group"]<-"134~185"
cyto.ring[cyto.ring$Mass>105000&cyto.ring$Mass<135000,"MW_Group"]<-"106~134"
cyto.ring[cyto.ring$Mass>66000&cyto.ring$Mass<86000,"MW_Group"]<-"67~85"
cyto.ring[cyto.ring$Mass>21000&cyto.ring$Mass<45000,"MW_Group"]<-"22~44"
for (i in 1:nrow(FEA)) {
  cyto.ring<-cbind(cyto.ring, rep(0,nrow(cyto.ring)))
  colnames(cyto.ring)[4+i]<-FEA$ID[i]
  cyto.ring[tolower(cyto.ring$GeneSymbol)%in%tolower(c(unlist(str_split(FEA[i,"geneID"],"/")))),
    4+i]<-1
}
write.csv(cyto.ring,file = "Cyto.ring.csv",row.names = FALSE)
BarChartCol<-"#490C41&quot;,&quot;#43AE69&quot;,&quot;#159FE8&quot;,&quot;#4EFEBF&quot;,&quot;#130986&quot;,&quot;#832554&quot;,&quot;#B9CE5D&quot;,&quot;#B1D125&quot;,&quot;#A64552&quot;,&quot;#D7B725&quot;,&quot;#D6C97F&quot;,&quot;#D7372C&quot;,&quot;#DE4616&quot;,&quot;#356B28&quot;,&quot;#80E5CF&quot;,&quot;#4B6D3F&quot;,&quot;#D398E0&quot;,&quot;#E29479&quot;,&quot;#4AD68A&quot;,&quot;#CAC817&quot;,&quot;#B3CA4F&quot;,&quot;#C60E47&quot;,&quot;#5D3D5E&quot;,&quot;#B6B538&quot;,&quot;#0F1C9E&quot;,&quot;#15677B&quot;,&quot;#96A960&quot;,&quot;#0D66D9&quot;,&quot;#881E58&quot;,&quot;#ED054D&quot;,&quot;#15FA99&quot;,&quot;#4A12BB&quot;,&quot;#E631E2&quot;,&quot;#A1DB27&quot;,&quot;#A72A5D&quot;,&quot;#BD8BAC&quot;,&quot;#02E3B5&quot;,&quot;#4791B4&quot;,&quot;#217B43&quot;,&quot;#AC11DC&quot;,&quot;#6EB158&quot;,&quot;#13C74D&quot;,&quot;#D15EBB&quot;,&quot;#995B07&quot;,&quot;#7F8B41&quot;,&quot;#321846&quot;,&quot;#2473D6&quot;,&quot;#5B0D6B&quot;,&quot;#A331DC&quot;,&quot;#D8A016&quot;,&quot;#4212B4&quot;,&quot;#47DC20&quot;,&quot;#83135C&quot;,&quot;#797BD1&quot;,&quot;#CBDD7F&quot;,&quot;#313B74&quot;,&quot;#44BD49&quot;,&quot;#F23E2A&quot;,&quot;#BDBE70&quot;,&quot;#2F9015&quot;,&quot;#291FA5&quot;,&quot;#CE1513&quot;,&quot;#E09FD7&quot;,&quot;#61221C&quot;,&quot;#AF6B2A&quot;,&quot;#7CC898&quot;,&quot;#6FBC29&quot;,&quot;#4AAB01&quot;,&quot;#724B19&quot;,&quot;#EB88C7"
BarChartCol<-unlist(str_split(BarChartCol, '&quot;,&quot;'))
#Plotting FEA
gplot<-FEA                                                               %>%
  ggplot(aes(x=List.Hit, y=Term, colour=p.adjust, size=Count)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="List Hits (%)", y="", 
       title="",
       colour="P value", size="Gene Number")+
  scale_y_discrete(limits=c(FEA$Term))
for (i in 1:nrow(FEA)) {
  gplot<-gplot+annotation_custom(
    grob = rectGrob(gp=gpar(fill=BarChartCol[i],lwd=0)),
    xmin = -0.5,
    xmax = -0.1,
    ymin = 0.6+1*(i-1),
    ymax = 1.4+1*(i-1)
  )
}
svg(filename = "FEA.plot.svg", width = 10, height = 10)
print(gplot)
dev.off()

# Seurat codes of mine https://github.com/ximbao/AI_BrCa_Stemness/blob/main/scripts/seurat_objects_QCs_clusterization/Reference_dataset_breast_GSE176078/Ref_dset_standart_seurat_pipeline.r


# 1. How to make a random matrix

paises.data <- matrix(c(runif(63,0,100)),ncol=7,byrow=TRUE)
colnames(paises.data) <- c("Pop.","IDH","Ext.terr.", "homens","mulheres", "alfabetizados", "Idade.do.pais")
rownames(paises.data) <- c("brasil","EUA","Irlanda", "Noruega", "Jamaica", "Congo",
                           "México", "Egito", "China")

# 2. entering and manipulate data tips

#http://131.111.177.41/statistics/R/enteringdata.html


# 3. acesscing elements in a dataframe 

#https://rstudio-pubs-static.s3.amazonaws.com/48099_9f2a7ab501114ecaa49749eea8ba4881.html


# 4. to compare missing characters between two vectors 

vector.X[!(vector.X %in% vector.y)] # atencao nao ordem dos vetores !!!
#or
setdiff(vector.X,vector.y)
#or
vector.y %in% vector.X
#NOTE: I used this when I needed to compare a relabels at the 109 line of the expression code 




# 5. how to change a single rowname

rownames(teste)[rownames(teste) == "SC11.014BEB.133.5.6.11"] <- "EB"
#"SC11.014BEB.133.5.6.11" is a rowname and "EB" is the change we want to do



# 6. hot to occult lines of code (click at <-> ou v right of the number of the line (42))

#---------------------




#---------------------

# 7. download dataframe direct into R

#https://www.datacamp.com/community/tutorials/r-data-import-tutorial?utm_source=adwords_ppc&utm_campaignid=10267161064&utm_adgroupid=102842301792&utm_device=c&utm_keyword=&utm_matchtype=b&utm_network=g&utm_adpostion=&utm_creative=332602034358&utm_targetid=dsa-473406573035&utm_loc_interest_ms=&utm_loc_physical_ms=1001763&gclid=EAIaIQobChMIo6WkjrvB6gIVCQqRCh0syAKHEAAYASAAEgKO9_D_BwE#csv


# 8. how to set a working directory or check the current directory

#setwd(“…”) will set the current working directory to a specific location
#getwd() will print out the current directory.


# 9. how to know how many elements are the same between two vectors. NOTE:
# the intersect fucntion dont reveal if the elements are in the same order.



# 9.1 how to know the order of obs among datas

identical(rownames(data1),rownames(data2))



# 10. how to check NA values in a matrix or dataframe

apply(is.na(data),2,sum)
# 'is.na(data)' is the data argument; the '2' means rows and colun; sum is the sum of all values in each colun
# if we put just '1' at the second argument, we could check vectors instead

# 10. how to check the total count of NA in a dataframe

table(is.na(data))

# 11. how to scale the values of a matrix between 0 and 1
aiai<- matrix(c(runif(10,0,100)),ncol=2,byrow=TRUE)
min(aiai)
aiai.1 <- aiai - min(aiai)
aiai
aiai.2 <- aiai.1 / max(aiai.1)
# to visualise the diference, print aiai matrix, aiai.1 and aiai.2 and compare the diferences



# 12. How to generate a randon number 

sample(1:6,2) # it could be considered a roll of two six face dice
# (interval,amout)


# 13. How to select random observations in a data set
sample(1:length(data$colum),length(data$colum)*0.7)
# in this case, we are geting 70% from our data randomly


# 14. how to subset a dataframe through another smaller (but intersected) dataframe -> objective: ordering by obs

data.sub <- data1[rownames(data1)%in%rownames(data2),]
# data1 == larger dataframe (more obs)


# 15. getwd e list.filelist -> revela o atual diretório; revela as files do atual diretorio


# 16. como checar quantidades das diversas categorias/levels de uma coluna

#table(data$colum )❤


# 17. apagando colunas em um dataframe

#data[,c(2,3,4)] <- NULL # apagando colunas 2, 3 e 4

# 18. como usar string para extrair caracteres dos elementos de um vetor 

#colnames(data.complete.sub) <- str_split_fixed (as.character(colnames(data.complete.sub)), "_", 2)[,1] # "corta em 2 e fica com a primeira parte" 


# 19. plotar 2 graficos side-by-side

# 1º par(mfrow=c(1,2)) 2º plot.1  3º plot.2


# 20. bom pra salvar plot (ggplot)
# ggsave(file=caminho,width = 6, height = 4, dpi = 300) # salvar sempre em pdf
version



# 21. save headtmap
#save:
#pdf(hm.NPCs.diif.met)="/media/hd/maycon/meu.projeto.tcc/teste.estatistico.9vs126smp/hm.NPCs.diif.met", width=1000, height=600)
#draw(hm.plot, annotation_legend_side = "right", heatmap_legend_side = "right")
#dev.off()




# 22. mudando colnames
names(data)[names(data) == "old colname"] <- "new colname"



# 23. filtrando valores de uma coluna

df %>% filter(<Column Name Here> != "Colum value here") 
filtro <- subset(df, !(<Column Name Here>) %in% c("Colum values here")) # prefiro este

# 24. como salvar o "enviroment" do prompt de comando do PC do lab 
#save(list=ls(), file="/path/filname.Rda")


# 25. vizuliar jobs no terminal do Rstudio
# Ir no terminal e digitar htop


# 26. ordernar classes para melhorar vizualizaçao do plot
mydf$colun <- factor(mydf$colun, levels = c())

# 27. removendo varias colunas de uma unica vez de acordo com a terminacao da suas strings
library(dplyr)
data.CPM <- dplyr::select(data.CPM, -ends_with(c("ção.","back."))) 


# 28. checando tamanho do arquivo 
file.info("pathway")$size #just to check file size

# 29. outra forma de baixar dados suplementares do GEO (Renan que me passou)
sfiles <- getGEOSuppFiles("GSE104210", baseDir ="/media/hd/maycon/MESTRADO/projeto.principal/amostras/Pediatric_fetal_brain_tumors_GEO" ) #Download das SuppFiles
fnames <- rownames(sfiles) #Listagem das SuppFiles, checar qual é a de interesse
GSE104210.450K <- read.delim(fnames[1],header=TRUE) #Importar os dados 


# 30. para checar pData fora do R, acessar link a baixo e inserir o GSE no search
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE119834&o=acc_s%3Aa



# 31. ler file xlsx (excel)
library(openxlsx)
data.readed <- read.xlsx("directory")




# 32. exemplo de ajuste de legendas do eixo X
ggplot(pdata.merged[!pdata.merged$TCGA.subtype %in% NA,], aes(x=TCGA.subtype,y=mDNAsi.GLASS.Neuron)) + geom_boxplot(fill=c("purple","red","green","orange","yellow","blue","cyan3","purple","red","green","orange","yellow","blue","cyan3"),outlier.colour = NA)+
  xlab("Glioma subtype") + ylab("Neuron mDNAsi") + theme(axis.title.x = element_text(colour="black", size = 12), axis.title.y = element_text(colour="black", size = 12)) + ggtitle("Neuron mDNAsi by GLASS samples") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle = 45, size = 10 ,vjust=1,hjust=1)) + facet_wrap(recurrent_status~.)


# 33. como substituir caracter em uma string
gsub("out", "in", vector) #ex
gsub("R","D",s$TESTE)


# 34. teste estatistico 
t.test(df[df$group %in% "group1", "stemness.index.column"], df[df$group %in% "group2", "stemness.index.column"]) # tathi sugestão
# web site com mais exemplos



# 35. teste estatistco + plot
ggboxplot(tcga.glioma.merged[!tcga.glioma.merged$Subtype_Selected %in% "NA",], x = "Subtype_Selected", y = "mDNAsi.TCGA.data.Neuron", color = "Subtype_Selected", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(tcga.glioma.merged$mDNAsi.TCGA.data.Neuron), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1.2)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")
#fonte: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/



# 36. download de dados pelo Terminal
# 1. ssh maycon@143.107.202.128
# 2. temp (password)
# 3. cd "directory" 
# 4. wget "http..." (NÃO entrar no R!!)



# 37. como salvar um dataframe em xlsv (excel)

library(writexl)
write_xlsx(probes.Annot.GSC,"/media/hd/maycon/MESTRADO/InfiniumMethylation_probeprospect/probes.Annot.GSC.xlsx")


# 38. criando matrix 
mymat <- matrix(nrow=30, ncol=30)


# 39 paciente vs mDNAsi em ordem crescente 
ggplot(ss,aes(reorder(id.16, +mDNAsi, sum), mDNAsi))+
  geom_col(width=1)

# 40. reodenando levels (trocando ordem das categorias/levels)
mydf$Subtype_Selected <- factor(mydf$Subtype_Selected, levels =c("GBM_LGG.Codel", "GBM_LGG.G-CIMP-high","GBM_LGG.G-CIMP-low", "GBM_LGG.Classic-like","GBM_LGG.Mesenchymal-like","GBM_LGG.LGm6-GBM","GBM_LGG.PA-like", "GBM_LGG.NA")) #reordenando os leveis


# 41. alterando nome dos levels (trocando nome das categorias/levels)
mydf$Subtype_Selected <- mapvalues(mydf$Subtype_Selected, from = c("GBM_LGG.Codel", "GBM_LGG.G-CIMP-high","GBM_LGG.G-CIMP-low", "GBM_LGG.Classic-like","GBM_LGG.Mesenchymal-like","GBM_LGG.LGm6-GBM","GBM_LGG.PA-like", "GBM_LGG.NA"), to = c("codel","G-CIMP-high","G-CIMP-low","classic","mesenchymal","LGG.LGm","PA-like","NA")) #alterando o nome dos leveis



# 42. confusion.matrix 

install.packages("cvms")
library(cvms)
# seguindo as ideias de https://rdrr.io/cran/cvms/man/plot_confusion_matrix.html

# fazendo um filtro apenas com as colunas de interesse
data.p.matrix.conf<- (pData.BRCA.dog.[,c("paper_BRCA_Subtype_PAM50","subtype")])


# transformando em character os lvs das colunas
data.p.matrix.conf$paper_BRCA_Subtype_PAM50 <- data.p.matrix.conf$paper_BRCA_Subtype_PAM50 %>% as.character()

data.p.matrix.conf$subtype <- data.p.matrix.conf$subtype %>% as.character()

# removendo NA
data.p.matrix.conf <- na.omit(data.p.matrix.conf)
eval <- evaluate(
  data = data.p.matrix.conf,
  target_col = "paper_BRCA_Subtype_PAM50",
  prediction_cols = "subtype",
  type = "multinomial"
)

# Inspect confusion matrix tibble
eval[["Confusion Matrix"]][[1]]

# plotando confusion.matrix
plot_confusion_matrix(eval[["Confusion Matrix"]][[1]])
ggsave(file="/media/hd/maycon/MESTRADO/Dog_data/matrix.confusion.pam50.pdf",width = 6, height = 4, dpi = 300)


# 43. gerando sensibilidade e especificidade entre 2 classificoes (esses valores sao mais importantes/praticos do que o plot de confusion.matrix)
confusionMatrix(data= (as.factor(data$predction)), reference = as.factor(data$targed))

# exemplo de output
#          Class: Basal Class: Her2   Class: LumA Class  
#Sensitivity  0.9897       0.92683        0.4462     
#Specificity  0.9732       0.96726         0.9943   


# 44. explicando tcga barcode
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

# 45. tcga repositorio amostras
# https://portal.gdc.cancer.gov/repository

# 46. como fazer uma lista de characters pra colar em outros aplicativos/plataformas/softwares 
ENSEMBL.BRCA #vetor com id de genes
ENSEMBL.BRCA <- matrix(ENSEMBL.BRCA) #transformando em matrix
write.csv(ENSEMBL.BRCA, file ="/media/hd/maycon/MESTRADO/Dog_data/ENSEMBL.BRCA.csv", row.names=FALSE) #salvando em .csv para acessar a coluna com os nomes -> dar Ctrl+C na coluna -> colar essa lista em algum aplicativo ou software fora do R.

# 47. como inserir 'virgula' ou qualquer outra coisa entre os elementos de um vetor
library('stringr')
list.to.paste.ENSEMBL.BRCA<- str_c(ENSEMBL.BRCA, collapse=',')    
#'ENSEMBL.BRCA' é o vetor com characteres 


# 48. baixando pData (dados clinicos) de amostras do TCGA
GDCquery_clinic(project = "", type = "clinical", save.csv = FALSE)
# https://rdrr.io/bioc/TCGAbiolinks/man/GDCquery_clinic.html


# 49. baixando idat (pra normalizar) de amostras do TCGA
query <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
                  legacy = TRUE,
                  data.category = "Raw microarray data",
                  platform = c("Illumina Human Methylation 450"),
                  file.type = "idat")
                  #sample.type = "Solid Tissue Normal")
GDCdownload(query)



# 50. como excluir colunas de um dataframe 
matrix.meth.global[,!(colnames(matrix.meth.global) %in% GBM_LGG.pData$barcode)]


# 51. como saber a frequencia de cada classe em uma coluna 
count(data$colum)
library(plyr)


# 52. remove colunas com characteres finais especificos
testete <- select(data, -ends_with(c("terminacao"))) 
# +ends_with mantem as colunas com o final espeficio no nome


# 53. como separar strings por ";" "," etc -> em colunas {excel}
#https://support.microsoft.com/en-us/office/split-text-into-different-columns-with-the-convert-text-to-columns-wizard-30b14928-5550-41f5-97ca-7a3e9c363ed7
# Excel -> dados -> texto para coluna ...

# 54. como apagar celulas em branco no excel 
#https://www.techtudo.com.br/noticias/2019/03/como-excluir-todas-as-linhas-em-branco-no-excel-de-uma-vez.ghtml
# selecionar celulas preenchidas e celulas em branco -> ferramenta 'localizar e selecionar' -> 'ir para especial' -> seleciona as celulas em branco -> excluir

# 55. como juntar todos os elementos em uma so coluna no excel
# faço na mao mesmo... copio as colunas e vou colando em baixo da primeira.


# 56. como identificar rownames de amostras de acordo com 1 level de uma coluna
X <- rownames(data)[which(data$"column name" == "level")]


# 57. removendo elementos especificos de um vetor
X <- setdiff(vetor,elementos a remover)
GSCsig.probes <- setdiff(GSCsig.probes, c("cg06007645","cg03643137","cg01291854"))


# 58. complexHeatmap 
https://www.biostars.org/p/286187/#286507

# 59. sugestao de scrip para complexHeatmap
  hm_EMT <- Heatmap(heat.EMT.sig,
                    cluster_columns=T,
                    cluster_rows = T,
                    clustering_method_rows="complete",  
                    show_row_names = T,
                    row_names_gp = gpar(fontsize = 10),
                    show_column_names = T,
                    name = "Gene expression",
                    col= jet.colors(75),
                    row_title = "genes",
                    #row_names_gp = gpar(fontsize = 12),
                    column_title = "EMT genes signature 63/76 (Mak, MP et al. 2016)",
                    #split=label.EMT.genes$E_ou_M,
                    #row_names_side="left",
                    
                    #top_annotation = top.anno,
                    #row_names_gp = gpar(fontsize = 8),
                    #column_names_gp = gpar(fontsize = 12),
                    heatmap_legend_param = list(
                      color_bar = 'continuous',
                      legend_direction = 'vertical',
                      legend_width = unit(8, 'cm'),
                      legend_height = unit(5.0, 'cm'),
                      title_position = 'leftcenter-rot',
                      title_gp=gpar(fontsize = 12, fontface = 'bold'),
                      labels_gp=gpar(fontsize = 12, fontface = 'bold')))

RowAnn <- HeatmapAnnotation(df = label.EMT.genes, col=list("E_ou_M" = c("Epi" = "green", "Mesq" = "yellow")),
                            show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
                            na_col= "white",which = "row",show_legend =T)


draw(hm_EMT + RowAnn) 



# 60. removendo duplicatas
df[!duplicated(df$colunm),]


# 61. como cortar string do inicio pro final ou do final para o incio
library(stringr)
# faça o primeiro corte (first.splt) -> faca o segundo corte (second.splt)
second.splt <- str_split_fixed(as.character(teste$first.splt), "[0-9][0-9][0-9]", 2)[,]
# caso precise cortar no meio da string -> fazer 2 cortes


# 62. alterando nome de leveis/lv/categorias
library(plyr)
mydf$Subtype_Selected <- mapvalues(mydf$Subtype_Selected, from = c("GBM_LGG.Codel", "GBM_LGG.G-CIMP-high","GBM_LGG.G-CIMP-low", "GBM_LGG.Classic-like","GBM_LGG.Mesenchymal-like","GBM_LGG.LGm6-GBM","GBM_LGG.PA-like", "GBM_LGG.NA"), to = c("codel","G-CIMP-high","G-CIMP-low","classic","mesenchymal","LGG.LGm","PA-like","NA")) #alterando o nome dos leveis



# 63. como mudar a ordem das colunas
teste[,c(1,2,3,4)] # ordem normal
teste[,c(2,1,3,4)] # ordem que voce quer



# 64. como saber a distribuicao de NA por linha (ou probe) em uma matrix
############## Descobrindo a distribuicao de probes:NA para saber quantas selecionar pro kNN ##############
colSums(is.na(volcano.NAs))
rowSums(is.na(volcano.NAs)) # esse é melhor pra usar nesse caso
NA.maior.5 <- data.frame(rowSums(is.na(volcano.NAs)) > 5)
table(NA.maior.5$rowSums.is.na.volcano.NAs.....5) 
#FALSE  TRUE 
#17620  1173 

NA.maior.5.T <- subset(NA.maior.5,rowSums.is.na.volcano.NAs.....5 %in% TRUE)
table(duplicated(rownames(NA.maior.5.T)))
#FALSE 
#1173

histogram(rowSums(is.na(volcano.NAs[rownames(NA.maior.5.T),])))
########################################## descobri que nao precisa fazer isso - dá pra rodar todas essas 18mil probes no kNN  SOOORRY :/ #############################

# Mas fica de exemplo pra se algum dia eu precisar desse raciocinio !



# 65. como criar filtro por intervalo
library(dplyr)
df$filtro <- "blabla"
b.b <- df%>%filter(between(column, "de tanto", "a tanto")) 
df[rownames(b.b),"filtro"] <- "nome da classe"



# 66. pra quando se quer informacao de amostras especifias de um dataframe
df[df$coluna para selecionar as amostras %in% categoria das amostras,]$coluna com a informacao que vc quer dessas amostras

# Ex: scores2[scores2$Class.added %in% "neural.stem.cells",]$barcode



# 67. maneira de normalizar seus dados (de 0-1) em uma matrix 
normalize <- function(x) { 
  x <- sweep(x, 2, apply(x, 2, min)) 
  x <- sweep(x, 2, apply(x, 2, max), "/") 
  2*x - 1
}

normalize(df)



#  68. teste estatistico que tathi me passou uma vez 
t.test(df[df$group %in% "group1", "stemness.index.column"], df[df$group %in% "group2", "stemness.index.column"]) 




#  69. Processamento de IDAT file 
######
# 07/06/2021 (segunda-feira)
# Normalizacao dos GBM, LGG e non-tumor(amostras) do TCGA
library(BiocManager)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(shinyMethyl)
library(bumphunter)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)

setwd("/media/hd/maycon/MESTRADO/rascunho.objetos/Normalizacao_GBM.LGG_TCGA")

untar("/media/hd/maycon/MESTRADO/projeto.principal/amostras/GSE109399_GSC_and_others/?acc=GSE109399&format=file")
library(GEOquery)
idatFiles <- list.files("/media/hd/maycon/MESTRADO/projeto.principal/amostras/GSE92462_GSC_and_others", pattern = "idat.gz$", full = TRUE, recursive = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

### Step 1: Read idat files GBM samples ###

baseDir <- "/media/hd/maycon/MESTRADO/rascunho.objetos/Normalizacao_GBM.LGG_TCGA/TCGA.idat.files.to.process/GBM_tcga.idat.together"  #define the path where the idat files are located
idat.file <- read.metharray.exp(base = baseDir, targets = NULL, force = TRUE, recursive = T) 

# controle de qualidade dos idat.files (so pra vizualizarmos)
summary.idat <- shinySummarize(idat.file)
runShinyMethyl(summary.idat) 



### Step 2: Calculate p-values ###
#verifica estatisticamente se os sinais de metilacao sao confiaveis 

detP.idat <- detectionP(idat.file, type = "m+u") 
table(detP.idat > 0.05)
#FALSE     TRUE
#75091366   162994  

# FALSE sao as probes com p value adequado (pvalue > 0.05 no caso) e TRUE sao as probes com p value insatisfatorio. 


### Step 3: Preprocess ###

proc.idat.file <- preprocessNoob(idat.file, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference")  


### Step 4: mask probes that failed p-value detection ###
# valores com p values insatisfatorios

proc.idat.file.r <- ratioConvert(proc.idat.file)

is.na(assays(proc.idat.file.r)$Beta) <- (detP.idat[rownames(proc.idat.file.r), colnames(proc.idat.file.r)] > 0.05)


### Step 5: Get beta-values ###

beta.idat.file <- assays(proc.idat.file.r)$Beta
beta.idat.file  <- as.data.frame(beta.idat.file )  
save(beta.idat.file , file = "/media/hd/maycon/MESTRADO/rascunho.objetos/Normalizacao_GBM.LGG_TCGA/beta.idat.GBM.Rda") 

# # #    E n d    o f   p r o c e s s i n g    N o r m a l i z a t i o n    # # #





#  70. traducao de gene ID por biomart
library(biomaRt)

# de ensembl ->  gene_symbol
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
traduction <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters ='ensembl_gene_id', values =rownames(data), mart = ensembl)  #rownames(data) é o que vai ser convertido
dim(traduction) #815   2
data$ensembl_gene_id <- rownames(data)
data <- merge(data,traduction,by.y='ensembl_gene_id')
rownames(data) <- data$hgnc_symbol

# de gene_symbol ->  ensembl
library(biomaRt)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
traduction <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters ='hgnc_symbol', values = rownames(raw_GSE163899), mart = ensembl)


#  71. ggplot exemplo 1
ggplot(mydf[!mydf$Subtype_Selected %in% "NA",], aes(x=Subtype_Selected,y=mDNAsi.y)) + geom_boxplot(fill=c("purple","red","green","orange","yellow","blue","cyan3"),outlier.color = NA)+ geom_jitter (alpha=0.2)  +
  xlab("by Subtype_tumor") + ylab("mDNAsi   GSC") + theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), axis.title.y = element_text(colour="black", size = 12), axis.text.x = element_text(angle = 45, vjust= 0.7)) + ggtitle("Gliomas Stemness index - Glioma Stem Cell Signature") 


ggsave("/media/hd/maycon/MESTRADO/rascunho.objetos/plots/mDNAsi.GSC.TCGA.subtype.pdf",width = 6, height = 4)


#  72. como separar uma sequencia de letras
vector <- sapply(1:169, function(i) { substr(dara$Tom20.de.A..fumigatus,i,i) }
vec <- sapply('quantidade de letras',function(i) {substr(lista,i,i)}    
              
              
#  73. como criar dataframe com 2 vetores de comprimento distintos
df <- data.frame(vetor.MAIOR)
df$"outra coluna" <- NA
df$"outra coluna" <- vetor.MENOR



#  74. selecionando amostras aleatorias por categoria
library(caret)
selecting.randon.samples <- createDataPartition(pData.gliomas.Capper$Class.added, p=0.1, list=FALSE, times=1) # p = porcentagem do total de amostras que irao ser selecionadas / categoria na coluna
gliomas.GSC.BK <- pData.gliomas.Capper[selecting.randon.samples,]
dim(gliomas.GSC.BK)




# 75. arrumando tabulacao e extraindo as informacoes de strings em linhas (probe anno day)
#################################################################################
# Para as colunas 'transcriptIDs' e 'distToTSS' das probes NSC
load("/media/hd/maycon/MESTRADO/probe.anno/prob.anno.10_08_2021/probe.map.NSC.Rda")
# 'transcriptIDs'
library(splitstackshape)
View(cSplit(data.frame(probe.map.NSC[,"transcriptIDs"]), c("probe.map.NSC....transcriptIDs.."), sep=";", type.convert=FALSE))

test <- cSplit(data.frame(probe.map.NSC[,"transcriptIDs"]), c("probe.map.NSC....transcriptIDs.."), sep=";", type.convert=FALSE)

t.test <- t(test)
vect <- paste(t.test,1:954)
vect <- data.frame(vect)
library(stringr)
stringr::str_subset(vect, pattern = "NA_", negate = T)

NSC.trasncriptID <- data.frame(vect[!grepl("NA",vect$vect),])
str_split_fixed (as.character(NSC.trasncriptID$vect..grepl..NA...vect.vect....), "[:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][.][:alnum:]", 2)[,2]

NSC.trasncriptID <- str_split_fixed (as.character(NSC.trasncriptID$vect..grepl..NA...vect.vect....), "[ ]", 2)[,1]

NSC.trasncriptID
length(NSC.trasncriptID) #3848

# 'distToTSS'
library(splitstackshape)
test <- cSplit(data.frame(probe.map.NSC[,"distToTSS"]), c("probe.map.NSC....distToTSS.."), sep=";", type.convert=FALSE)

t.test <- t(test)

vect <- paste(t.test,1:954)
vect <- data.frame(vect)


NSC.distToTSS <- data.frame(vect[!grepl("NA",vect$vect),])
NSC.distToTSS<- str_split_fixed (as.character(NSC.distToTSS$vect..grepl..NA...vect.vect....), "[ ]", 2)[,1]

NSC.distToTSS
length(NSC.distToTSS) #3848

NSC.transcriptID.distToTSS <- data.frame(NSC.trasncriptID,NSC.distToTSS)
#########################################################################





# 76. Como criar um matrix
met <- matrix(rep(0,15),ncol = 5)
colnames(met) <- c("Sample1",
                   "Sample2", 
                   "Sample3",
                   "Sample4",
                   "Sample5")
rownames(met) <- c("cg26928153","cg16269199","cg13869341")


# 77. Como listar files em um diretorio que contém um padrao especifico no nome
dir(path = "result", pattern = "getMethdiff")  


# 78. como plotar graficos xxxxx
load("/media/hd/maycon/MESTRADO/rascunho.objetos/Carrapato_thales/carrapatos.Rda")
dim(carrapatos)
teste <- t(carrapatos)
dim(teste) #16 184
teste <- teste[-c(15,16),]
teste <- data.frame(teste)
teste$Month <- rownames(teste)
library(reshape2)
teste <- melt(teste, id.vars = "Month", measure.vars = c(2:184))

library(ggplot2)
ggplot(teste, aes(value, fill = Month)) + 
  geom_bar(position = 'identity', alpha = .3)



# 79. selecionar elementos unicos (remover duplicatas) em um vetor
unique(vector)


# 80. frequencia de repeticao de elementos (ex: genes em uma coluna/vetor)
GSC.gene.list <- data.frame(table(GSC.gene.list))


# 81. Como acessar a lista de linhas do heatmap
rownames(carrapatos[row_order(hm.plot),])
my_sample_col[column_order(hm_unique_probset), ] #ordem das amostras clusterizadas. Lembre-se de deixar os ronames do metadata usado no complex heatmap na forma de index.

#or
lista.heatmap <- data.frame(rownames(carrapatos[row_order(hm.plot),]))
colnames(lista.heatmap) <- "Animais" #trocando nome de coluna apenas
rownames(lista.heatmap) <- lista.heatmap$Animais
lista.heatmap <- carrapatos[rownames(lista.heatmap),]


# 82. como ordenar um dataframe/dataset a partir de uma de suas colunas
coursera_data[order(coursera_data$course_students_enrolled,decreasing = TRUE),]



# 83. como selecionar amostras randomicamente 
library(caret)
selecting.randon.samples<- createDataPartition(just.GSC$Class.added, p=0.66, list=FALSE, times=1)
just.GSC <- just.GSC[selecting.randon.samples,]
dim(just.GSC)



# 84. toop annotation do complexheatmap [como colorir sem deixar no aleatorio]
top.anno = HeatmapAnnotation(df = label, col=list("Stemcell.or.not" = c("Stem.cell" = "purple", "not.Stem.cell" = "black"),"samples.status" = c("bulk.brain" = "gold1", "GBM_tcga.proc" = "firebrick4", "LGG_tcga.proc" = "firebrick2", "GSC.not.xe.validation" = "dodgerblue4", "GSC.xe.validation" = "dodgerblue2", "NSC" = "darkolivegreen1", "NSC.validation" = "darkolivegreen3")),
                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
                             na_col= "white")


# 85. como add informacoes no dataframe (criando coluna)
meta$filtro <- "M05"
C <- meta[meta$sample %in% c("C1","C2","C3"),] 
meta[rownames(C),"filtro"] <- "C"
CF41 <- meta[meta$sample %in% c("CF41.01","CF41.02","CF41.03"),] 
meta[rownames(CF41),"filtro"] <- "CF41"



# 86. removendo colunas que tenham NA 
teste <- teste[ , colSums(is.na(teste)) == 0]



# 87. hisogramas justapostos ou sobrepostos
# Histogramas sobrepostos
par(mfrow = c(1,1))
hist(df_2$G2,
     col = rgb(1, 0, 0, 0.5),
     add=F)

hist(df_2$G3,
     col = rgb(0, 1, 0, 0.5),
     add=T)


hist(df_2$G4,
     col = rgb(0, 0, 1, 0.5),
     add=T)

# Hisogramas justapostos
par(mfrow = c(1,3))
hist(df_2$G2,
     col = rgb(1, 0, 0, 0.5),
     add=F)

hist(df_2$G3,
     col = rgb(0, 1, 0, 0.5),
     add=F)

hist(df_2$G4,
     col = rgb(0, 0, 1, 0.5),
     add=F)



# 88. como adicionar número de amostras por level 
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

ggplot(diamonds, aes(cut, price)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text")




# 89. como add uma string em outra
save(ss, file = paste(colnames(ss), ".Rda", sep=" "))
# add '.Rda' no final de uma string



# 90. UMAP plot (nunca plotei) (é tipo o tSNE, só que tem sido mais utilizado)
library(umap)
UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


# 91. selecionar colunas que contenham uma certa string
data[,grepl("string", colnames(data))]


# 92. ajustes de download pro Rstudio | conection
# (ajuste de download)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2) # turning to 
getOption('timeout') # 60
options(timeout=300) # turning to 300 



# 93. como plotar uma tabela em figura
library(gridExtra)
x <- df[,]
xx <- tableGrob(x)
grid.arrange(xx)


# 94. acessando pData do repositório GEO
library(GEOquery)
setwd("/media/hd/maycon/MESTRADO/Qualificação/Curva_sobrevida")
# For get pData from GSE92462.450K dataset 
gset.GSE61160.450K <- getGEO("GSE61160", GSEMatrix =TRUE, destdir = "/media/hd/maycon/MESTRADO/Qualificação/Curva_sobrevida")
dim(pData(gset.GSE61160.450K [[1]])) # 51 43
head(pData(gset.GSE61160.450K [[1]])[, 1:3])
pData.GSE61160.450K <- pData(gset.GSE61160.450K [[1]])
save(pData.GSE61160.450K, file="/media/hd/maycon/MESTRADO/projeto.principal/pData/pData.GSE61160.450K.Rda")


# 95. Contar nº de NAs em um dataframe/matrix
sum(is.na(df))


# 96. Visualização em tabela das estatisticas descritivas do df
library(skimr)
skin(df)


# 97. Looping pra gerar beta-values a partir de dados ja processados (unmetilated a methylated)
library(dplyr)
for(i in seq(1,150,2)) {
  teste <- mutate(teste, new_column = teste[,i+1]/(teste[,i]+teste[,i+1]+100))
  names(teste)[names(teste) == "new_column"] <- gsub("Signal_A", "Beta_value", names(teste)[i])
} 

dim(teste) #485512    225 
names(teste)[151:225]
teste[1:3,151] 
beta_values_part1 <- teste[,151:225] # selecionando os 75 beta values pois foram add no mesmo df que o dos sinais
# teste é o dataframe
# seq(1,122,2) -> o os argumentos 1 e 122 são o cromprimento do dataframe em que o script irá se repetir; o argumento 2 está dizendo para repetir o script de duas em duas colunas
# mutate é pra gerar a coluna automaticamente
# a fórmula do new_column é para gerar os beta-values (methylated/unmeth+methylated+100)
# o i e i+1 foi para a lógica de ir aplicando o looping de duas em duas colunas
# o names(teste) foi pra trocar o nome das colunas 


# 97. Checar quantidade de levels/classes/categorias de uma coluna a partir de um filtro
table(LGG_GBM_Kaplan_Meier[!LGG_GBM_Kaplan_Meier$`Vital.status.(1=dead)` %in% NA,]$IDH.status)


# 98. Transformando variavel quanti. em quali. 
df$column <- ifelse(df$column > 0,
                    yes = 1,
                    no = 0)
# Caso deseje fazer uma variável binária pra um determinado cutoff


# 99. Como extrair de uma string pela posição das letras
str_sub(pdata$glassID, start = 1, end = 15) #do primeiro elemento da string até o 15th


# 100. curva de sobrevida

############################
# OS Curve - Kaplan Meier" #
############################
library(survMisc)
pData$cut <- cut(pData$mDNAsi , breaks = 2) #transform variable quanti. into quali
library(survival)
km.model <- survfit(Surv(pData$os_months,
                         pData$vital_satus`) ~ 
                      pData$cut,type="kaplan-meier") #build de model 
summary(km.model)

plot(km.model, conf.int=F, xlab="OS.(months)",
     ylab="Death", main="KM-Model-Gliomas_XXX (XXXmodel stemness index)",
     col=c("green","red"), las=1, lwd=2, mark.time = F)

legend(18,0.95, legend=c("low-stemness","high-stemness"),lty=1,lwd=2,col=c("green","red"), bty="",cex=1)

survdiff(Surv(pData$os_months,pData$vital_satus) ~ pData$cut) #extract pvalue 

library(survminer) 
surv_median(km.model) #extract the median values

##################################################################################################



# 101. Algumas funções ninjas para manipular dataframe
#select() chooses columns.
#gather() combines many columns into two, one column to label the data and one to store the value
#group_by() makes subsets of the data, one subset for each value of the chosen column(s).
#summarize() does calculations on the subsets created by group_b


# 102. dummyzação 
df$column <- ifelse(wt.glass.NSC.500bp$column == "level",
                                          yes = 1,
                                          no = 0) #transformando em dummy
# Existe um pacote que dammyza varias colunasde uma vez



# 103. gerar figura da tabela para não precisar ficar olhando toda hora no view
library(knitr)
library(kableExtra)
df %>% #dataframe
  select(everything()) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 12)

# 104. contar numero de ocorrencias de um valor especifico em uma coluna 
length(which(df$column_name==value))



# 105. cotando frequencia das ocorrencias em uma coluna 
genes <- names(table(NSC_2000bp_promotor_ANNO$genesUniq))
freq <- as.data.frame(table(NSC_2000bp_promotor_ANNO$genesUniq))
rownames(freq) <- genes
freq
freq[freq$Var1 %in% "AC007731.4;SCARF2",]$Freq # 5 duplicatas
length(which(NSC_2000bp_promotor_ANNO$genesUniq=="AC007731.4;SCARF2")) # 5 duplicatas -> conferido!



# 106. plotar boxplot com apenas uma variavel 
library(ggplot2)
ggplot(pData[,], aes(x="",y=NSC_hypo_island_model)) + geom_boxplot()



# 107. como renomar um específico level na coluna
levels(pData.prediction_4801$Class.added)[levels(pData.prediction_4801$Class.added)=='Glio.neural'] <- 'Glio_neuronal'



# 108. Outra forma + sofisticada para  plotar sobrevida 
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(pData_wt$`Survival.(months)`,
                        pData_wt$`Vital.status.(1=dead)`) ~ 
                     pData_wt$cut,type="kaplan-meier"), data = pData_wt,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
# PS: o cut é dividido pela mediana dos valores (ex: high e low-stemness)



# 109. Criando um dataframe
employee <- c('John Doe','Peter Gynn','Jolie Hope')
salary <- c(21000, 23400, 26800)
startdate <- as.Date(c('2010-11-1','2008-3-25','2007-3-14'))
df <- data.frame(employee, salary, startdate)



# 110. Veen Diagram
library(ggVennDiagram)
x <- list(Hsi_MDA_MB_231_table04 = Hsi_MDA_MB_231_table04$Description, Hsi_cell_lines_table02 = Hsi_cell_lines_table02$Description); x
ggVennDiagram(x[1:2]) 
# Caso tiver um terceiro conjungo, é só inserir mais um vetor e printar x[1:3]



# 111. Identificar strings parecidas 
library(tidyverse)
library(stringdist)
table02 <- Hsi_cell_lines_table02$Description
table02 <- map_dfr(table02, ~ {
  i <- which(stringdist(., table02, "jw") < 0.10) # quanto menor essa porcentagem, o algoritimo agrupa strings muuuito semelhantes
  tibble(index = i, title = table02[i])
}, .id = "group") %>%
  distinct(index, .keep_all = T) %>% 
  mutate(group = as.integer(group)) %>% 
  as.data.frame()


# 112. Filtrando missing values (NA) por uma coluna
df %>%
  filter(is.na("column") == TRUE) %>%
  count()



# 113. rank bar plot (amostra por amostra, Tathi script)
barplot <- df[order(df$stemness.index),] #ordena

barplot$sample.id <- factor(barplot$sample.id, levels=as.character(barplot$sample.id)) #troca os levels DAS AMOSTRAS

library(ggplot2)

ggplot(barplot, aes(sample.id, stemness.index)) +
  
  geom_bar(stat="identity", aes(fill=factor(group), position="fill")) +
  
  scale_fill_manual(values=c("color1", "color2", "color3","color4"), name="Legend") +
  
  ggtitle("Stemness Index") +
  
  ylab("Stemness Index") + 
  
  xlab("bla bla") + 
  
  theme(axis.text.x = element_blank())


# 114. Como contabilizar frequencia de probes anotados pra genes em ilhas CpG 
HM450.hg38.gencode.v36 <- read.table(gzfile("/media/hd/maycon/MESTRADO/InfiniumMethylation_probeprospect/HM450.hg38.manifest.gencode.v36.tsv.gz"),header = T) # carregando as anotações do manifesto da illumina

suas_probes_info <- subset(HM450.hg38.gencode.v36,HM450.hg38.gencode.v36$probeID %in% suas_probes) #filtrando apenas info. pra sua lista de probes
suas_probes_info <- suas_probes_info[suas_probes_info$CGIposition %in% "Island",] #filtrando apenas ilhas CpG
suas_probes_info <- suas_probes_info [!suas_probes_info$genesUniq %in% NA,] #removendo NA


gene_ID <- names(table(suas_probes_info$genesUniq)) #nome dos genes
freq <- as.data.frame(table(suas_probes_info$genesUniq)) #dataframe com as frequencias (linha=gene;valor=nº probes)
rownames(freq) <- gene_ID #adicionando os gene_ID na tabela de frequencia

dim(freq[freq$Freq >= 2,]) #nº genes anotados para 2 probes ou mais
dim(freq[freq$Freq >= 3,]) #nº genes anotados para 3 probes ou mais ....


# 115. Para add linha de regressao no ggplot
geom_smooth(method = "lm", se = FALSE, colour="black")



# 116. Survival analysis 
# Filtrando apenas as amostras que quero plotar na curva kaplan-meier 
library(survival)
load("/media/hd/maycon/MESTRADO/Qualificação/Pipeline_quali/mDNAsi_GSC_refinado_TCGA_gliomas.Rda")
pData_wt <- pData[pData$IDH.status %in% "WT",]
# Transformando minha variável continua em categoria atráve da mediana
pData_wt$cut <- NA
pData_wt[pData_wt$NSC_hypo_island_model >= median(pData_wt$NSC_hypo_island_model),]$cut <- "high-stemness"
pData_wt[!pData_wt$NSC_hypo_island_model >= median(pData_wt$NSC_hypo_island_model),]$cut <- "low-stemness"
km.model <- survfit(Surv(pData_wt$`Survival.(months)`,
                         pData_wt$`Vital.status.(1=dead)`) ~ 
                      pData_wt$cut, stype=1, ctype=1) #build de model

# plotando a curva de sobrevida
plot(km.model, conf.int=F, xlab="Tempo desde o diagnóstico (meses)",
     ylab="Sobrevida  (%)", main="",
     col=c("firebrick1","seagreen4"), las=1, lwd=2, mark.time = F)
legend("topright", cex=1,
       legend=c(
         paste("(n=",km.model$n[1],")"),
         #paste("Intermediate-frequency (n=",fit$n[2],")",sep=""),
         paste("(n=",km.model$n[2],")")),
       col=c("firebrick1", "seagreen4", "dodgerblue", "grey"),
       ext.width=1,
       lty=c(1, 1,1),
       lwd=3,
       title="TCGA",
       box.lwd=3,
       bg="white")
legend(18,0.95, legend=c("baixo-stemness","alto-stemness"),lty=1,lwd=2,col=c("seagreen4","firebrick1"), bty="",cex=1)

# Estatísticas 
survdiff(Surv("tempo", "status") ~ "category variable", data = "your dataframe")  # qui quadrado que gera o p-valor da curva de kaplan meier
surv_median(survfit(Surv("tempo", "status") ~ "category variable",type="kaplan-meier") # da a diferença das medianas da curva de kaplan meier
coxph(Surv("tempo", "status") ~ "variable 1" + "variable 2", data = "your dataframe") # pode adicionar so a variável desejada ou todas as variaveis que você gostaria de testar algum tipo de colinearidade (essa função é uma regressão cox)

# END #



# 117. quando gunzip nao funcionar, tente gzfile()
gzfile("HM450.hg38.manifest.tsv.gz")


# 118. loop basic structure
for(i in dice) {
  if (teste[i,] == 6) {
    print(paste("The dice number is", i, "Yahtzee!"))
  } else {
    print(paste("The dice number is", i, "Not Yahtzee"))
  }
}


# 118. criar diretorio 
# Checa se já existe pasta com mesmo nome
# Caso não, este script cria o diretório e entra na pasta
if(file.exists("/media/hd/maycon/MESTRADO/Single_Cell/GSE176078")){
  setwd("/media/hd/maycon/MESTRADO/Single_Cell/GSE176078")
}else{
  dir.create("/media/hd/maycon/MESTRADO/Single_Cell/GSE176078")
  setwd("/media/hd/maycon/MESTRADO/Single_Cell/GSE176078")
}


# 119. baixando e dezipando arquivos
library(HelpersMG)
library(GEOquery)
# Aumentando limite de tempo para download
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100000000) # turning to 
getOption('timeout') # 60
options(timeout=100000000) # turning to 300 

# Baixando os de single cell  do GEO
wget("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz")

# Dezipando arquivo
untar("/media/hd/maycon/MESTRADO/Single_Cell/GSE176078/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz")

sc_matrix <- list.files("/media/hd/maycon/MESTRADO/Single_Cell/GSE176078", pattern = "txt.gz$", full = TRUE, recursive = TRUE)
sapply(sc_matrix, gunzip, overwrite = TRUE)

# OBS: as vezes temos que usa a seguinte função pra dezipar 
df <- read.table(gzfile("/media/hd/maycon/MESTRADO/InfiniumMethylation_probeprospect/HM450.hg38.manifest.gencode.v36.tsv.gz"),header = T)

# Carregando matriz de expressão
exp_matrix_GSE176078 <- read.delim("",header=TRUE)


# 120. dezipando arquivos
library(plyr)
zipF <- list.files(path = "/media/hd/maycon/DataScience/basket_dataset_kaggle", pattern = "*.zip", full.names = TRUE)
ldply(.data = zipF, .fun = unzip, exdir = "/media/hd/maycon/DataScience/basket_dataset_kaggle")

            
# 121. just a boxplot script/code 
ggplot(df[,], aes(x=celltype_minor,y=ss), fill = stemness) + geom_boxplot(outlier.color = NA)+ geom_jitter (alpha=0.01)  +
  xlab("") + ylab("") + theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), axis.title.y = element_text(colour="black", size = 12), axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1)) + ggtitle("") 


# 122. seurat tips
# coloring plots
## https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html


# 123. seurat plot edit
#ggplot + scale_color_viridis(discrete=TRUE) is a really efficient straightfoward method
library(viridis)
n_colors <- length(table(Sobj_sub_pdata@meta.data$CellType))
pal <- viridis(n = n_colors, option = "G", direction = -1)

FeatureScatter(object = Sobj_sub_pdata, feature1 = "NE_score", feature2 = "stemness")+    scale_fill_manual(values=c(pal), name="Legend") +
  ggtitle("Cells by patient") +
  ylab("Number of cells") + 
  xlab("Patient ID") + 
  theme(axis.title.x = element_text(colour="black", size = 18), 
        axis.title.y = element_text(colour="black", size = 18),
        legend.text = element_text(colour="black", size = 16),
        legend.title = element_text(colour="black", size = 18)) + 
  #geom_jitter(alpha=0.4) +
  theme(axis.text.x = element_text(angle = 45, size = 10 ,vjust=1,hjust=1),
        axis.text.y = element_text(angle = 0, size = 16)) +
  guides(fill = FALSE)  


# 124. upset plot
library(UpSetR)
listInput <- list(cancer_basal_DEGs = cancer_basal_DEGs$genes, 
                  cancer_cycling_DEGs = cancer_cycling_DEGs$genes, 
                  cycling_tcells_DEGs = cycling_tcells_DEGs$genes)

upset(fromList(listInput), order.by = "freq", nsets = 6)

# extracting elements from upsetplot 
# from: https://stackoverflow.com/questions/65027133/extract-intersection-list-from-upset-object
x <- upset(fromList(listInput), nsets = 6)
x$New_data[1:5, 1:5]
dim(x$New_data) #901 genes
x1 <- unlist(listInput, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #901 genes

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 6 ]
# Or 
Reduce(intersect, listInput)
# selecting features across specific groups

# Still need to be reviewed ,,,
# x <- upset(fromList(listInput), nsets = 14)
# df <- x$New_data %>% as.data.frame()
# df <- df[df$Melanocytic_probes %in% 1 & 
#            df$NSC_probes %in% 1 &
#            df$Haematopoietic_probes %in% 0 & 
#            df$glioblastoma_probes %in% 0 & 
#            df$gliomaIDH_probes %in% 0 &
#            df$embryonal_probes %in% 0 &
#            df$ependymal_probes %in% 0 &
#            df$glioneural_probes %in% 0 &
#            df$Mesenchymal_probes %in% 0 &
#            df$Nerve_probes %in% 0 &
#            df$OtherGlioma_probes %in% 0 &
#            df$Pineal_probes %in% 0 &
#            df$Plexus_probes %in% 0 &
#            df$Sella_probes %in% 0, ]
# 
# x1 <- unlist(listInput, use.names = FALSE)
# x1 <- x1[ !duplicated(x1) ]
# Mela_NSC_probes <- x1[as.numeric(rownames(df))] 




# 125. density plot 
# Density plot
dens_plot_1 = data.frame(stemness = Sobj_Ref@meta.data$all_stemness,
                         stemness_version = "all_stemness")

dens_plot_2 = data.frame(stemness = Sobj_Ref@meta.data$noNormal_stemness,
                         stemness_version = "noNormal_stemness")

dens_plot = plyr::rbind.fill(dens_plot_1, dens_plot_2)

ggplot(dens_plot, aes(x = stemness, fill = stemness_version)) + 
  geom_density(alpha = 0.5) + 
  theme_classic() 


# 126. quantile
quantile(pdata[pdata$celltype_major %in% "iPS",]$stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# 127. grep a row by a string
table(grepl("string",vector))
table(grepl("IPS",s[s$stemness >= 0.7055146,]$barcode))

# 128. merge more than one dataframe at time
pdata = list_ss %>% reduce(merge, by="barcode") #merging all daframes on a list by "barcode" column

# 129. Pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)

## NSC_hypo_island gene_list
load("/mnt/plummergrp/maycon/Masters/3_months_last/figures_thesis/fig_5/NSC_hypo_island_gene_list_to_GO.rda")
gene_list = data.frame(gene_list)
gene_list = gene_list[!duplicated(gene_list$genesUniq), ]

# Enrichment 
ego2 <- enrichGO(gene = gene_list,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2, drop=TRUE, showCategory=12, main = "")
clusterProfiler::dotplot(ego2) + ggtitle("NSC_hypo_island genes annotated to probes")


# 130. Combine same daframes in a list into one dataframe 
pdata_melted = do.call("rbind",list_x)


# 131. How to count files in sub directories
fils <- list.files("/mnt/plummergrp/maycon/Kim_thymus_fetal", pattern = ">0.1", full = TRUE, recursive = TRUE)
table(duplicated(fils)) # no duplicated file
dirs <- dirname(fils)
dirs <- sapply(dirs,digest::digest) 
freq_files = as.data.frame(table(dirs))
dirs_with_few_files = freq_files[!freq_files$Freq %in% 14, ]

df = data.frame(dirs)
df$dir_ext = names(dirs)
table(df[df$dirs %in% dirs_with_few_files$dirs, ]$dir_ext)
# One of those two are the directory it self  (/mnt/plummergrp/maycon/Kim_thymus_fetal) and the other one is the only directory with different number of files (/mnt/plummergrp/maycon/Kim_thymus_fetal/Enrichment_Reference_Adrenal) with only 11 files instead of 14.


# 132. Merging new meta.data to Seurat Object meta.data 



# i) AddMetaData 
# Obs: pData_sub need to have a column with all cell barcodes you want (same dimension)
#Sobj = AddMetaData(
#object = Sobj_Valid_Merged,
#metadata = pData_sub)
# Double check it !!


# ii) merge both meta.data and sign it to Seurat Object meta.deta - this code can merge or filter cells while it is merging meta.datas 
# First merge, then cbind
# meta_data = Sobj_Val_dsmp@meta.data
# meta_data_toadd = merge(meta_data, pData_sub) 
# 
# identical(meta_data_toadd$Barcode, meta_data$Barcode) 
# meta_data_toadd_final = cbind(meta_data, meta_data_toadd)
# 
# identical(Sobj_Val_dsmp@meta.data$Barcode,
#           meta_data_toadd_final$Barcode) #TRUE - - check if the meta.data we are about to sign over Sobj meta.data are in the same sequence
# 
# Sobj_Val_dsmp@meta.data = meta_data_toadd_final
# head(Sobj_Val_dsmp@meta.data) # okay 


# 133. How to FILTER cells on Seurat Object 
# ii) Subset Sobj colnames by cell barcodes 
# Sobj_subset = Sobj[, colnames(Sobj) %in% "cells barcode you want to keep"]


# 134. Convert to upper/lower case all letter in a string to 
toupper("your vector")
tolower("your vector")


# 135. stack barplot (requeries seurat object)
# https://bioconductor.org/packages/devel/bioc/vignettes/dittoSeq/inst/doc/dittoSeq.html
library(dittoSeq)
dittoBarPlot(Sobj_Ref, "celltype_major", group.by = "Patient_ID") #simple like that


# 136. intersect without need to check dimension 
all(rownames(metadata) %in% colnames(matrix_exp)) # TRUE


# 137. intersect elements in more than two vectors at once
# BECAREFUL !! To get the same intersection as in upsetplot we need to discount from "all intersections". 
all_intersects = Reduce(intersect, list(CRISPR_pos_final_table = CRISPR_pos_final_table$Description,
                                        CRISPR_neg_final_table = CRISPR_neg_final_table$Description,
                                        RIME_SKO_final_table = RIME_SKO_final_table$Description,
                                        RIME_DKO_final_table = RIME_DKO_final_table$Description
))

RIME_SKO_DKO_CRISPR_neg_intersects = Reduce(intersect, list(
  CRISPR_neg_final_table = CRISPR_neg_final_table$Description,
  RIME_SKO_final_table = RIME_SKO_final_table$Description,
  RIME_DKO_final_table = RIME_DKO_final_table$Description
))

RIME_SKO_DKO_CRISPR_neg_intersects = setdiff(RIME_SKO_DKO_CRISPR_neg_intersects,all_intersects) #removing all_intersects part of intersection 


# 138. save pheatmap .pdf
p = pheatmap(as.matrix(test), main = "Similarity Ontology Score")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "/mnt/plummergrp/maycon/Figures_to_paper/Figures/Figure_3/Similarity_Ontology_CancerBasal_Dataset_1_2.pdf")

# 139. Delete column by index (position)
df[, -(c(2, 5:10))]


# 140. update theme in ggplot globally
library(ggplot2); theme_set(theme_classic())


# 141. order groups/levels in a column to plot it
test$normal_complex_simple <- ordered(test$normal_complex_simple,
                                      levels = c("normal", "Complex carcinoma", "Simple carcinoma", "other subtypes"))

# 142. ANOVA/Tukey test to the plot visualization
# But can't get a good pvalue plot from all comparisons at once 
# Step 1: perform Tukey test 
tukey_result <- TukeyHSD(aov(stemness_index ~ normal_complex_simple, data = meta_binded))

# Step 2: Extract p-values from Tukey test
tukey_pvalues <- tukey_result[[1]][, 4]

# Step 3: Create the ggplot
library(ggplot2); theme_set(theme_classic())
my_plot <- ggplot(meta_binded[, ], aes(x=normal_complex_simple,y=stemness_index)) + 
  geom_boxplot(fill= jet.colors(colors), outlier.color = NA)+ 
  geom_jitter (alpha=0.2)  +
  xlab("Groups") + 
  ylab("Stemness") + 
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(angle = 45, vjust= 1, size = 12, hjust = 1)) +
  #facet_wrap(~ sample_char) +
  ggtitle("GSE119810 stemness - Normal and Complex/Simple Carcinoma") 


# Step 4: Add Tukey p-values to the ggplot using geom_signif()
my_plot_with_signif <- my_plot +
  geom_signif(comparisons = list(c("normal", "Complex carcinoma"),
                                 c("normal", "Simple carcinoma"),
                                 c("normal", "other subtypes"),
                                 c("Complex carcinoma", "Simple carcinoma"),
                                 c("Complex carcinoma", "other subtypes"),
                                 c("Simple carcinoma", "other subtypes")),
              map_signif_level = TRUE,
              y_position = 1) +
  labs(title = "GSE119810 stemness - Normal and Complex/Simple Carcinoma")

# Step 5: Print the plot
print(my_plot_with_signif)




# 143. Stack barplot (general dataset)
# df structure
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
# code 
library(ggplot2); theme_set(theme_classic())
library(viridis)
df_HR$HR_genes_per_f100 = df_HR$HR_genes_per * 10
ggplot(df_HR[df_HR$HR_OS >1, ], aes(reorder(x=CellType, +HR_genes_per_f100, sum), y=HR_genes_per_f100, fill=genes)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  xlab("Cell type") + 
  ylab("HR > 1 genes (%)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))


# 144. replace NA 
library(impute)
pca.data <- impute.knn(as.matrix(pca.data), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]



# 145. get gene promoter DNAmet probes annotation
# *************** Get probes at gene promoters ***************
# packages
library(GenomicRanges)
library(rtracklayer)
library(HelpersMG)
# Genconde > .gtf https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz


# NSC gene markers from Panglao database
panglao_merkers = read.table(file = '/media/hd/maycon/MESTRADO/Reta_final/3_months_last/PanglaoDB_markers_27_Mar_2020.tsv', sep = '\t', header = TRUE)
NSC_panglao_merkers = panglao_merkers[panglao_merkers$cell.type %in% 'Neural stem/precursor cells',]
NSC_genes = NSC_panglao_merkers$official.gene.symbol
length(NSC_genes) #58
NSC_genes[1] #SOX1

library(GenomicRanges)
library(rtracklayer)


gencode <- readGFF("/media/hd/maycon/MESTRADO/Reta_final/3_months_last/map_gene_promoter/gencode.v42.annotation.gtf.gz") #3413841      26 (this file has annotation of the entire genome)
gencode.s <- subset(gencode, type %in% "gene")

# map to 450K probes (check which probes are located in the promoter regions of your list of genes)
hm450.hg38.anno <- read.table(gzfile("/media/hd/maycon/MESTRADO/InfiniumMethylation_probeprospect/HM450.hg38.manifest.gencode.v36.tsv.gz"),header = T)
# PROBE ANOTATION 485578. 12
hm450.hg38.anno <- hm450.hg38.anno[!is.na(hm450.hg38.anno$CpG_chrm),] #remove probes with no chromosome annotation
rownames(hm450.hg38.anno) <- NULL
hm450.hg38.anno$order <- rownames(hm450.hg38.anno)
hm450.hg38.anno.gr <- makeGRangesFromDataFrame(hm450.hg38.anno, keep.extra.columns=F, start.field = "CpG_beg", end.field = "CpG_end", seqnames.field = "CpG_chrm", strand.field="probe_strand") #866857 

probes_mapped_list = list()

for(i in 1:length(NSC_genes)) {
  #62696 26  #select only "gene" category
  gencode.gene <- subset(gencode.s, gene_name %in% NSC_genes[i])
  rownames(gencode.gene) <- NULL
  gencode.gene$order <- rownames(gencode.gene)
  
  # get promoter location of these genes:
  gene.promoters <- promoters(makeGRangesFromDataFrame(gencode.gene, keep.extra.columns=T), upstream = 2000, downstream = 200) #get the window from the TSS, 2000 upstream and 200 downstream (this window can be changed, but this is well accepted in the literature)
  
  A <- findOverlaps(gene.promoters, hm450.hg38.anno.gr, ignore.strand=TRUE, select="all") #check the overlapping between probes in the promoters of your list of genes - based on the genomic location
  B.df <- as.data.frame(A) #9 2 (number of probes that overlap)
  
  # merge information of enhancers locations and overlapped probes
  aux <- as.data.frame(gene.promoters) #1
  aux2 <- merge(B.df, aux, by.x="queryHits", by.y="order", all.x=T, sort=F) #9 #column 'subjectHits' represents the index of the probe that overlap. In other words, it is the number of the row of the data frame
  # Identify probes that overlap gene promoters
  probe <- hm450.hg38.anno[aux2$subjectHits,] #probes that overlapp with gene promoters  #94
  probe.promoter <- merge(probe, aux2, by.x="order", by.y="subjectHits", all.x=T) #add info of list of genes
  probes_mapped_list[[i]] = probe.promoter
  print(paste0(NSC_genes[i]," promoter probes have been mapped."))
  
}

NSC_probes_mapped = do.call(rbind.data.frame, probes_mapped_list)
save(NSC_probes_mapped, file = "/media/hd/maycon/MESTRADO/Reta_final/3_months_last/map_gene_promoter/NSC_probes_mapped.rda")



# 147 - pheatmap 
load("/media/hd/maycon/Glioma_classifier/trainingdata.Rda")
my_sample_col <- data.frame(row.names = rownames(trainingdata), Subtypes = trainingdata$Supervised.DNA.Methylation.Cluster, DNAmet_cluster = trainingdata$IDH.specific.DNA.Methylation.Cluster) # patient ID in row names; any columns are for sample information

dim(RFtrain) #881 probes, 430 samples (DNA matrix)
library(pheatmap)
p = pheatmap(RFtrain, 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'IDHwt - IDHwt subtypes (881outof914)')

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "/media/hd/maycon/Glioma_classifier/heatmap_IDHwt_IDHwt_subtypes_(881outof914).pdf")




# 148. transforming comma-separated values in a cell into multiple rows
library(tidyverse)
df <- hypo_hyper_anno %>%
  separate_rows(distToTSS, sep = ";") %>% as.data.frame()


# 149. MY FUNCTION - Differential DNAmet 
# create a function --------------
compare_groups <- function(DNAmtx, metadata, group_var) {
  # Selecting groups to compare
  brain_tissue_meta <- metadata[metadata$Class.added %in% 'brain_CONTR', ]
  tumor_tissue_meta <- metadata[metadata$Class.added %in% group_var, ]
  
  # Removing NAs and preparing the DNAmtx input
  DNAmtx_brain <- DNAmtx[, colnames(DNAmtx) %in% brain_tissue_meta$barcode]
  DNAmtx_tumor <- DNAmtx[, colnames(DNAmtx) %in% tumor_tissue_meta$barcode]
  DiffMeth <- cbind(DNAmtx_brain, DNAmtx_tumor)
  #saveRDS(DiffMeth, file= '/media/hd/maycon/NSC_model_rethink/DNAmet_diff_analysis/CNS_specific_NSC_model/DiffMeth.rds')
  DiffMeth <- as.matrix(DiffMeth)
  #DNAmtx <- na.omit(DNAmtx)
  library(impute)
  DiffMeth <- impute.knn(as.matrix(DiffMeth), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
  values <- t(DiffMeth) #transpose the dataframe
  values <- data.frame(values)
  
  require(exactRankTests)
  require(parallel)
  w.p.values <- unlist(mclapply(values,
                                function(probe) {
                                  zz <- wilcox.exact(probe[1:  dim(DNAmtx_brain)[2]],
                                                     probe[c(dim(DNAmtx_brain)[2]+1) : dim(DiffMeth)[2]], exact=F) # excat = T nos da um resultado mais preciso, mas demora mais pra rodar o looping. Entao, exact = F é melhor custo-beneficio. . .
                                  z <- zz$p.value
                                  return(z)
                                }, mc.cores= 6))
  w.p.values.adj <- p.adjust(w.p.values, method = "BH")
  save(w.p.values, w.p.values.adj, file= paste0("/media/hd/maycon/NSC_model_rethink/DNAmet_diff_analysis/CNS_specific_NSC_model/", group_var,"_vs_brain_w.p.values_w.p.values.adj.rda"))
  
  DiffMeth[1:10, 1:10]
  DiffMeth <- as.data.frame(DiffMeth)
  DiffMeth$meanM1 <- apply(DiffMeth[, 1:  dim(DNAmtx_brain)[2]],1,mean,na.rm=T)  #group 1 (brain)
  DiffMeth$meanM2 <- apply(DiffMeth[, c(dim(DNAmtx_brain)[2]+1) : dim(DiffMeth)[2]],1,mean,na.rm=T) #group 2 (tumor)
  DiffMeth$DiffMean <- DiffMeth$meanM1 - DiffMeth$meanM2 #calculate the difference of mean (Diff mean)
  
  if (identical(rownames(DiffMeth), names(w.p.values.adj)) == TRUE) {
    # Do something
  } else {
    stop("rownames(DiffMeth) and names(w.p.values.adj are not identical (same order)")
  }
  
  DiffMeth$p.val <- w.p.values.adj
  DiffMeth$threshold <- "1" # threshold info will be used to define DMP and to color volcano plots
  b <- DiffMeth[DiffMeth$p.val < 0.01 & DiffMeth$DiffMean < -0.30,] #hypo NSC probes
  DiffMeth[rownames(b),"threshold"] <- "2"
  #c <- subset(DiffMeth, p.val < 0.001 & DiffMean > 0.6) #hyper in not.NPC
  c <- DiffMeth[DiffMeth$p.val < 0.01 & DiffMeth$DiffMean > 0.30,] #hyper NSC probes 
  DiffMeth[rownames(c),"threshold"] <- "3"
  print(table(DiffMeth$threshold))
  
  DiffMeth$probe_ID <- rownames(DiffMeth)
  hypo <- DiffMeth %>% 
    subset(threshold %in% c("2")) %>%
    select(probe_ID) 
  
  hyper <- DiffMeth %>% 
    subset(threshold %in% c("3")) %>%
    select(probe_ID) 
  
  save(hypo, hyper, file = paste0('/media/hd/maycon/NSC_model_rethink/DNAmet_diff_analysis/CNS_specific_NSC_model/', group_var, '_hyper_hypo_probes.rda'))
}

# In case you have multiple groups to compare
for(i in 3:14) { #2 had already been done; 1 is normal brain and should be skipped 
  compare_groups(DNAmtx =  Capper_CNS_DNAmtx,
                 metadata = Capper_CNS_metadata,
                 group_var = names(table(Capper_CNS_metadata$Class.added))[i])
  #rm(list = ls())
}


# 149. Changing theme plot on ggplot
library(ggplot2); theme_set(theme_classic())


# 150. bar plot side by side ('dodge')
# print(probe_orien_plot)
# n_probes DNAmet_orien    probeset
# 53        hyper          Mesenchymal_NSC_probes
# 60        hyper          NSC_probeset
# 57        hypo           Mesenchymal_NSC_probes
# 50        hypo           NSC_probeset

library(ggplot2); theme_set(theme_classic())
library(ggrepel)
ggplot(probe_orien_plot, aes(x = probeset, y = n_probes, fill = DNAmet_orien)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("hyper" = "red", "hypo" = "blue")) +
  labs(x = "Probeset", y = "Number of intersected probes", fill = "") +
  geom_text(aes(label = n_probes), vjust = -0.5) +
  #geom_text_repel(aes(label = n_probes), vjust = +2) +
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(angle = 45, vjust= 1, 
                                   size = 12, hjust = 1))


# 151. subset/filter vectors by vectors
listInput_two_by_two[[4]][listInput_two_by_two[[4]] %in% hyper$probe_ID]


# 152. subset/filter by a pattern on the column's elements of your interest
dmp[!grepl("NA", dmp$probeID), ] # regecting all element with a "NA" on it


# 154. Barplot side by side to compare DNA methylation orientation 
# Create a input for side by side barplot ('dodge')
#i <- 1
list_df_plot <- list()
for(i in 1:length(names(listInput_two_by_two))) {
  # 1. select tumor type
  tumor_type <- names(listInput_two_by_two)[i]
  tumor_type <- gsub('_NSC_probes', '', tumor_type)
  # 2. create a empty data frame for temporary storage 
  probe_orien_df <- data.frame(Probeset = 'NA',
                               DNAmet_hyper = 0,
                               DNAmet_hypo = 0)
  
  # 3. fill it up with tumor DNAmet orientation
  load(paste0("/media/hd/maycon/NSC_model_rethink/DNAmet_diff_analysis/CNS_specific_NSC_model/", tumor_type, "_hyper_hypo_probes.rda")) # it's all Differential Methyalted probes
  probe_orien_df[1, 'Probeset'] <- paste0(tumor_type, '_NSC_probes')
  probe_orien_df[1, 'DNAmet_hyper'] <- length(intersect(listInput_two_by_two[[i]], hyper$probe_ID))
  probe_orien_df[1, 'DNAmet_hypo'] <- length(intersect(listInput_two_by_two[[i]], hypo$probe_ID))
  
  # 4. fill it up with NSC DNAmet orientation
  load("/media/hd/maycon/NSC_model_rethink/DNAmet_diff_analysis/NSC_vs_Capper_Brain_trainset_DNAmetProbes.Rda") # it's all Differential Methyalted probes
  probe_orien_df[2, 'Probeset'] <- 'NSC_probeset'
  probe_orien_df[2, 'DNAmet_hyper'] <- length(intersect(listInput_two_by_two[[i]], hyper$probe_ID))
  probe_orien_df[2, 'DNAmet_hypo'] <- length(intersect(listInput_two_by_two[[i]], hypo$probe_ID))
  
  # 5. put it in a proper data structure for side by side barplots ('dodge')
  probe_orien_plot <- data.frame(n_probes = c(probe_orien_df$DNAmet_hyper, probe_orien_df$DNAmet_hypo),
                                 DNAmet_orien = c('hyper', 'hyper', 'hypo', 'hypo'),
                                 probeset = c(probe_orien_df$Probeset, probe_orien_df$Probeset),
                                 tumor = tumor_type)
  # 6. store it in a list of dataframes
  list_df_plot[[i]] <- probe_orien_plot
  
}

# Combine all dataframes from the list
df_plot <- do.call('rbind', list_df_plot)

# Plot by tumor saving one at a time
for(i in 1:length(names(listInput_two_by_two))) {
  # 1. select tumor type
  tumor_type <- names(listInput_two_by_two)[i]
  tumor_type <- gsub('_NSC_probes', '', tumor_type)
  
  # 2. plot and save this barplot to the given tumor
  library(ggplot2); theme_set(theme_classic())
  library(ggrepel)
  ggplot(df_plot[df_plot$tumor %in% tumor_type, ], 
         aes(x = probeset, y = n_probes, fill = DNAmet_orien)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("hyper" = "red", "hypo" = "blue")) +
    labs(x = "Probeset", y = "Number of intersected probes", fill = "") +
    geom_text(aes(label = n_probes), vjust = -0.5) +
    #geom_text_repel(aes(label = n_probes), vjust = +2) +
    theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
          axis.title.y = element_text(colour="black", size = 12), 
          axis.text.x = element_text(angle = 45, vjust= 1, size = 12, hjust = 1)) +
    ggtitle(paste0(tumor_type, " differential DNAmet compared to NSC's")) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file= paste0('/media/hd/maycon/NSC_model_rethink/DNAmet_diff_analysis/CNS_specific_NSC_model/', tumor_type,'_intersected_NSC_DNAmet_orien.pdf'), width = 6, height = 6, dpi = 300) 
}


# 155. primitive way of taking a look on a matrix distribution (feature x samples)
hist(matrix$sample1) 



# 156. calculate DiffMean on DNAmet data in the mode of "one-versus-all" or even "one-versus-another"
group_names <- names(table(metadata$Ceccarelli_six_subtypes))
calculate_diff_mean <- function(data, 
                                metadata, 
                                condition_one, 
                                threshold) {
  list_diffmean_dfs <- list()
  data <- as.data.frame(data)
  condition_all_but_one <- group_names[!group_names %in% condition_one]
  ###condition_all_but_one <- 'Mesenchymal-like' #TESTING only. 
  for(i in 1:length(condition_all_but_one)) {
    print(paste0("Processing group: ", condition_all_but_one[i]))
    
    try({ #it ignores any error that might prevent the code of going forward
      #I did it because of the absence of hyper or hypo probes were making the loop stop
      loop_data <- data
      condition_1_sampleID <- metadata[metadata$Ceccarelli_six_subtypes %in% condition_one, ]$Case
      condition_2_sampleID <- metadata[metadata$Ceccarelli_six_subtypes %in% condition_all_but_one[i], ]$Case
      
      print(paste0("Length of condition_1_sampleID: ", length(condition_1_sampleID)))
      print(paste0("Length of condition_2_sampleID: ", length(condition_2_sampleID)))
      
      loop_data$meanM1 <- apply(loop_data[, condition_1_sampleID], 1, mean, na.rm = TRUE)
      loop_data$meanM2 <- apply(loop_data[, condition_2_sampleID], 1, mean, na.rm = TRUE)
      loop_data$DiffMean <- loop_data$meanM1 - loop_data$meanM2
      loop_data$Comparison <- paste0(condition_one, '_', condition_all_but_one[i], '_', condition_one, '_', 'Orientation')
      
      
      
      loop_data$DNAmet_orientation <- NA
      #loop_data$DNAmet_orientation <- as.character(loop_data$DNAmet_orientation) #this time it's not necessary
      if(dim(loop_data[loop_data$DiffMean > threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean > threshold, ]$DNAmet_orientation <- 'hyper'
      } else {
        # do nothing
      }
      
      if(dim(loop_data[loop_data$DiffMean < -threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean < -threshold, ]$DNAmet_orientation <- 'hypo'
      } else {
        # do nothing
      }
      
      loop_data[loop_data$DNAmet_orientation %in% NA, ]$DNAmet_orientation <- 'not_diff'
      
      
      
      loop_data$probeID <- rownames(loop_data)
      
      list_diffmean_dfs[[i]] <- loop_data[, c('DiffMean', 'Comparison', 'DNAmet_orientation', 'probeID')]
      print(paste0(list_diffmean_dfs[[i]]$Comparison[1], ' has been stored into the list.'))
    },  silent = FALSE)
  }
  return(list_diffmean_dfs) 
}


list_diffmean_dfs <-calculate_diff_mean(data = DNAmtx,
                                        metadata = metadata,
                                        condition_one = group_names[1], #"Classic-like" 
                                        threshold = 0.1)

Classiclike_markers <- do.call('rbind', list_diffmean_dfs)
head(Classiclike_markers)
length(names(table(Classiclike_markers$Comparison))) # 5 (it should be 5) 
# Attention 2: it's okay to have more probes than the ~28k MVF because in this dataframe should be 5 comparisons. So it's the sum of 5*~28k probes at total.

list_diffmean_dfs <- calculate_diff_mean(data = DNAmtx,
                                         metadata = metadata,
                                         condition_one = group_names[2], #"Codel" 
                                         threshold = 0.3)

Codel_markers <- do.call('rbind', list_diffmean_dfs)
head(Codel_markers)
length(names(table(Codel_markers$Comparison))) # 5 (it should be 5)


#########GWAS##############
## R script to carry out a GWAS analysis
## kinship matrix used to account for population structure in the data
## input: Plink .raw and .map files + phenotype file
# run as Rscript --vanilla gwas.R genotype_file=path_to_genotypes snp_map=path_to_map phenotype_file=path_to_phenotypes trait=trait_name_in_phenotype_file trait_label=label_to_use_for_trait

#setwd("/Users/adri/Nextcloud/Documents/Rumigen/Scripts/epiGWAS/")
setwd("C:/Users/alopez.catalina/Nextcloud/Documents/Rumigen/Scripts/epiGWAS/")
#setwd("C:/Users/alopez.catalina/OneDrive - Universidad Politécnica de Madrid/Documentos/Rumigen/Scripts/epiGWAS")

## libraries

library("qqman")
library("gMatrix")
library("data.table")
library("factoextra")
library(stringr)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
txdb <- TxDb.Btaurus.UCSC.bosTau9.refGene
library(clusterProfiler)
library(viridis)
library(org.Bt.eg.db)
library(biomaRt)
library(fuzzyjoin)

## include common functions

source("scripts/gwas.r")
source("scripts/emma.r")

# Read imputed CpGs
#CpG_imputed = fread("/Users/adri/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/Misisng_imputation/imputed/Imputed_test_short.tsv")
CpG_imputed = fread("C:/Users/alopez.catalina/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/Misisng_imputation/imputed/Imputed_test_short.tsv")
#CpG_with_NA = fread("/Users/adri/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/Misisng_imputation/withNA/NA_removed_standard_BetaValues_1_and_2_envio.txt")
CpG_with_NA = fread("C:/Users/alopez.catalina/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/Misisng_imputation/withNA/NA_removed_standard_BetaValues_1_and_2_envio.txt")
rownames(CpG_imputed) = CpG_with_NA$Probe_ID
CpG_imputed_to_merge = as.data.frame(t(CpG_imputed))
colnames(CpG_imputed_to_merge) = rownames(CpG_imputed)
CpG_imputed_to_merge$MUESTRA = rownames(CpG_imputed_to_merge)

# Read and trim metadata

#namescoded1 = fread("/Users/adri/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/01_MUESTRAS/relacion.csv")
namescoded1 = fread("C:/Users/alopez.catalina/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/01_MUESTRAS/relacion.csv")
#namescoded2 = fread("/Users/adri/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/02_MUESTRAS/relacion.csv")
namescoded2 = fread("C:/Users/alopez.catalina/Nextcloud/Datos/RUMIGEN/EpiChip/RUMIGEN_SEPT_2024/02_MUESTRAS/relacion.csv")
namescoded = rbind(namescoded1,namescoded2)

#gestacion = fread("/Users/adri/Nextcloud/Documents/Rumigen/Metadata Terneras/metadatos_previo.csv")
gestacion = fread("C:/Users/alopez.catalina/Nextcloud/Documents/Rumigen/Metadata Terneras/metadatos_previo.csv")

gestacion$fnac <- as.Date(sub(" .*", "", gestacion$fnac), format = "%d/%m/%Y")
gestacion$fcontrol <- as.Date(sub(" .*", "", gestacion$fcontrol), format = "%d/%m/%Y")
gestacion$fpar <- as.Date(sub(" .*", "", gestacion$fpar), format = "%d/%m/%Y")
gestacion = gestacion %>% dplyr::select(-c(Envio,Observaciones,Ganaderia,Nombre,Fecha_envio,
                                    Teléfono, `Correo electrónico`,Registro, 
                                    `Tipo Muestra`, Chip, numeromuestra))


gestacion = gestacion %>%
  mutate(Dam = case_when(
    numpar == 1 & fpar == fnac & edad_ternera_al_parto == 0 ~ "PP",
    TRUE ~ "MP"
  ),
  tiempo_al_control = as.numeric(fnac - fcontrol))


subset_gestacion <- gestacion %>%
  filter(
    Dam == "PP" | 
    (Dam == "MP" & tiempo_al_control <= 290 & tiempo_al_control >= 270)
  ) %>%
  distinct() %>%  
  group_by(RegistroGenealogico) %>%
  filter(
    Dam == "PP" & tiempo_al_control == min(tiempo_al_control) | 
    (Dam == "MP" & row_number() == 1)
  ) %>%
  ungroup() %>%
  group_by(RegistroGenealogico) %>%
  slice_max(bhb, with_ties = FALSE) %>%  # Keep the row with the highest bhb
  ungroup()


designmatrix = merge(subset_gestacion,namescoded, by.x="cib",by.y="IDENTIFICACION")
designmatrix$MUESTRA = paste0("Sample_",designmatrix$MUESTRA)


# Merge metadata and methylation information
ewasfile = merge(designmatrix,CpG_imputed_to_merge,by="MUESTRA")
ewasfile = ewasfile %>% dplyr::select(-c(MUESTRA,cib,RegistroGenealogico_madre,fnac,fpar,fcontrol,
  numpar,leche,grasa,proteina,rcs,bhb,
  year,edad_ternera_al_parto,tiempo_al_control,`TIPO DE MUESTRA`,
  CHIP,POSICION))

ewasfile <- ewasfile %>%
  mutate(Dam = factor(Dam, 
                                levels = c("PP", "MP"),
                                labels = c(0, 1)))
ewasfile = ewasfile[-which(duplicated(ewasfile$RegistroGenealogico)),]

###################################
## read arguments 
###################################
#dataset = basename(CpG_file)
## READING DATA


#snpMatrix <- fread("/Users/adri/Nextcloud/Documents/Rumigen/Metadata\ Terneras/GenotiposIMP_Newmap202409_RUMI/Genos_Clean.txt")
snpMatrix = fread("C:/Users/alopez.catalina/Nextcloud/Documents/Rumigen/Metadata\ Terneras/GenotiposIMP_Newmap202409_RUMI/Genos_Clean.txt")

print(paste(nrow(snpMatrix),"records read from the genotype file",sep=" "))

snp_vec = snpMatrix$V1
design_vec <- ewasfile$RegistroGenealogico

common_pos = intersect(snp_vec,design_vec)
  


## kinship matrix
print("Calculating the kinship matrix")
X = snpMatrix[which(snpMatrix$V1 %in% common_pos),-1]
#X <- as.data.frame(snpMatrix[,-1])
#colnames(X) <- gsub("\\_[A-Z]{1}$","",colnames(X))
rownames(X) <- unname(unlist(snpMatrix[which(snpMatrix$V1 %in% common_pos),1]))

#vec <- rownames(X) %in% CpG_file$RegistroGenealogico
#X = X[vec, ]

#setwd("/Users/adri/Nextcloud/Datos/RUMIGEN/EpiChip/EWAS")
setwd("C:/Users/alopez.catalina/Nextcloud//Datos/RUMIGEN/EpiChip/EWAS")
CpG_INFO<-data.frame(str_extract(colnames(X), "(?<=chr)[A-Za-z0-9]+"),colnames(X),0,str_extract(colnames(X), "(?<=_)[0-9]+"))
names(CpG_INFO) <- c("Chr","SNP","cM","Pos")
CpG_INFO$Chr[CpG_INFO$Chr=="MT"]<-30
CpG_INFO$Chr[CpG_INFO$Chr=="X"]<-31
CpG_INFO$Chr<-as.numeric(CpG_INFO$Chr)
CpG_INFO$Pos<-as.numeric(CpG_INFO$Pos)

SNP_INFO<-CpG_INFO

K <- gVanRaden.2(X)
fwrite(K,"kinship_matrix.tsv",col.names=T,row.names=T,sep="\t")

##Check for pop stratification visually
res.pca<-prcomp(K, scale = FALSE)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             label = "none"
)

#vec <- colnames(K) %in% CpG_file$RegistroGenealogico
#K <- K[vec,vec]

############ Nombre
dataset = "ewasimputed"

print("writing out the kinship matrix ...")
fname = paste(dataset,".kinship",sep="")
write.table(K, file=fname, quote = FALSE, row.names = FALSE)

print("producing the heatmap kinship matrix ...")
pdf(paste(dataset,"_kinship_heatmap",".pdf",sep=""))
heatmap(K,col=rev(heat.colors(75)))
dev.off()

###################
## Running the epiGWAS
###################
filtered_CpG_file = ewasfile[which(ewasfile$RegistroGenealogico %in% common_pos),]
CpG_matrix = as.data.frame(filtered_CpG_file[,-c(1:2)])
rownames(CpG_matrix) = ewasfile[which(ewasfile$RegistroGenealogico %in% common_pos),]$RegistroGenealogico
Y <- as.data.frame(as.vector(as.numeric(as.factor(filtered_CpG_file$Dam)) - 1 ))
rownames(Y) <- filtered_CpG_file$RegistroGenealogico

#Xtrim = X[which(rownames(X) %in% rownames(Y)),]
#Ytrim = as.data.frame(Y[which(rownames(Y) %in% rownames(K)),])
#Ktrim = as.data.frame(K[which(rownames(K) %in% rownames(Y)),])

#manifest = fread("/Users/adri/Nextcloud/Documents/Rumigen/RUMIGENEpichip_ Manifest Files20240311095639/Rumigen_Epichip_SS_20042758X389925_A1.csv")
manifest = fread("C:/Users/alopez.catalina/Nextcloud/Documents/Rumigen/RUMIGENEpichip_ Manifest Files20240311095639/Rumigen_Epichip_SS_20042758X389925_A1.csv")
SNP_INFO = manifest %>% dplyr::select(Probe_ID,CHR,MAPINFO)
SNP_INFO = SNP_INFO[which(SNP_INFO$Probe_ID %in% colnames(CpG_matrix)),]
colnames(SNP_INFO) = c("SNP","Chr","Pos")
SNP_INFO$Chr[SNP_INFO$Chr=="MT"]<-30
SNP_INFO$Chr[SNP_INFO$Chr=="X"]<-31
SNP_INFO$Chr<-as.numeric(SNP_INFO$Chr)
SNP_INFO$Pos<-as.numeric(SNP_INFO$Pos)

print("Running epiGWAS ...")
res <- amm_gwas(Y = Y, X = CpG_matrix, K = K, use.SNP_INFO = T)




###########
### RESULTS
###########


print("writing out results and figures ...")
gwasResults <- res[,c("SNP","Chr","Pos","Pval")]



# add FDR value - multiple testing correction
gwasResults$FDR <- p.adjust(gwasResults$Pval, method="fdr")
gwasResults$bonferroni <- p.adjust(gwasResults$Pval, method = "bonferroni")

names(gwasResults) <- c("SNP","CHR","BP","P","Padj_FDR", "Padj_bonf")

#Write results in a file
fname <- paste(dataset,"GWAS.results", sep="_")
fwrite(x = gwasResults, file = fname)
#Write results with FDR<0.05 in a file
gwasResults_sig_fdr <- gwasResults[gwasResults$Padj_bonf < 0.05,]
fname_sig <- paste(dataset,"GWAS.resultsBonferroni", sep="_")
fwrite(x = gwasResults_sig_fdr, file = fname_sig)


#Manhattan plot
gwasResults$CHR[gwasResults$CHR=="MT"]<-30
gwasResults$CHR[gwasResults$CHR=="X"]<-31
gwasResults$CHR<-as.numeric(gwasResults$CHR)
gwasResults$Pos<-as.numeric(gwasResults$BP)
gwasResults = gwasResults[!grepl("^NKL", gwasResults$CHR), ]
gwasResults <- gwasResults[!is.na(gwasResults$CHR) & grepl("^[0-9]+$", gwasResults$CHR), ]

fwrite(gwasResults,"gwasResults.tsv",sep="\t")
png(paste(dataset,"manhattan.png",sep="_"), width = 1200, height = 600, res = 100)
manhattan(gwasResults, chr = "CHR", bp = "BP", p = "P", snp = "SNP",col = c("blue4", "orange3"))
dev.off()

#qq-plot for pvalues
png(paste(dataset,"qqplot.png",sep="_"), width = 600, height = 600)
qq(gwasResults$P)
dev.off()




bonferronithreshold = unname(unlist(gwasResults[which.min(abs(gwasResults$Padj_bonf - 0.05)),"P"]))
fdrthreshold = unname(unlist(gwasResults[which.min(abs(gwasResults$Padj_FDR - 0.05)),"P"]))

gwasResults = gwasResults %>% 
  filter(CHR!=0)


bonferronithreshold = unname(unlist(gwasResults[which.min(abs(gwasResults$Padj_bonf - 0.05)),"P"]))
fdrthreshold = unname(unlist(gwasResults[which.min(abs(gwasResults$Padj_FDR - 0.05)),"P"]))

gwasResults = gwasResults %>% 
  filter(CHR!=0)

don <- gwasResults %>% 
  #mutate(CHR = as.character(CHR)) %>%
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  mutate(is_highlight=ifelse(Padj_bonf <= 0.05, "yes", "no"))

axisdf = don %>%
  group_by(CHR) %>%
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )



# Assuming 'don' is your data frame and it has columns: BPcum, P, CHR, is_highlight

# Create the chrom_color column
don <- don %>%
  mutate(chrom_color = case_when(
    is_highlight == "yes" & CHR %% 2 == 1 ~ "#597445",    # Highlighted and odd CHR
    is_highlight == "yes" ~ "#e68e3c",                   # Highlighted and even CHR or character
    CHR == 31 ~ "#FFE6A9",                       # CHR is a character
    CHR %% 2 == 1 ~ "#729762",                           # Odd CHR
    TRUE ~ "#FFE6A9"                                     # Default (even CHR)
  ),
  alpha_case = case_when(
    is_highlight == "yes" & CHR %% 2 == 1 ~ 1,    # Highlighted and odd CHR
    is_highlight == "yes" ~ 1,                   # Highlighted and even CHR or character
    CHR == 31 ~ 0.3,                       # CHR is a character
    CHR %% 2 == 1 ~ 0.3,                           # Odd CHR
    TRUE ~ 0.3),
  size_case = case_when(
    is_highlight == "yes" & CHR %% 2 == 1 ~ 2.5,    # Highlighted and odd CHR
    is_highlight == "yes" ~ 2.5,                   # Highlighted and even CHR or character
    CHR == 31 ~ 2.3,                       # CHR is a character
    CHR %% 2 == 1 ~ 2.3,                           # Odd CHR
    TRUE ~ 2.3 ),
    shape_case = case_when(
      is_highlight == "yes" ~ 18,
      TRUE ~ 16             # Highlighted and even CHR or character
      )
    )

# Ensure CHR is treated as a factor
don$CHR <- as.factor(don$CHR)

# Create the plot
ggplot(don, aes(x = BPcum, y = -log10(P))) +
  # Show all points with custom colors
  geom_point(aes(color = chrom_color, alpha = alpha_case, shape = shape_case, size = size_case)) +
  scale_color_identity() +  # Use the colors directly from the chrom_color column
  scale_shape_identity() +
  scale_size_identity() +
  # Custom X axis:
  scale_x_continuous(
    labels = c(1:29, "MT"),
    breaks = axisdf$center
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove space between plot area and x axis

  # Custom the theme:
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(size=25),
    axis.text.x = element_text(size = 10)
  ) +
  geom_hline(yintercept = -log10(bonferronithreshold), linetype = "dotdash", color = "black") +
  geom_hline(yintercept = -log10(fdrthreshold), linetype = "dotted", color = "black") +
  xlab("Chromosome") +
  ylim(0, 18)



print("#########")
print("## END ##")
print("#########")

gwasResultsp05 = gwasResults %>% filter(Padj_bonf < 0.05)

ewas_results_chipseeker = data.frame(
  Chr = paste0("chr", gwasResultsp05$CHR),
  Start = gwasResultsp05$Pos - 1,
  End = gwasResultsp05$Pos + 1,
  r = paste0("r", 1:nrow(gwasResultsp05))
)

write.table(ewas_results_chipseeker, "Ewas_ChipSeeker.tsv",quote=F, col.names = F,row.names=F,sep="\t")

files = "Ewas_ChipSeeker.tsv"
peak <- readPeakFile(files)
covplot(peak)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix)

plotPeakProf2(peak = peak, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb,ignore_strand = T)


peakAnno <- annotatePeak(files,
   tssRegion=c(-3000, 3000),
   annoDb="org.Bt.eg.db",
   TxDb=txdb)

plotAnnoBar(peakAnno) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  ggtitle("") +
  labs(fill="")

genes = as.data.frame(peakAnno)$geneId
names(genes) = sub("_", "\n", names(genes))
write.table(unique(as.data.frame(peakAnno)$SYMBOL), "Genes_in_peak_ewas.tsv",sep="\t",quote=F,col.names=T,row.names=F)

compGO<- compareCluster(geneCluster   = genes,
                        OrgDb = org.Bt.eg.db,
                        fun           = "enrichGO",
                        ont = "ALL",
                        pool = T,
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        readable = T)
                        
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")

# Hacer el enrichment con Webgestalt


# Script Oscar

ensembl_dataset = "btaurus_gene_ensembl" # ARS UCD 1.3
window = 50000 # number of bases to search upstream and downstream the SNP position
ensembl = biomaRt::useEnsembl(biomart="ensembl",dataset=ensembl_dataset)

results = gwasResultsp05
rownames(results) <- results$SNP
results$CHR[results$CHR=="30"]<-"MT"
results$CHR[results$CHR=="31"]<-"X"
genes = list()

for (snp_name in rownames(results)) {
  snp = results[snp_name,]
  genes[[snp_name]] = biomaRt::getBM(c('ensembl_gene_id',
                                       'entrezgene_id',
                                       'external_gene_name',
                                       'start_position',
                                       'end_position',
                                       'uniprotsptrembl',
                                       'uniprotswissprot'),  
                                     filters = c("chromosome_name","start","end"),
                                     values=list(snp$CHR,snp$BP-window,snp$BP+window),
                                     mart=ensembl)
}

gwas_genes <- ldply(genes, function(x) {
  rbind.data.frame(x)
})

gwas_genes <- gwas_genes[!is.na(gwas_genes$external_gene_name) & gwas_genes$external_gene_name != "",]

gwas_genes <- do.call(rbind,genes)$external_gene_name

View(gwas_genes)
gwas_genes <-unique(gwas_genes)
gwas_genes<-gwas_genes[!is.na(gwas_genes) & gwas_genes!=""]

## write out file
write.table(gwas_genes,file = "epigwas_genes.tsv", sep = "\t", col.names = FALSE,row.names = FALSE, quote=FALSE)


# Test ensembl

ensembl = useEnsembl(biomart="ensembl", dataset="btaurus_gene_ensembl")
View(listFilters(ensembl))
##T

bos_genes <- getBM(attributes=c('ensembl_gene_id',
'ensembl_transcript_id','hgnc_symbol',"external_gene_name",'chromosome_name','start_position','end_position'),  mart = ensembl)

bos_genes_filtered <- bos_genes %>%
  filter(str_starts(external_gene_name, "DGAT")) %>%
  distinct()

bos_genes_dgat1 = bos_genes_filtered[1,] %>%
  dplyr::select(external_gene_name,chromosome_name,start_position,end_position)

gwasResultsp05$CHR[gwasResultsp05$CHR=="30"]<-"MT"
gwasResultsp05$CHR[gwasResultsp05$CHR=="31"]<-"X"


CpGtraspose = as.data.frame(t(CpG_imputed_to_merge))
CpGtraspose$Probe_ID = rownames(CpGtraspose)

manifest_clean = manifest %>% dplyr::select(Probe_ID,CHR,MAPINFO)
CpGmapped = merge(CpGtraspose,manifest_clean,by="Probe_ID")

CpGmapped_dgat = CpGmapped %>% 
  dplyr::filter(CHR %in% bos_genes_dgat1$chromosome_name)

window = 10000
filtered_data <- CpGmapped_dgat %>%
  left_join(bos_genes_dgat1, by = c("CHR" = "chromosome_name"))  %>%
  filter(MAPINFO >= start_position-window & MAPINFO <= end_position+window) 

### Estimate average methylation
for(i in 1:nrow(filtered_data)){
  assign(filtered_data$Probe_ID[i],filtered_data[i,])
  write.table(t(filtered_data[i,]),file = paste0("blupf90/DGAT1_",filtered_data$Probe_ID[i],".tsv"),
              quote=F,col.names = F,row.names=F,sep="\t")
}
CpGtraspose_means = CpGtraspose %>%
  rownames_to_column(var = "rowname") %>%
  filter(!(rowname == "MUESTRA" | grepl("^ctl", rowname))) %>%
  column_to_rownames(var = "rowname") %>%
  mutate(across(-ncol(.), as.numeric))

meth_means = apply(CpGtraspose_means[,-ncol(CpGtraspose_means)],2,mean)


long_data <- filtered_data %>%
  pivot_longer(
    cols = -c(Probe_ID, start_position, end_position, CHR, MAPINFO, external_gene_name),
    names_to = "Sample_ID",
    values_to = "Value"
  ) %>%
  dplyr::select(Probe_ID,Value)
ggplot(long_data, aes(y = Value, colour = Probe_ID)) +
  geom_density() +
  labs(title = "Density Plot of SNP Values",
      x = "Value",
      y = "Density",
      colour = "SNP (Probe_ID)") +
  theme_minimal()

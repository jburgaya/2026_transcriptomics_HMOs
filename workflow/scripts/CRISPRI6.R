library(readxl) # reading Excel files
library(DESeq2) # differential enrichment analysis
library(ggplot2) # easy pretty plots
library(ggrepel) # adding labels to plots
library(ggtext) # customizing label typesetting on plots
library(scales) # making pretty units for numbers
library(RColorBrewer) # for color palettes
library(ggbeeswarm) # making beeswarm plots
library(reshape2) # easy table format conversion
library(patchwork) # combine multiple plots into one frame
library(pheatmap) # quick heatmaps
library(viridis) # colorblind-friendly continuous color palette
library(plotly) # make interactive version of plots
library(writexl) # write to Excel files
library(ComplexHeatmap)
library(limma)


setwd("~/Project_Crispr_Kathi")
path.out <-  "~/Project_Crispr_Kathi/6batch"

#PATH_counts <- "~/Project_Crispr_Kathi/4batch/compiled.csv" #the sgRNA counts
#PATH_coldata <- "~/Project_Crispr_Kathi/20231004_kb_sample_data_2.xlsx" #the coldata, with information on each column, or biological replicate sample, in the sgRNA counts file
#PATH_sgRNAinfo <- "~/Project_Crispr_Kathi/targets_operon_manual-curation_20220401.xlsx" #a table with for each sgRNA (per row) information including at least the target gene names

PATH_counts <- "~/Project_Crispr_Kathi/6batch/compiled_6th.csv" #the sgRNA counts
PATH_coldata <- "~/Project_Crispr_Kathi/6batch/metadata_6th.xlsx" #the coldata, with information on each column, or biological replicate sample, in the sgRNA counts file
PATH_sgRNAinfo <- "~/Project_Crispr_Kathi/targets_operon_manual-curation_20220401.xlsx" #a table with for each sgRNA (per row) information including at least the target gene names


#The hypothesis test parameters / thresholds to use:

thr_lfc <- 1 # default 0: any change
thr_alpha <- 0.05 # default 0.1

raw <- read.csv(file = PATH_counts, 
                check.names = FALSE, 
                sep = ",")

temp <- as.data.frame(t(raw))
colnames(temp) <- temp["#Feature",]
temp <- temp[-1,]
temp <- apply(temp,c(1,2), function(x) as.numeric(x))
temp <- as.data.frame(temp)
temp$group <- sub("_L1.*", "\\1", sub("kbr6th_*", "\\1", rownames(temp)))

temp2 <- aggregate(x = temp[ ,colnames(temp) != "group"],
                   by = list(temp$group),
                   FUN = sum)
rownames(temp2) <- temp2$Group.1
temp2 <- temp2[,-1]

cts <- as.matrix(t(temp2))


#cts <- as.matrix(raw[, -1]) # take off the 1st column with the sgRNA names
#rownames(cts) <- raw[, 1] # set sgRNA names as rownames instead


coldata_raw <- read_excel(PATH_coldata)
coldata <- as.data.frame(coldata_raw) # force class data frame
#coldata$generations_induced_scaled <- scale(coldata$generations_induced)
#cts <- cts[,as.character(coldata$ID)]
rownames(coldata) <- colnames(cts) # check if samples in counts and coldata are matched

###Important: make sure the rownames and the sample ID match!

#Remove samples with Glucose
#cts <- cts[, -c(57:80)]
#coldata <- coldata[-c(57:80),]

# find correct order cts columns to match coldata rows:
#cts_ID <- sub("kb_", "\\1", sub("_L1_.*", "\\1", colnames(cts)))
cts_ID <- colnames(cts)
# store reordered cts:
cts_sort <- cts[, match(coldata$ID, cts_ID)]
# reassign names to check if match between cts and coldata:
rownames(coldata) <- colnames(cts_sort)

#Lastly, make truly separate variables for strain and inducer:
coldata$induced <- ifelse(coldata$inducer == "no", "no", "yes")


#subset to samples with biological replicates
comb_rep <- interaction(coldata$media, 
                        coldata$`condition`, 
                        coldata$strain, coldata$inducer)
true_rep <- sapply(comb_rep, function(x){sum(x == comb_rep)})
coldata_sel <- coldata[true_rep > 1, ]


cts_sel <- cts_sort[, true_rep > 1]


#target list
targets_raw <- read_excel(PATH_sgRNAinfo)
targets <- targets_raw[match(rownames(cts), targets_raw$sgRNA_name), ]

# check how many sgRNAs of the counts table are in the targets table:
table(rownames(cts) %in% targets$sgRNA_name)
# vice versa
table(targets$sgRNA_name %in% rownames(cts))

write.table(cts_sel, file=paste(path.out, "/sgRNA_counts", ".csv", sep=""), 
            quote=F, sep = ",", row.names = T, col.names = T)



#Define the model
dds <- DESeqDataSetFromMatrix(countData = cts_sel, 
                              colData = coldata_sel, 
                              design = ~ condition * induced)


#Compare induced with uninduced, so make sure latter is the reference group:
dds$induced <- relevel(dds$induced, ref = "no")

#Normalize
# standard DESeq2 normalized counts (akin to transcripts per million):
dds <- estimateSizeFactors(dds)
cts_norm <- counts(dds, normalized = TRUE)
# regularized log:
rld <- rlog(dds)

####################
#Quality control

#heatmap
p1 <- pheatmap(t(log10(1+cts_norm)), 
               show_colnames = FALSE, 
               labels_row = paste0("inducer: ", colData(dds)$inducer, 
                                   ", condition: ", colData(dds)$condition, 
                                   ", strain:", colData(dds)$strain), 
               main = "log10(normalized count + 1)", 
               col = viridis(100, option = "magma"),
               legend=F)



pdf(paste(path.out,"/Heatmap_all_samples", ".pdf", sep=""), width=20, height=12)
plot(p1)
dev.off()


png(paste(path.out,"/Heatmap_all_samples", ".png", sep=""), width=1600, height=750)
plot(p1)
dev.off()

#PCA
pcaData <- plotPCA(rld, intgroup = c("condition", "strain", "inducer", "ID"), 
                   returnData = TRUE, ntop = nrow(rld))
percentVar <- round(100 * attr(pcaData, "percentVar"))

# plot
ggplot(pcaData, 
       aes(x = PC1, 
           y = PC2, 
           fill = condition, 
           shape = interaction(inducer, strain))) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(size = 3) +
  # scale_color_brewer(palette = "Dark2") + 
  scale_shape_manual(values = 21:25) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = paste0("PC1: ", format(percentVar[1], digits = 3), "% variance"), 
       y = paste0("PC2: ", format(percentVar[2], digits = 3), "% variance"), 
       fill = NULL, 
       shape = NULL) +
  coord_fixed() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        legend.direction = "vertical",
        legend.box = "horizontal")


p2 <- ggplot(pcaData, 
             aes(x = PC1, 
                 y = PC2, 
                 shape = interaction(inducer, strain), 
                 fill = condition)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(size = 4) +
  # scale_color_brewer(palette = "Dark2") + 
  scale_shape_manual(values = 21:25) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = paste0("PC1: ", format(percentVar[1], digits = 3), "% variance"), 
       y = paste0("PC2: ", format(percentVar[2], digits = 3), "% variance"), 
       fill = NULL, 
       shape = NULL) +
  coord_fixed() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        legend.direction = "vertical",
        legend.box = "horizontal")



pdf(paste(path.out,"/PCA_all_samples", ".pdf", sep=""), width=15, height=7)
plot(p2)
dev.off()


png(paste(path.out,"/PCA_all_samples", ".png", sep=""), width=600, height=550)
plot(p2)
dev.off()


#PCA D39
cts1 <- cts_sel[, c(1:72)]
coldata1 <- coldata_sel[c(1:72),]
#Define the model
dds2 <- DESeqDataSetFromMatrix(countData = cts1, 
                               colData = coldata1, 
                               design = ~ condition * induced)

dds2$induced <- relevel(dds2$induced, ref = "no")

dds2 <- estimateSizeFactors(dds2)
cts_norm <- counts(dds2, normalized = TRUE)
rld2 <- rlog(dds2)


pcaData <- plotPCA(rld2, intgroup = c("condition", "inducer", "ID"), 
                   returnData = TRUE, ntop = nrow(rld))
percentVar <- round(100 * attr(pcaData, "percentVar"))


p2 <- ggplot(pcaData, 
             aes(x = PC1, 
                 y = PC2, 
                 shape = inducer, 
                 fill = condition)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(size = 4) +
  # scale_color_brewer(palette = "Dark2") + 
  scale_shape_manual(values = 21:25) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = paste0("PC1: ", format(percentVar[1], digits = 3), "% variance"), 
       y = paste0("PC2: ", format(percentVar[2], digits = 3), "% variance"), 
       fill = NULL, 
       shape = NULL) +
  coord_fixed() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        legend.direction = "vertical",
        legend.box = "horizontal")

pdf(paste(path.out,"/PCA_D39", ".pdf", sep=""), width=15, height=7)
plot(p2)
dev.off()



#PCA not induced only
cts1 <- cts_sel[, c(1:4,9:12,17:20,25:28,33:36,41:44,49:52,57:60,65:68,73:76,81:84,89:92,97:100)]
coldata1 <- coldata_sel[c(1:4,9:12,17:20,25:28,33:36,41:44,49:52,57:60,65:68,73:76,81:84,89:92,97:100),]
#Define the model
dds2 <- DESeqDataSetFromMatrix(countData = cts1, 
                               colData = coldata1, 
                               design = ~ condition)

dds2 <- estimateSizeFactors(dds2)
cts_norm <- counts(dds2, normalized = TRUE)
rld2 <- rlog(dds2)


pcaData <- plotPCA(rld2, intgroup = c("condition", "strain","ID"), 
                   returnData = TRUE, ntop = nrow(rld))
percentVar <- round(100 * attr(pcaData, "percentVar"))


p2 <- ggplot(pcaData, 
             aes(x = PC1, 
                 y = PC2, 
                 shape = strain, 
                 fill = condition)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(size = 4) +
  # scale_color_brewer(palette = "Dark2") + 
  scale_shape_manual(values = 21:25) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = paste0("PC1: ", format(percentVar[1], digits = 3), "% variance"), 
       y = paste0("PC2: ", format(percentVar[2], digits = 3), "% variance"), 
       fill = NULL, 
       shape = NULL) +
  coord_fixed() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        legend.direction = "vertical",
        legend.box = "horizontal")

pdf(paste(path.out,"/PCA_notInduced", ".pdf", sep=""), width=15, height=7)
plot(p2)
dev.off()


#########################
#Differential enrichment

#25mM Glucose (pre-growth Glucose)
#25mM Lactose (pre-growth Glucose)
#25mM Lactose (pre-growth Galactose)
#25mM Galactose (pre-growth Galactose)
#25mM 6SL (pre-growth Galactose)
#25mM 3SL (pre-growth Galactose)
#25mM LNT (pre-growth Galactose)
#100mM 3FL (pre-growth Galactose)


temp_col <- coldata_sel[,-c(2,6,8)]
temp_col$Col <- 1:104

###EXPERIMENTS
#For R6
#1 experiment:(pre-growth Glucose) Glucose vs.Lactose
   exp1_col <- temp_col[temp_col$strain=="R6" & 
               temp_col$pregrowth=="Glucose" & 
               temp_col$media %in% c("Glucose","Lactose"),"Col"]
   
#2 experiment:(pre-growth Galactose) Galactose vs Lactose
   exp2_col <- temp_col[temp_col$strain=="R6" & 
                          temp_col$pregrowth=="Galactose" & 
                          temp_col$media %in% c("Galactose","Lactose"),"Col"]
   
#For D39
#3 experiment:(pre-growth Glucose) Glucose vs Lactose 
   exp3_col <- temp_col[temp_col$strain=="D39" & 
                          temp_col$pregrowth=="Glucose" & 
                          temp_col$media %in% c("Glucose","Lactose"),"Col"]
   
#4 experiment:(pre-growth Galactose) Galactose (OD 0.4 or ID 457-464) vs Lactose
   exp4_col <- temp_col[temp_col$strain=="D39" & 
                          temp_col$pregrowth=="Galactose" & 
                          temp_col$OD==0.40 &
                          temp_col$media %in% c("Galactose","Lactose"),"Col"]
   
#5 experiment:Galactose (OD 0.4 or ID 457-464) vs 6SL
   exp5_col <- temp_col[temp_col$strain=="D39" & 
                          temp_col$pregrowth=="Galactose" & 
                          temp_col$OD==0.40 &
                          temp_col$media %in% c("Galactose","6SL"),"Col"]
   
#6 experiment:Galactose (OD 0.4 or ID 457-464) vs 3SL
   exp6_col <- temp_col[temp_col$strain=="D39" & 
                          temp_col$OD==0.40 & 
                          temp_col$media %in% c("Galactose","3SL"),"Col"]
   
#7 experiment:Galactose (OD 0.12 or ID 489-496) vs LNT
   exp7_col <- temp_col[temp_col$strain=="D39" & 
                          temp_col$OD==0.12 & 
                          temp_col$media %in% c("Galactose","LNT"),"Col"]
   
#8 experiment:Galactose (OD 0.12 or ID 489-496) vs 3FL
   exp8_col <- temp_col[temp_col$strain=="D39" & 
                          temp_col$OD==0.12 & 
                          temp_col$media %in% c("Galactose","3FL"),"Col"]


######XXXXXXX#############################################
#X experiment

exp_list <- c("exp1_col","exp2_col","exp3_col","exp4_col",
              "exp5_col","exp6_col","exp7_col","exp8_col")

for (i in 1:8){
  exp_col <- eval(parse(text = exp_list[i]))
  cts1 <- cts_sel[, exp_col]
  coldata1 <- coldata_sel[exp_col,]
   
  #Define the model
  dds <- DESeqDataSetFromMatrix(countData = cts1, 
                                colData = coldata1, 
                                design = ~ condition * induced)
     
  #Compare induced with uninduced, so make sure latter is the reference group:
  dds$induced <- relevel(dds$induced, ref = "no")
     
  #Normalize
  # standard DESeq2 normalized counts (akin to transcripts per million):
  dds <- estimateSizeFactors(dds)
  cts_norm <- counts(dds, normalized = TRUE)
  #regularized log:
  rld <- rlog(dds)
  
  dds <- DESeq(dds)

  resultsNames(dds)
     
  levels(dds$induced)
  levels(dds$condition)
     

  res_1S <- results(dds, name = "induced_yes_vs_no", 
                    lfcThreshold = thr_lfc, alpha = thr_alpha)
  summary(res_1S)
     
  res_2S <- results(dds, 
                    contrast = list(c("induced_yes_vs_no", resultsNames(dds)[4])), 
                    lfcThreshold = thr_lfc, alpha = thr_alpha)
  summary(res_2S)
     
  res_condition1Svs2S <- results(dds, name = resultsNames(dds)[4], 
                                    lfcThreshold = thr_lfc, alpha = thr_alpha)
  summary(res_condition1Svs2S)
     
  res_condition1Svs2S_gg <- data.frame(sgRNA = rownames(dds), 
                                       target = targets$target_operon_names, 
                                       log2FC_1S = res_1S$log2FoldChange, 
                                       log2FC_2S = res_2S$log2FoldChange, 
                                       padj = res_condition1Svs2S$padj)
  
  pX_title <- paste(coldata1$strain[1], " 1S:",
                       levels(dds$condition)[1], " vs. 2S:",
                       levels(dds$condition)[2], sep="")
     
  pX <- ggplot(res_condition1Svs2S_gg, 
               aes(log2FC_1S, log2FC_2S, 
               label = ifelse(padj < thr_alpha, target, NA),
               size = -log10(padj), 
               fill = padj < thr_alpha)) + 
       geom_vline(xintercept = 0) + 
       geom_hline(yintercept = 0) + 
       geom_abline(slope = 1, intercept = 0) +
       geom_point(shape = 21) + 
       geom_text_repel(fontface = "italic", size = 3) +
       scale_fill_manual(values = c("grey", "darkred")) +
       labs(title = pX_title) +
       theme_classic(base_size = 16)
     
     
  pdf(paste(path.out,"/", "Exp", i, "_",coldata1$strain[1],"_",
            coldata1$media[1], "_vs_",coldata1$media[16], ".pdf", sep=""), 
            width=15, height=7)
  plot(pX)
  dev.off()
}
   
###############################
############<xxxxxxxxxxxxxxxxxxxxxxxx#######
###########################################   
   
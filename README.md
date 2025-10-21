#Load Library
```{r}
library(EnhancedVolcano)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(collapse)
library(edgeR)
library(DEFormats)
library(GSA)
library(PPInfer)
library(plotly)
library(readxl)
library(biomaRt)
```


#Set directory and load
```{r}
YH1 <- read_excel("Directory")
```

#Data processing
```{r}
colnames(YH1)
target = c(2,3,4:7)
selected = YH1[,target]
YH1_collapse = collap(selected, ~`Gene_Symbol`, FUN=max)
#isUnique(YH1_collapse[1])
genename = YH1_collapse[1]
YH1_MT = as.matrix(YH1_collapse[2:ncol(YH1_collapse)])
rownames(YH1_MT) <- genename$`Gene_Symbol`
YH1_final = YH1_MT

colnames(YH1_final)

colsna = c("YHIM1053_1", "YHIM1053_2", "YHIM1035_1", "YHIM1035_2")
YH1_final_col = YH1_final
colnames(YH1_final_col) = colsna
```

```{r}
YH1_mart4 <- YH1_final_col
```

# 1053 vs 1035
```{r}
head(YH1_mart4)
1053 = YH1_mart4[,1:2]
1035 = YH1_mart4[,3:4]



counts =as.matrix(cbind(1053, 1035))
counts = round(counts)
dim(counts)
group = c(rep(c("A"),ncol(1053)),rep(c("B"),ncol(1035)))
dge = DGEList(counts, group = group)
dge
```

```{r}
dds = as.DESeqDataSet(dge)
identical(dge, as.DGEList(dds)) # if it is FALSE, please process as below
dds <- dds[rowSums(counts(dds)) > 1,]



dds1 = DESeqDataSetFromMatrix(counts, data.frame(condition=group), ~ condition)
dds2 = DESeqDataSetFromMatrix(counts, data.frame(condition=group), ~ condition)
identical(dds1, dds2) #should be TRUE

#DEseq2
dds <- DESeq(dds, betaPrior=FALSE)
res1 <- results(dds,
                contrast = c('group','A','B'))
res1 <- lfcShrink(dds,
                  contrast = c('group','A','B'), 
                  res=res1,
                  type = "normal")
```


#dependency
```{r}
options(bitmapType='cairo')
```

#Fig. 5B
```{r, fig.height=18, fig.width=16}
res1 <- na.omit(res1)  # 결측값이 있는 행을 제거

res1$color <- ifelse(res1$log2FoldChange > 0 & res1$pvalue < 0.05, "red", "black")
# NA 값이 있을 경우 기본 색상("black")으로 채우기
res1$color[is.na(res1$color)] <- "black"


pdf("Refractory_SDPD_23samples.pdf", width=12, height = 8)
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'pd vs pr',
                pCutoff = 0.05,
                ylim = c(0,15),
                FCcutoff = 0.25,
                pointSize = 3.0,
                drawConnectors = TRUE,
                colConnectors = 'grey30',
                boxedLabels = T,
               labSize = 6.0,
selectLab = c("MET","AXL","HGF","KRAS","EGFR"))

dev.off()
```


```{r, fig.height=18, fig.width=16}
PI3K_AKT_MTOR_SIGNALING <- c(
  "ACACA","ACTR2","ACTR3","ADCY2","GRK2","AKT1","AKT1S1","AP2M1","ARF1","ARHGDIA",
  "ARPC3","ATF1","CAB39","CAB39L","CALR","CAMK4","CDK1","CDK2","CDK4","CDKN1A",
  "CDKN1B","CFL1","CLTC","CSNK2B","CXCR4","DAPP1","DDIT3","DUSP3","E2F1","ECSIT",
  "EGFR","EIF4E","FASLG","FGF17","FGF22","FGF6","GNA14","GNGT1","GRB2","GSK3B",
  "HRAS","HSP90B1","IL2RG","IL4","IRAK4","ITPR2","LCK","MAP2K3","MAP2K6","MAP3K7",
  "MAPK1","MAPK10","MAPK8","MAPK9","MAPKAP1","MKNK1","MKNK2","MYD88","NCK1",
  "NFKBIB","NGF","NOD1","PAK4","PDK1","PFN1","PIK3R3","PIKFYVE","PIN1","PITX2",
  "PLA2G12A","PLCB1","PLCG1","PPP1CA","PPP2R1B","PRKAA2","PRKAG1","PRKAR2A",
  "PRKCB","PTEN","PTPN11","RAC1","RAF1","RALB","RIPK1","RIT1","RPS6KA1","RPS6KA3",
  "RPTOR","SFN","SLA","SLC2A1","SMAD2","SQSTM1","STAT2","TBK1","THEM4","TIAM1",
  "TNFRSF1A","TRAF2","TRIB3","TSC2","UBE2D3","UBE2N","VAV3","YWHAB"
)




MAPK_SIGNALING_PATHWAY<- c(
  "AKT3","PLA2G4B","RASGRP1","RASGRP2","CACNG3","CACNG2","TAB1","MAP3K2","GADD45G","DUSP14","MAP4K1","DUSP10","CHP1","CHUK","RASGRP4","PLA2G4E","MAP3K8","ATF2","CRK","CRKL","MAPK14","DAXX","GADD45A","DDIT3","DUSP1","DUSP2","DUSP3","DUSP4","DUSP5","DUSP6","DUSP7","DUSP8","DUSP9","EGF","EGFR","ELK1","ELK4","AKT1","AKT2","MECOM","FGF1","FGF2","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGF10","FGF11","FGF12","FGF13","FGF14","FGFR1","FGFR3","FGFR2","FGFR4","RRAS2","MRAS","TAB2","FLNA","MAPK8IP3","FLNB","FLNC","FOS","MAPK8IP2","RASGRP3","PLA2G2D","FGF20","FGF21","FGF22","CACNG5","CACNG4","RPS6KA6","GNA12","MKNK2","GRB2","PLA2G2E","NR4A1","HRAS","HSPA1A","HSPA1B","HSPA1L","HSPA2","HSPA6","HSPA8","HSPB1","FAS","IKBKB","IL1A","IL1B","IL1R1","FASLG","JUN","JUND","KRAS","PLA2G2C","STMN1","ARRB1","ARRB2","MAPT","MAX","MEF2C","MAP3K1","MAP3K3","MAP3K4","MAP3K5","MAP3K11","MOS","MYC","GADD45B","ATF4","NF1","NFATC2","NFATC4","NFKB1","NFKB2","NGF","NRAS","NTF3","NTF4","NTRK1","NTRK2","PLA2G3","PAK1","PAK2","ECSIT","TAOK3","PDGFA","PDGFB","PDGFRA","PDGFRB","NLK","MAP3K20","PLA2G1B","PLA2G2A","PLA2G4A","PLA2G5","PPM1A","PPM1B","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP3R2","PPP5C","PRKACA","PRKACB","PRKACG","PRKCA","PRKCB","CACNA2D3","PRKCG","MAPK1","MAPK3","GNG12","MAPK7","MAPK8","MAPK11","MAPK9","MAPK10","MAPK13","MAP2K1","MAP2K2","MAP2K3","MAP2K5","MAP2K6","MAP2K7","PRKX","TAOK1","PTPN7","PTPRR","MAP4K2","RAC1","RAC2","RAC3","RAF1","RAP1A","RAP1B","RASA1","RASA2","RASGRF1","RASGRF2","CACNG8","CACNG7","CACNG6","RELA","RELB","RPS6KA1","RPS6KA2","RPS6KA3","RRAS","BDNF","MAPK12","CHP2","MAP2K4","PLA2G2F","SOS1","SOS2","SRF","BRAF","STK3","STK4","MAP3K7","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TNF","TNFRSF1A","TP53","TRAF2","TRAF6","CACNA1A","CACNA1B","CACNA1C","CACNA1D","CACNA1E","CACNA1F","MAP3K12","CACNA1S","CACNA2D1","CACNB1","CACNB2","CACNB3","CACNB4","IL1R2","CACNG1","MAPKAPK3","FGF23","DUSP16","PLA2G12A","CASP3","PLA2G6","PLA2G10","PLA2G12B","PTPN5","MAP4K3","IKBKG","MAPKAPK5","MKNK1","LAMTOR3","JMJD7-PLA2G4B","FGF18","FGF17","FGF16","CACNA1I","CACNA1H","CACNA1G","RPS6KA4","MAP3K14","MAP3K6","MAP3K13","RPS6KA5","CACNA2D2","MAPKAPK2","CD14","TAOK2","CACNA2D4","MAP4K4","MAPK8IP1","RAPGEF2","CDC25B","FGF19","CDC42"
)


CELL_CYCLE <- c(
  "CDK2","CDK4","CDK6","CDK7","CDKN1A","CDKN1B","STAG1","CDKN1C","CDKN2A","CDKN2B",
  "CDKN2C","CDKN2D","ANAPC10","MAD2L2","STAG2","PTTG2","GADD45G","DBF4","YWHAQ",
  "CHEK1","CHEK2","CREBBP","GADD45A","E2F1","E2F2","E2F3","E2F4","E2F5","EP300",
  "ORC6","ORC3","CDC26","ABL1","ANAPC13","SMC1B","SFN","GSK3B","ANAPC2","ANAPC4",
  "HDAC1","HDAC2","MAD2L1","SMAD2","SMAD3","SMAD4","MCM2","MCM3","MCM4","MCM5",
  "MCM6","MCM7","MDM2","MYC","GADD45B","ATM","WEE2","ORC1","ORC2","ORC4","ORC5",
  "PCNA","FZR1","ANAPC5","ANAPC7","ANAPC11","PLK1","ATR","PRKDC","RAD21","RB1",
  "RBL1","RBL2","CCND1","ANAPC1","SKP1","SKP2","BUB1","BUB1B","TFDP1","TFDP2",
  "TGFB1","TGFB2","TGFB3","TP53","TTK","SKP1P2","WEE1","YWHAB","YWHAE","YWHAG",
  "YWHAH","YWHAZ","ZBTB17","SMC1A","CDC7","CDC45","MAD1L1","CUL1","CCNB3",
  "CDC14B","CDC14A","CDC23","CDC16","CCNA2","CCNA1","CCNB1","CCND2","CCND3",
  "CCNE1","CCNH","PKMYT1","SMC3","CCNB2","CCNE2","BUB3","PTTG1","ESPL1","CDK1",
  "CDC6","CDC20","CDC25A","CDC25B","CDC25C","CDC27","RBX1")

PATHWAYS_IN_CANCER<- c(
  "AKT3","CDK2","CDK4","CDK6","CDKN1A","CDKN1B","CDKN2A","APC2","CDKN2B","LAMC3",
  "TFG","PIAS3","CEBPA","RALBP1","RASSF1","FZD10","EGLN2","EGLN3","CHUK","CKS1B",
  "COL4A1","COL4A2","COL4A4","COL4A6","CREBBP","CRK","CRKL","CSF1R","CSF2RA",
  "CSF3R","CTBP1","CTBP2","CTNNA1","CTNNA2","CTNNB1","DAPK1","DAPK3","DCC",
  "DVL1","DVL2","DVL3","E2F1","E2F2","E2F3","EGF","EGFR","EP300","EPAS1","ERBB2",
  "AKT1","AKT2","ETS1","MECOM","FGF1","FGF2","FGF3","FGF4","FGF5","FGF6","FGF7",
  "FGF8","FGF9","FGF10","FGF11","FGF12","FGF13","FGF14","FGFR1","FGFR3","FGFR2",
  "FH","VEGFD","LAMB4","FOXO1","FLT3","FLT3LG","FN1","FOS","PIK3R5","DAPK2",
  "CBLC","MTOR","ABL1","FZD2","APPL1","FGF20","FGF21","FGF22","STK36","GLI1",
  "GLI2","GLI3","LAMA1","GRB2","CTNNA3","GSK3B","GSTP1","MSH6","HDAC1","HDAC2",
  "HGF","HIF1A","APC","HRAS","BIRC2","BIRC3","XIAP","BIRC5","HSP90AA1","HSP90AB1",
  "IGF1","IGF1R","KLK3","FAS","IKBKB","FASLG","IL6","CXCL8","ITGA6","AR","ITGA2",
  "ITGA2B","ITGA3","ITGAV","ITGB1","ARAF","JAK1","JUN","JUP","KIT","KRAS","RHOA",
  "LAMA2","LAMA3","LAMA4","LAMA5","LAMB1","LAMB2","LAMB3","LAMC1","LAMC2","ARNT",
  "SMAD2","SMAD3","SMAD4","MAX","MDM2","MET","KITLG","MITF","MLH1","MMP1","MMP2",
  "MMP9","MSH2","MSH3","MYC","NFKB1","NFKB2","NFKBIA","NKX3-1","NOS2","NRAS",
  "NTRK1","LEF1","WNT16","PDGFA","PDGFB","PDGFRA","PIAS4","PDGFRB","SUFU","PGF",
  "PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PLCG1","PLCG2","PLD1",
  "PML","CYCS","WNT4","EGLN1","PPARD","PPARG","PRKCA","PRKCB","PRKCG","MAPK1",
  "MAPK3","MAPK8","MAPK9","MAPK10","MAP2K1","MAP2K2","BAD","PTCH1","PTEN","PTGS2",
  "PTK2","BAX","RAC1","RAC2","RAC3","RAD51","RAF1","RALA","RALB","RALGDS","RARA",
  "RARB","RB1","CCND1","BCL2","RELA","RET","BCL2L1","BCR","RXRA","RXRB","RXRG",
  "BID","HHIP","SHH","BMP2","SKP2","SLC2A1","BMP4","SMO","SOS1","SOS2","SPI1",
  "BRAF","BRCA2","STAT1","STAT3","STAT5A","STAT5B","STK4","ELOC","ELOB","TCF7",
  "TCF7L2","TGFA","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TP53","TPM3","TPR",
  "HSP90B1","TRAF1","TRAF2","TRAF3","TRAF5","TRAF6","VEGFA","VEGFB","VEGFC","VHL",
  "WNT1","WNT2","WNT3","WNT5A","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT10B",
  "WNT11","WNT2B","WNT9A","WNT9B","ZBTB16","PAX8","FZD5","FZD3","CCDC6","NCOA4",
  "WNT10A","FGF23","WNT5B","AXIN1","AXIN2","FZD1","FZD4","FZD6","FZD7","FZD8",
  "FZD9","TCF7L1","RASSF5","CASP3","CASP8","CASP9","CUL2","PIK3R3","IKBKG","PIAS1",
  "RUNX1","RUNX1T1","PTCH2","CBL","CBLB","FADD","FGF18","FGF17","FGF16","CCNA1",
  "WNT3A","CCNE1","PIAS2","CCNE2","TRAF4","ARNT2","FGF19","RBX1","CDC42","CDH1"
)


ADHERENS_JUNCTION <- c(
  "WASF2","BAIAP2","SORBS1","WASF3","SSX2IP","CREBBP","CSNK2A1","CSNK2A2","CSNK2B",
  "CTNNA1","CTNNA2","CTNNB1","CTNND1","EGFR","EP300","ERBB2","FER","FGFR1","FYN",
  "NECTIN3","CTNNA3","IGF1R","INSR","RHOA","LMO7","SMAD2","SMAD3","SMAD4","MET",
  "AFDN","LEF1","NLK","ACP1","MAPK1","MAPK3","PARD3","PTPN1","PTPN6","PTPRB",
  "PTPRF","PTPRJ","PTPRM","NECTIN1","NECTIN2","RAC1","RAC2","RAC3","ACTB","SNAI2",
  "SNAI1","SRC","MAP3K7","TCF7","TCF7L2","TGFBR1","TGFBR2","TJP1","ACTG1","VCL",
  "WAS","YES1","ACTN4","NECTIN4","TCF7L1","ACTN1","ACTN2","IQGAP1","ACTN3","WASF1",
  "WASL","FARP2","CDC42","CDH1"
)


NON_SMALL_CELL_LUNG_CANCER <- c(
  "AKT3","CDK4","CDK6","CDKN2A","RASSF1","E2F1","E2F2","E2F3","EGF","EGFR",
  "ERBB2","AKT1","AKT2","FHIT","FOXO3","PIK3R5","GRB2","HRAS","ARAF","KRAS",
  "NRAS","PDPK1","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PLCG1",
  "PLCG2","PRKCA","PRKCB","PRKCG","MAPK1","MAPK3","MAP2K1","MAP2K2","BAD",
  "RAF1","RARB","RB1","CCND1","RXRA","RXRB","RXRG","SOS1","SOS2","BRAF","STK4",
  "TGFA","TP53","RASSF5","CASP9","PIK3R3"
)

MET_PROMOTES_CELL_MOTILITY <- c(
  "ITGA3","HGF","LAMC3","LAMA3","LAMC2","COL11A1","COL5A3","LAMB1","CRKL",
  "LAMA1","MET","RAPGEF1","COL1A1","GAB1","LAMA4","FN1","RAP1A","DOCK7","RAP1B",
  "COL5A1","LAMA5","TNS4","LAMC1","TNS3","RAC1","COL2A1","ITGB1","ITGA2","COL1A2",
  "CRK","COL3A1","PTK2","COL24A1","LAMB2","GRB2","LAMA2","COL27A1","LAMB3","SRC",
  "COL11A2","COL5A2"
)

EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE <- c(
  "AKT3","BCL2L11","EGF","EGFR","EIF4E","EIF4EBP1","ERBB2","ERBB3","AKT1","AKT2",
  "FGF2","FGFR3","FGFR2","RRAS2","MRAS","FOXO3","MTOR","GAB1","SHC2","GAS6","GRB2",
  "GSK3B","HGF","NRG1","HRAS","IGF1","IGF1R","IL6","IL6R","ARAF","JAK1","JAK2","KDR",
  "KRAS","SHC4","MET","MYC","NF1","NRAS","PDGFA","PDGFB","PDGFRA","PDGFRB","PDPK1",
  "PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","PLCG1","SHC3","PLCG2","PRKCA","PRKCB",
  "AXL","PRKCG","MAPK1","MAPK3","PDGFC","MAP2K1","MAP2K2","BAD","PTEN","BAX","RAF1",
  "CCND1","BCL2","BCL2L1","RPS6","RPS6KB1","RPS6KB2","RRAS","SHC1","SOS1","SOS2",
  "SRC","BRAF","STAT3","TGFA","VEGFA","PDGFD","PIK3R3","EIF4E2","NRG2"
)

SIGNAL_TRANSDUCTION_BY_L1<- c(
  "ITGA2B","CSNK2A2","FGFR1","NRP1","MAPK1","CSNK2A1","MAPK3","MAP2K2","RAC1",
  "ITGAV","ITGA9","EGFR","PAK1","NCAM1","ITGB1","VAV2","ITGA5","MAP2K1","L1CAM",
  "CSNK2B","ITGB3"
)


gene_sets3 <- list(
EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE=EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE,
MAPK_SIGNALING_PATHWAY=MAPK_SIGNALING_PATHWAY,
PI3K_AKT_MTOR_SIGNALING =PI3K_AKT_MTOR_SIGNALING,
NON_SMALL_CELL_LUNG_CANCER=NON_SMALL_CELL_LUNG_CANCER,
CELL_CYCLE=CELL_CYCLE,
PATHWAYS_IN_CANCER=PATHWAYS_IN_CANCER,
ADHERENS_JUNCTION=ADHERENS_JUNCTION,
SIGNAL_TRANSDUCTION_BY_L1=SIGNAL_TRANSDUCTION_BY_L1)


```

#DEG patwhay_geneset_1053
```{r}
Post_translational_protein_phosphorylation <- c("AMBN", "APOE", "C4A", "C4B", "CHGB", "CP", "ENAM", "F5", "FGA", "FGG", "IGFBP5", "ITIH2", "MEPE", "MGAT4A", "NOTUM", "PROC", "SCG3")
Interaction_of_cadherin  <- c("CDH10", "CDH12", "CDH17", "CDH6", "CDH7", "CDH9")
NRG2_activated_ERBB4_binds_GABA_receptors_through_GABRA1 <- c("ERBB4", "GABRA1", "GABRB2", "GABRG2")

epithelial_cell_differentiation <- c(
  "AQP3", "CAMSAP3", "CASP14", "CDH3", "CDSN", "CERS3", "CES1", "COL18A1", "CSTA", "DSP",
  "EVPL", "F2RL1", "FGFR2", "FLG", "FOXN1", "FZD7", "GATA6", "HEY2", "IRF6", "KDF1",
  "KRT13", "KRT15", "KRT17", "KRT19", "KRT6C", "KRT7", "KRT78", "KRT79", "LHX1",
  "MACROH2A2", "NELL1", "NRG1", "OVOL1", "PAX6", "PLAAT1", "PRKX", "RAB25", "RHCG",
  "S1PR3", "SALL1", "SMO", "ST14", "STC1", "SULT1B1", "SULT2B1", "TBX1", "VDR",
  "WNT5B", "WNT7A", "WNT7B", "ZBED2"
)

DNA_binding_TF_activity <- c(
  "ADCY1", "ALK", "ALX1", "ALX3", "ARNT2", "ARNTL2", "BARX2", "BEX2", "BMP7", "BNC1",
  "CAMK2A", "CARD11", "CX3CL1", "DLX6", "DMRT2", "DMRTA2", "E2F5", "ENPP1", "FOSB", "FOSL1",
  "FOXD1", "FOXE1", "FOXF2", "FOXL1", "FOXL2", "FOXN1", "FOXQ1", "GATA6", "GDNF", "HES2",
  "HEY2", "HOXC13", "HOXD13", "IRF6", "IRX2", "IRX4", "JMY", "LEF1", "LHX1", "LHX5",
  "LMX1B", "MDFI", "MSC", "MSX1", "MSX2", "NHLH2", "NKX1-2", "NKX2-8", "NKX3-1", "NKX3-2",
  "NLRP2", "NTSR1", "OVOL1", "PAX6", "PLAG1", "POU3F1", "POU4F1", "PRDM16", "PRKCZ", "PRNP",
  "PRRX2", "RORB", "SALL1", "SLCO3A1", "SMO", "SNAI2", "SOX15", "SOX7", "SP9", "TBX1",
  "TFAP2C", "TRIML2", "VAX1", "VAX2", "VDR", "ZBED2", "ZFP28", "ZFP30", "ZFY", "ZIC2",
  "ZIC5", "ZNF100", "ZNF135", "ZNF215", "ZNF253", "ZNF257", "ZNF347", "ZNF415", "ZNF429", "ZNF43",
  "ZNF443", "ZNF486", "ZNF516", "ZNF528", "ZNF568", "ZNF573", "ZNF578", "ZNF607", "ZNF626", "ZNF66",
  "ZNF665", "ZNF667", "ZNF675", "ZNF681", "ZNF682", "ZNF701", "ZNF724", "ZNF726", "ZNF730", "ZNF736",
  "ZNF737", "ZNF750", "ZNF793", "ZNF799", "ZNF808", "ZNF813", "ZNF818P", "ZNF829", "ZNF844", "ZNF85",
  "ZNF90", "ZNF91", "ZNF93", "ZSCAN18"
)

CEll_Cell_adhesion <- c(
  "ADAM19", "ADD2", "ADTRP", "ALOX12", "ASS1", "BMP7", "CAMSAP3", "CARD11", "CCL21", "CD177",
  "CD99", "CDH1", "CDH3", "CDSN", "CERCAM", "CLDN2", "CLDN8", "COL13A1", "CSTA", "CX3CL1",
  "DMTN", "DSC3", "DSP", "EGFR", "FAT2", "FUT3", "FXYD5", "IGF2", "IGSF11", "IGSF9",
  "LEF1", "LGALS9B", "LRRC4", "MAGI1", "MDGA1", "MMRN1", "MPZL2", "NECTIN1", "NECTIN4", "NEO1",
  "NLGN1", "NLGN4Y", "NT5E", "NTNG2", "PERP", "PKP1", "PKP2", "PKP3", "PRKCZ", "PRNP",
  "PTPRG", "ROBO3", "SDK1", "SPINT2", "STXBP6", "TENM1", "TENM2", "TINAGL1", "TRPV4", "UBASH3B",  "WNT7B", "ZDHHC2")

epidermal_cell_differentiation <- c(
  "AQP3", "CASP14", "CDH3", "CDSN", "CERS3", "CSTA", "DSP", "EVPL", "FLG", "GATA6",
  "KDF1", "KRT6C", "KRT7", "KRT78", "KRT79", "MACROH2A2", "OVOL1", "ST14", "SULT2B1",
  "VDR", "ZBED2"
)

Extracellular_matrix_organization<- c(
  "ADAM19", "ADAMTS16", "BMP7", "CAPN11", "CAPNS2", "CDH1", "COL13A1", "COL17A1", "COL18A1","COL6A1", "COL6A2", "COL6A3", "EMILIN3", "ITGB4", "LAMA1", "LAMB3", "LAMC3", "MATN3", "MMP1", "MMP15", "MMP17", "P3H2", "PDGFA", "PXDN", "SCUBE3")


Degradation_of_the_extracellular_matrix <- c(
  "ADAMTS16", "CAPN11", "CAPNS2", "CDH1", "COL13A1", "COL17A1", "COL18A1",
  "COL6A1", "COL6A2", "COL6A3", "LAMB3", "MMP1", "MMP15", "MMP17", "SCUBE3")

Hippo_signaling_pathway <- c(
  "AMOT", "BMP7", "CCND1", "CDH1", "FZD7", "GDF7", "LEF1", "PPP2R2C",
  "PRKCZ", "SNAI2", "WNT3A", "WNT5B", "WNT7A", "WNT7B", "WWC1"
)

regulation_deveolpmental <- c(
  "ABCG1", "ADGRV1", "AGT", "ALOX5", "APOE", "ASB4", "ASCL1", "ASIC2", "ATOH8", "ATP10A",
  "BMP5", "BMP6", "BRINP2", "BRINP3", "BRSK1", "CALCA", "CD36", "CDKN2A", "CEACAM1", "CFTR",
  "CHI3L1", "CHN1", "CHRNA3", "CHRNB2", "CNMD", "CR1", "CYP1B1", "DCC", "ENAM", "EPHA3",
  "ERBB4", "FBN2", "FGA", "FGB", "FGF13", "FGG", "GPC6", "GRID2", "GRIN3A", "HES6",
  "HES7", "HGF", "IHH", "INSM1", "ISL1", "ISLR2", "KIF1A", "KLHL41", "KLKB1", "LEP",
  "LFNG", "LGALS12", "LILRB1", "LILRB2", "LILRB3", "LMO3", "LZTS1", "MAPT", "METRN", "MSTN",
  "MTURN", "MYLK3", "NCMAP", "NLGN3", "NOTUM", "NRXN1", "OLFM4", "PAK3", "PCK1", "PCP4",
  "PDE3B", "PKHD1", "PPFIA2", "PROC", "PTN", "PTPRD", "RANBP3L", "RARRES2", "RBP4", "RFLNA",
  "RHOH", "RND2", "RPS6KA6", "SARM1", "SCIN", "SEMA6B", "SEZ6", "SFRP4", "SLITRK2", "SLITRK4",
  "SPOCK2", "SRCIN1", "ST7", "SYT1", "SYT17", "SYT2", "SYT4", "TENM4", "TENT5C", "TGM2",
  "TMEM100", "TMEM176A", "TMEM176B", "TMEM178A", "TRIM46", "UNC13A", "ZC4H2", "ZDHHC15"
)
```


#Fig.5C
```{r}
library(GSVA)
```

```{r}
gene_sets5 <- list(
epithelial_cell_differentiation = epithelial_cell_differentiation,
DNA_binding_TF_activity = DNA_binding_TF_activity,
CEll_Cell_adhesion = CEll_Cell_adhesion,
epidermal_cell_differentiation=epidermal_cell_differentiation,
Extracellular_matrix_organization=Extracellular_matrix_organization,
Degradation_of_the_extracellular_matrix = Degradation_of_the_extracellular_matrix,
Hippo_signaling_pathway =Hippo_signaling_pathway)

```

```{r, fig.width= 8, fig.height= 5}
gsva_results <- gsva(YH1_final, gene_sets5, method = "gsva")

# 4. Heatmap 생성
# 스코어링 결과를 정규화 및 시각화
Heatmap(gsva_results,
        name = "Score", 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Gene Sets",
        column_title = "Samples",
        heatmap_legend_param = list(title = "GSVA Score"))

```

#Fig.5D

```{r}
library(GSVA)
```

```{r}
gene_sets6 <- list(
Post_translational_protein_phosphorylation=Post_translational_protein_phosphorylation,
regulation_deveolpmental = regulation_deveolpmental,
epithelial_cell_differentiation = epithelial_cell_differentiation,
CEll_Cell_adhesion = CEll_Cell_adhesion,
Extracellular_matrix_organization=Extracellular_matrix_organization,
Degradation_of_the_extracellular_matrix = Degradation_of_the_extracellular_matrix)
```

```{r, fig.width= 8, fig.height= 5}
gsva_results <- gsva(YH1_final, gene_sets6, method = "gsva")

# 4. Heatmap 생성
# 스코어링 결과를 정규화 및 시각화
Heatmap(gsva_results,
        name = "Score", 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Gene Sets",
        column_title = "Samples",
        heatmap_legend_param = list(title = "GSVA Score"))

```

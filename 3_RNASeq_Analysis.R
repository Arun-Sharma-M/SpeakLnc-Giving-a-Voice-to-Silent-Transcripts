#Input_Entire gene list
##              ##              ##              ##              ##              ##              ##              ##              ##
countData<-read.table("SRA_CHOL_lncRNA_readcounts_forDeseq.tsv",sep="\t",header=T,row.names=1)
##              ##              ##              ##              ##              ##              ##              ##              ##
dim(countData)
#[1] 127802     34
countData[1:5,1:3]

#Input_gene length # Make sure of lncRNA and gene lengths file
##              ##              ##              ##              ##              ##              ##              ##              ##
length<-read.table("lncRNA_start_end_FULL.csv",sep=",",header=T,row.names=1) 
##              ##              ##              ##              ##              ##              ##              ##              ##
dim(length)
# [1] 127802      2

RC2TPM <- data.frame(countData) # Creating data
######################################################################################################################################
############################################# [READCOUNTS TO TPM CONVERSION] #########################################################
######### Introducing gene_length_in_Kb Column by performing (END - START values by 1000) to get the gene length in Kb ###############
######################################################################################################################################
######################################################################################################################################
RC2TPM[,35] = (length[,2]-length[,1])/1000
colnames(RC2TPM)[35] <- "gene_length_in_Kb"
RC2TPM[RC2TPM == 0] <- 0.1
#cHANGE: 
# 1. RC2TPM$C TO RC2TPM$CC
# 2. RC2TPM$T TO RC2TPM$TC
########################### Introducing RPKT Columns ###########################
RC2TPM$RPKC1	<-	RC2TPM$CC1	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC1
RC2TPM$RPKC2	<-	RC2TPM$CC2	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC2
RC2TPM$RPKC3	<-	RC2TPM$CC3	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC3
RC2TPM$RPKC4	<-	RC2TPM$CC4	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC4
RC2TPM$RPKC5	<-	RC2TPM$CC5	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC5
RC2TPM$RPKC6	<-	RC2TPM$CC6	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC6
RC2TPM$RPKC7	<-	RC2TPM$CC7	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC7
RC2TPM$RPKC8	<-	RC2TPM$CC8	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC8
RC2TPM$RPKC9	<-	RC2TPM$CC9	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC9
RC2TPM$RPKC10	<-	RC2TPM$CC10	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC10
RC2TPM$RPKC11	<-	RC2TPM$CC11	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC11
RC2TPM$RPKC12	<-	RC2TPM$CC12	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC12
RC2TPM$RPKC13	<-	RC2TPM$CC13	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC13
RC2TPM$RPKC14	<-	RC2TPM$CC14	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC14
RC2TPM$RPKC15	<-	RC2TPM$CC15	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC15
RC2TPM$RPKC16	<-	RC2TPM$CC16	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC16
RC2TPM$RPKC17	<-	RC2TPM$CC17	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC17
RC2TPM$RPKC18	<-	RC2TPM$CC18	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC18
RC2TPM$RPKT1	<-	RC2TPM$TC1	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT1
RC2TPM$RPKT2	<-	RC2TPM$TC2	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT2
RC2TPM$RPKT3	<-	RC2TPM$TC3	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT3
RC2TPM$RPKT4	<-	RC2TPM$TC4	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT4
RC2TPM$RPKT5	<-	RC2TPM$TC5	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT5
RC2TPM$RPKT6	<-	RC2TPM$TC6	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT6
RC2TPM$RPKT7	<-	RC2TPM$TC7	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT7
RC2TPM$RPKT8	<-	RC2TPM$TC8	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT8
RC2TPM$RPKT9	<-	RC2TPM$TC9	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT9
RC2TPM$RPKT10	<-	RC2TPM$TC10	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT10
RC2TPM$RPKT11	<-	RC2TPM$TC11	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT11
RC2TPM$RPKT12	<-	RC2TPM$TC12	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT12
RC2TPM$RPKT13	<-	RC2TPM$TC13	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT13
RC2TPM$RPKT14	<-	RC2TPM$TC14	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT14
RC2TPM$RPKT15	<-	RC2TPM$TC15	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT15
RC2TPM$RPKT16	<-	RC2TPM$TC16	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT16

######################### Introducing PMRFT Columns ############################
RC2TPMC.1	<-	(sum(RC2TPM[,	'RPKC1'],	na.rm	=	TRUE))/10^6
RC2TPMC.2	<-	(sum(RC2TPM[,	'RPKC2'],	na.rm	=	TRUE))/10^6
RC2TPMC.3	<-	(sum(RC2TPM[,	'RPKC3'],	na.rm	=	TRUE))/10^6
RC2TPMC.4	<-	(sum(RC2TPM[,	'RPKC4'],	na.rm	=	TRUE))/10^6
RC2TPMC.5	<-	(sum(RC2TPM[,	'RPKC5'],	na.rm	=	TRUE))/10^6
RC2TPMC.6	<-	(sum(RC2TPM[,	'RPKC6'],	na.rm	=	TRUE))/10^6
RC2TPMC.7	<-	(sum(RC2TPM[,	'RPKC7'],	na.rm	=	TRUE))/10^6
RC2TPMC.8	<-	(sum(RC2TPM[,	'RPKC8'],	na.rm	=	TRUE))/10^6
RC2TPMC.9	<-	(sum(RC2TPM[,	'RPKC9'],	na.rm	=	TRUE))/10^6
RC2TPMC.10	<-	(sum(RC2TPM[,	'RPKC10'],	na.rm	=	TRUE))/10^6
RC2TPMC.11	<-	(sum(RC2TPM[,	'RPKC11'],	na.rm	=	TRUE))/10^6
RC2TPMC.12	<-	(sum(RC2TPM[,	'RPKC12'],	na.rm	=	TRUE))/10^6
RC2TPMC.13	<-	(sum(RC2TPM[,	'RPKC13'],	na.rm	=	TRUE))/10^6
RC2TPMC.14	<-	(sum(RC2TPM[,	'RPKC14'],	na.rm	=	TRUE))/10^6
RC2TPMC.15	<-	(sum(RC2TPM[,	'RPKC15'],	na.rm	=	TRUE))/10^6
RC2TPMC.16	<-	(sum(RC2TPM[,	'RPKC16'],	na.rm	=	TRUE))/10^6
RC2TPMC.17	<-	(sum(RC2TPM[,	'RPKC17'],	na.rm	=	TRUE))/10^6
RC2TPMC.18	<-	(sum(RC2TPM[,	'RPKC18'],	na.rm	=	TRUE))/10^6
RC2TPM.1	<-	(sum(RC2TPM[,	'RPKT1'],	na.rm	=	TRUE))/10^6
RC2TPM.2	<-	(sum(RC2TPM[,	'RPKT2'],	na.rm	=	TRUE))/10^6
RC2TPM.3	<-	(sum(RC2TPM[,	'RPKT3'],	na.rm	=	TRUE))/10^6
RC2TPM.4	<-	(sum(RC2TPM[,	'RPKT4'],	na.rm	=	TRUE))/10^6
RC2TPM.5	<-	(sum(RC2TPM[,	'RPKT5'],	na.rm	=	TRUE))/10^6
RC2TPM.6	<-	(sum(RC2TPM[,	'RPKT6'],	na.rm	=	TRUE))/10^6
RC2TPM.7	<-	(sum(RC2TPM[,	'RPKT7'],	na.rm	=	TRUE))/10^6
RC2TPM.8	<-	(sum(RC2TPM[,	'RPKT8'],	na.rm	=	TRUE))/10^6
RC2TPM.9	<-	(sum(RC2TPM[,	'RPKT9'],	na.rm	=	TRUE))/10^6
RC2TPM.10	<-	(sum(RC2TPM[,	'RPKT10'],	na.rm	=	TRUE))/10^6
RC2TPM.11	<-	(sum(RC2TPM[,	'RPKT11'],	na.rm	=	TRUE))/10^6
RC2TPM.12	<-	(sum(RC2TPM[,	'RPKT12'],	na.rm	=	TRUE))/10^6
RC2TPM.13	<-	(sum(RC2TPM[,	'RPKT13'],	na.rm	=	TRUE))/10^6
RC2TPM.14	<-	(sum(RC2TPM[,	'RPKT14'],	na.rm	=	TRUE))/10^6
RC2TPM.15	<-	(sum(RC2TPM[,	'RPKT15'],	na.rm	=	TRUE))/10^6
RC2TPM.16	<-	(sum(RC2TPM[,	'RPKT16'],	na.rm	=	TRUE))/10^6

####################### Introducing TPMT Columns ###############################
RC2TPM$TPMC1	<-	RC2TPM$RPKC1	/	RC2TPMC.1
RC2TPM$TPMC2	<-	RC2TPM$RPKC2	/	RC2TPMC.2
RC2TPM$TPMC3	<-	RC2TPM$RPKC3	/	RC2TPMC.3
RC2TPM$TPMC4	<-	RC2TPM$RPKC4	/	RC2TPMC.4
RC2TPM$TPMC5	<-	RC2TPM$RPKC5	/	RC2TPMC.5
RC2TPM$TPMC6	<-	RC2TPM$RPKC6	/	RC2TPMC.6
RC2TPM$TPMC7	<-	RC2TPM$RPKC7	/	RC2TPMC.7
RC2TPM$TPMC8	<-	RC2TPM$RPKC8	/	RC2TPMC.8
RC2TPM$TPMC9	<-	RC2TPM$RPKC9	/	RC2TPMC.9
RC2TPM$TPMC10	<-	RC2TPM$RPKC10	/	RC2TPMC.10
RC2TPM$TPMC11	<-	RC2TPM$RPKC11	/	RC2TPMC.11
RC2TPM$TPMC12	<-	RC2TPM$RPKC12	/	RC2TPMC.12
RC2TPM$TPMC13	<-	RC2TPM$RPKC13	/	RC2TPMC.13
RC2TPM$TPMC14	<-	RC2TPM$RPKC14	/	RC2TPMC.14
RC2TPM$TPMC15	<-	RC2TPM$RPKC15	/	RC2TPMC.15
RC2TPM$TPMC16	<-	RC2TPM$RPKC16	/	RC2TPMC.16
RC2TPM$TPMC17	<-	RC2TPM$RPKC17	/	RC2TPMC.17
RC2TPM$TPMC18	<-	RC2TPM$RPKC18	/	RC2TPMC.18
RC2TPM$TPMT1	<-	RC2TPM$RPKT1	/	RC2TPM.1
RC2TPM$TPMT2	<-	RC2TPM$RPKT2	/	RC2TPM.2
RC2TPM$TPMT3	<-	RC2TPM$RPKT3	/	RC2TPM.3
RC2TPM$TPMT4	<-	RC2TPM$RPKT4	/	RC2TPM.4
RC2TPM$TPMT5	<-	RC2TPM$RPKT5	/	RC2TPM.5
RC2TPM$TPMT6	<-	RC2TPM$RPKT6	/	RC2TPM.6
RC2TPM$TPMT7	<-	RC2TPM$RPKT7	/	RC2TPM.7
RC2TPM$TPMT8	<-	RC2TPM$RPKT8	/	RC2TPM.8
RC2TPM$TPMT9	<-	RC2TPM$RPKT9	/	RC2TPM.9
RC2TPM$TPMT10	<-	RC2TPM$RPKT10	/	RC2TPM.10
RC2TPM$TPMT11	<-	RC2TPM$RPKT11	/	RC2TPM.11
RC2TPM$TPMT12	<-	RC2TPM$RPKT12	/	RC2TPM.12
RC2TPM$TPMT13	<-	RC2TPM$RPKT13	/	RC2TPM.13
RC2TPM$TPMT14	<-	RC2TPM$RPKT14	/	RC2TPM.14
RC2TPM$TPMT15	<-	RC2TPM$RPKT15	/	RC2TPM.15
RC2TPM$TPMT16	<-	RC2TPM$RPKT16	/	RC2TPM.16

########################## CONFIRM TPM NORMALIZATION ###########################
if(sum(RC2TPM$TPMC1)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC1)}
if(sum(RC2TPM$TPMC2)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC2)}
if(sum(RC2TPM$TPMC3)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC3)}
if(sum(RC2TPM$TPMC4)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC4)}
if(sum(RC2TPM$TPMC5)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC5)}
if(sum(RC2TPM$TPMC6)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC6)}
if(sum(RC2TPM$TPMC7)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC7)}
if(sum(RC2TPM$TPMC8)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC8)}
if(sum(RC2TPM$TPMC9)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC9)}
if(sum(RC2TPM$TPMC10)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC10)}
if(sum(RC2TPM$TPMC11)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC11)}
if(sum(RC2TPM$TPMC12)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC12)}
if(sum(RC2TPM$TPMC13)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC13)}
if(sum(RC2TPM$TPMC14)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC14)}
if(sum(RC2TPM$TPMC15)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC15)}
if(sum(RC2TPM$TPMC16)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC16)}
if(sum(RC2TPM$TPMC17)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC17)}
if(sum(RC2TPM$TPMC18)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC18)}
if(sum(RC2TPM$TPMT1)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT1)}
if(sum(RC2TPM$TPMT2)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT2)}
if(sum(RC2TPM$TPMT3)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT3)}
if(sum(RC2TPM$TPMT4)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT4)}
if(sum(RC2TPM$TPMT5)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT5)}
if(sum(RC2TPM$TPMT6)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT6)}
if(sum(RC2TPM$TPMT7)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT7)}
if(sum(RC2TPM$TPMT8)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT8)}
if(sum(RC2TPM$TPMT9)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT9)}
if(sum(RC2TPM$TPMT10)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT10)}
if(sum(RC2TPM$TPMT11)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT11)}
if(sum(RC2TPM$TPMT12)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT12)}
if(sum(RC2TPM$TPMT13)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT13)}
if(sum(RC2TPM$TPMT14)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT14)}
if(sum(RC2TPM$TPMT15)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT15)}
if(sum(RC2TPM$TPMT16)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT16)}

write.csv(RC2TPM,"Results_SRA_CHOL_lncRNA\\1.Read_counts_to_TPM_SRA_CHOL_lncRNA.csv", row.names = TRUE)

#Making_Entire gene list without filter_DESeq result + TPM Values Control and Tumor
x <- RC2TPM
dim(x)
#[1] 127802    103
x[1:5,1:3]

x.heat <- x[70:103]############################################################
dim(x.heat)
#[1] 127802     34
write.csv(x.heat,"Results_SRA_CHOL_lncRNA\\2.Only_TPM_SRA_CHOL_lncRNA.csv", row.names = TRUE)
x.heat<-read.table("Results_SRA_CHOL_lncRNA\\2.Only_TPM_SRA_CHOL_lncRNA.csv", header = TRUE, sep = ",")
colnames(x.heat)[1]  <- "lncRNASymbol"    # change column name for x column

##              ##              ##              ##              ##              ##              ##              ##              ##
DESeq2_value<-read.table("DESeq2_Tumour_vs_Normal_SRA_CHOL_lncRNA.csv", header = TRUE, sep = ",")
##              ##              ##              ##              ##              ##              ##              ##              ##
colnames(DESeq2_value)[1]  <- "lncRNASymbol"    # change column name for x column
dim(DESeq2_value)
#[1] 127802      7
x.1 = merge(DESeq2_value, x.heat, ID="lncRNASymbol", sort = F)

# Average Control value as a new Variable (mean of Columns 8 to 12 is print in last column = 52)
x.1 <- transform(x.1, Mean_Contol_TPM = rowMeans(x.1[8:25], na.rm = TRUE))  #########################################################
# Average Control value as a new Variable (mean of Columns 13 to 41 is print in last column = 53)
x.1 <- transform(x.1, Mean_Tumour_TPM = rowMeans(x.1[26:41], na.rm = TRUE)) #########################################################

# Fold change (Tumour/average control)  (Fold change of Tumor from 43 to 59)
x.1 <- data.frame(x.1)                       # Creating example data
x.1$FC1  <- x.1$TPMT1  / x.1$Mean_Contol_TPM     # Creating new variable TPM1
x.1$FC2  <- x.1$TPMT2  / x.1$Mean_Contol_TPM     # Creating new variable TPM2
x.1$FC3  <- x.1$TPMT3  / x.1$Mean_Contol_TPM     # Creating new variable TPM3
x.1$FC4  <- x.1$TPMT4  / x.1$Mean_Contol_TPM     # Creating new variable TPM4
x.1$FC5  <- x.1$TPMT5  / x.1$Mean_Contol_TPM     # Creating new variable TPM5
x.1$FC6  <- x.1$TPMT6  / x.1$Mean_Contol_TPM     # Creating new variable TPM6
x.1$FC7  <- x.1$TPMT7  / x.1$Mean_Contol_TPM     # Creating new variable TPM7
x.1$FC8  <- x.1$TPMT8  / x.1$Mean_Contol_TPM     # Creating new variable TPM8
x.1$FC9  <- x.1$TPMT9  / x.1$Mean_Contol_TPM     # Creating new variable TPM9
x.1$FC10 <- x.1$TPMT10 / x.1$Mean_Contol_TPM # Creating new variable TPM10
x.1$FC11 <- x.1$TPMT11 / x.1$Mean_Contol_TPM   # Creating new variable TPM11
x.1$FC12 <- x.1$TPMT12 / x.1$Mean_Contol_TPM   # Creating new variable TPM12
x.1$FC13 <- x.1$TPMT13 / x.1$Mean_Contol_TPM   # Creating new variable TPM13
x.1$FC14 <- x.1$TPMT14 / x.1$Mean_Contol_TPM   # Creating new variable TPM14
x.1$FC15 <- x.1$TPMT15 / x.1$Mean_Contol_TPM   # Creating new variable TPM15
x.1$FC16 <- x.1$TPMT16 / x.1$Mean_Contol_TPM   # Creating new variable TPM16

write.csv(x.1,"Results_SRA_CHOL_lncRNA\\3.L2FC_SRA_CHOL_lncRNA.csv", row.names = TRUE)
x.1<-read.table("Results_SRA_CHOL_lncRNA\\3.L2FC_SRA_CHOL_lncRNA.csv", header = TRUE, sep = ",", row.names = 2)
x.1[1:5,1:6]
x.1 <- x.1[ -c(1) ] # Removes entire column be careful
x.1[1:5,1:6]
#log2 transform data
x.log <- log(x.1[43:58], 2)#################"x.log must have 16 Variables"####################################
x.log[x.log=="-Inf"]<-NA
colnames(x.log) <- c('L2FC1','L2FC2','L2FC3','L2FC4','L2FC5','L2FC6','L2FC7','L2FC8','L2FC9','L2FC10','L2FC11','L2FC12','L2FC13','L2FC14','L2FC15','L2FC16') #Rename the columns to L2FC
head(x.log)
write.csv(x.log,"Results_SRA_CHOL_lncRNA\\4.L2FC_only_SRA_CHOL_lncRNA.csv", row.names = TRUE)
z<-read.table("Results_SRA_CHOL_lncRNA\\4.L2FC_only_SRA_CHOL_lncRNA.csv", header = TRUE, sep = ",")
z[1:5,1:6]
colnames(z)[1]  <- "lncRNASymbol"    # change column name for x column
z[1:5,1:6]
y<-read.table("Results_SRA_CHOL_lncRNA\\3.L2FC_SRA_CHOL_lncRNA.csv", header = TRUE, sep = ",")
y[1:5,1:6]
y <- y[-c(1)]
y[1:5,1:6]
#******************************************************************************#
All_facotrs_merged_norow = merge(y, z, ID="lncRNASymbol", sort = F)            #
All_facotrs_merged <- All_facotrs_merged_norow[,-1]                            #
rownames(All_facotrs_merged) <- All_facotrs_merged_norow[,1]                   #
dim(All_facotrs_merged)                                                        #
#[1] 127802    74                                                             #
#******************************************************************************#
write.csv(All_facotrs_merged,"Results_SRA_CHOL_lncRNA\\5.All_factors_merged_TPM_CHOL_lncRNA.csv", row.names = TRUE)
#
#
#
#
## select features for clustering - for unsupervised analyses I often use highly expressed and variable.
l2fc.mean<-apply(x.log,1,mean,na.rm=TRUE)
l2fc <- round(l2fc.mean, digits = 2)
x.pvalue <- x.1[5] #########################################################
x.p <- round(x.pvalue, digits = 2)
x.pv<-apply(x.p,1,mean,na.rm=TRUE)
#print(x.pv)
#l2fc.pvalue<-apply(x.pvalue,1,sd,na.rm=TRUE)
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\1_CHOL_T_SRA_TPM_lncRNA.pdf")
plot(l2fc,x.pv,pch=".",ylab="pvalue",xlab="Log2 Fold change Mean Expression per lncRNA",main="CHOL_SRA_TPM_lncRNA")
dev.off()
#
#
#
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
                                     ##   UPREGULATED AND SIGNIFICANT lncRNA  ##                                        
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
#
#
#
#identify highly expressed and Significant genes
summary(l2fc>2 & x.pv<0.055)
#      Mode   FALSE    TRUE    NA's 
#   logical  124606    3196    5069 
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\2_CHOL_T_SRA_TPM_lncRNA.pdf")
plot(l2fc,x.pv,pch=".",main="CHOL_SRA_TPM_lncRNA_filtered",ylab="pvalue",xlab="Log2 Fold change Mean Expression per lncRNA",
     col=ifelse(l2fc>1.999& x.pv<0.055,"red","black"))
abline(v=2,h=0.055)
text(7.5,1.0,"n=3,196 genes")
dev.off()
#df2 <- data.frame(l2fc,x.pv)
#genes<-df2[df2$l2fc>2 & df2$x.pv<0.055,]  
#geneList<-na.omit(genes)
#dim(geneList)
#[1] 1580    2
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
#####################################   SELECTING ONLY GENES WITH L2FC > 2 EXPRESSION IN ALL SAMPLES   ###################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
x.l2fc_2g<- x.log[rowSums(x.log >= 2) == ncol(x.log), ]################### L2FC values >= 2 in all samples.###################
df <- data.frame(x.l2fc_2g)
#geneshortList<-df[df$x.l2fc_2g,]  
geneshortList<-na.omit(df)
dim(geneshortList)
#[1] 20   16
geneshortList[1:5,1:6]
## reduce matrix to list of genes that are both variable and expressed
x.heat1 <- data.frame(x.heat)
x.heat1 <- x.heat[,-1]
rownames(x.heat1) <- x.heat[,1]
l2fc1ly_TPM<-x.heat1[rownames(x.heat1) %in% rownames(geneshortList),]
dim(l2fc1ly_TPM)
#[1]  20   34
write.csv(l2fc1ly_TPM,"Results_SRA_CHOL_lncRNA\\6.l2fc1ly_TPM_SRA_CHOL_lncRNA.csv", row.names = TRUE)
forheatmap1<-read.table("Results_SRA_CHOL_lncRNA\\6.l2fc1ly_TPM_SRA_CHOL_lncRNA.csv", header = TRUE, sep = ",")
colnames(forheatmap1)[1]  <- "lncRNASymbol"    # change column name for x column
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
############################################## SELECTING ONLY GENES WITH pvalue <= 0.055 ####################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
## 
library(dplyr)
pvalue_filter <- All_facotrs_merged %>% filter(pvalue <= 0.055) 
dim(pvalue_filter)
#[1] 39240    74
write.csv(pvalue_filter,"Results_SRA_CHOL_lncRNA\\7.pvalue_sig_upr_CHOL_SRA_lncRNA.csv", row.names = TRUE)
#**********************************************************************************#
pvalue_0.055_filtered_TPM <- x.heat1[rownames(x.heat1) %in% rownames(pvalue_filter),]   #
dim(pvalue_0.055_filtered_TPM)                                                        #
#[1] 39240    34  *****************************************************************#
#**********************************************************************************#
write.csv(pvalue_0.055_filtered_TPM,"Results_SRA_CHOL_lncRNA\\7.pvalue_sig_upr_CHOL_SRA_lncRNA_TPM1ly.csv", row.names = TRUE)
#*#****************************************************************************************************************###########################################################################################################################################
##########################################################################################################################################
########################### SELECTING ONLY GENES WITH pvalue <= 0.055 and L2FC > 2 EXPRESSION IN ALL SAMPLES ################################
##########################################################################################################################################
##########################################################################################################################################
pvalue_0.055_L2FC_2_filter_intermediate <- pvalue_filter[c(59:74)]
dim(pvalue_0.055_L2FC_2_filter_intermediate)
#[1] 39240    16
pvalue_0.055_L2FC_2_filtered <- pvalue_0.055_L2FC_2_filter_intermediate[rowSums(pvalue_0.055_L2FC_2_filter_intermediate >= 2) == ncol(pvalue_0.055_L2FC_2_filter_intermediate), ]################ L2FC values >= 2 in all samples.#################
pvalue_0.055_L2FC_2_filtered <- na.omit(pvalue_0.055_L2FC_2_filtered)
dim(pvalue_0.055_L2FC_2_filtered)
#[1]  20    16
write.csv(pvalue_0.055_L2FC_2_filtered,"Results_SRA_CHOL_lncRNA\\8.pv_l2fc_upr_CHOL_SRA_lncRNA.csv", row.names = TRUE)
#**********************************************************************************************************#
pvalue_0.055_L2FC_2_filtered_TPM <- x.heat1[rownames(x.heat1) %in% rownames(pvalue_0.055_L2FC_2_filtered),]      #
dim(pvalue_0.055_L2FC_2_filtered_TPM)                                                                         #
#[1] 20  34  *********************************************************************************************#
#**********************************************************************************************************#
write.csv(pvalue_0.055_L2FC_2_filtered_TPM,"Results_SRA_CHOL_lncRNA\\8.pv_l2fc_upr_CHOL_SRA_lncRNA_TPM1ly.csv", row.names = TRUE)
#**********************************************************************************************************#
pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged <- All_facotrs_merged[rownames(All_facotrs_merged) %in% rownames(pvalue_0.055_L2FC_2_filtered),]
dim(pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged)
#[1]  20 74
write.csv(pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged,"Results_SRA_CHOL_lncRNA\\8.pvalue_pv_l2fc_upr_CHOL_SRA_lncRNA_AF.csv", row.names = TRUE)
#****************************************************************************************************************#
##########################################################################################################################################
##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.055 and L2FC > 2 and Mean_Tumour_TPM >= 1 EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
pvalue_0.055_l2fc_2_TPM_intermidiate_filtered <- All_facotrs_merged[rownames(All_facotrs_merged) %in% rownames(pvalue_0.055_L2FC_2_filtered_TPM),]
dim(pvalue_0.055_l2fc_2_TPM_intermidiate_filtered)
#[1]20  74
TumourTPM_filter_input <- pvalue_0.055_l2fc_2_TPM_intermidiate_filtered[c(25:40)]
pvalue_0.055_l2fc2_TPM_1_filtered <- TumourTPM_filter_input[rowSums(TumourTPM_filter_input >= 1) == ncol(TumourTPM_filter_input), ] 
dim(pvalue_0.055_l2fc2_TPM_1_filtered)
#[1]  9  16
write.csv(pvalue_0.055_l2fc2_TPM_1_filtered,"Results_SRA_CHOL_lncRNA\\9.pv_l2fc_TPM_upr_CHOL_SRA_lncRNA.csv", row.names = TRUE)
#****************************************************************************************************************#
pvalue_0.055_l2fc2_TPM_1_filtered_TPM <- x.heat1[rownames(x.heat1) %in% rownames(pvalue_0.055_l2fc2_TPM_1_filtered),]  #
dim(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)                                                                          #
#[1]   9  34 ****************************************************************************************************#
#****************************************************************************************************************#
write.csv(pvalue_0.055_l2fc2_TPM_1_filtered_TPM,"Results_SRA_CHOL_lncRNA\\9.pv_l2fc_TPM_upr_CHOL_SRA_lncRNA_TPM1ly.csv", row.names = TRUE)
##########################################################################################################################################
pvalue_0.055_l2fc2_TPM_1_filtered_TPM_all_factors_merged <- All_facotrs_merged[rownames(All_facotrs_merged) %in% rownames(pvalue_0.055_l2fc2_TPM_1_filtered),]
dim(pvalue_0.055_l2fc2_TPM_1_filtered_TPM_all_factors_merged)
#[1]  9 74
write.csv(pvalue_0.055_l2fc2_TPM_1_filtered_TPM_all_factors_merged,"Results_SRA_CHOL_lncRNA\\9.pvalue_pv_l2fc_TPM_upr_CHOL_SRA_lncRNA_AF.csv", row.names = TRUE)
#****************************************************************************************************************#
#*#****************************************************************************************************************#
ControlTPM_filter_input <- pvalue_0.055_l2fc2_TPM_1_filtered_TPM[c(1:18)]
pvalue_0.055_l2fc2_TPM_1_filtered_control <- ControlTPM_filter_input[rowSums(ControlTPM_filter_input <= 0.1) == ncol(ControlTPM_filter_input), ] 
dim(pvalue_0.055_l2fc2_TPM_1_filtered_control)
#[1]  0  18
write.csv(pvalue_0.055_l2fc2_TPM_1_filtered_control,"Results_SRA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_SRA_lncRNA.csv", row.names = TRUE)
#****************************************************************************************************************#
#****************************************************************************************************************#
pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM <- pvalue_0.055_l2fc2_TPM_1_filtered_TPM[rownames(pvalue_0.055_l2fc2_TPM_1_filtered_TPM) %in% rownames(pvalue_0.055_l2fc2_TPM_1_filtered_control),]  #
dim(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM)                                                                          #
#[1]   0  34 ****************************************************************************************************#
#****************************************************************************************************************#
write.csv(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM,"Results_SRA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_SRA_lncRNA_TPM1ly.csv", row.names = TRUE)
#****************************************************************************************************************#
pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_all_factors_merged <- All_facotrs_merged[rownames(All_facotrs_merged) %in% rownames(pvalue_0.055_l2fc2_TPM_1_filtered_control),]
dim(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_all_factors_merged)
#[1]  0 74
write.csv(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_all_factors_merged,"Results_SRA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_SRA_lncRNA_AF.csv", row.names = TRUE)
#****************************************************************************************************************#
#*#****************************************************************************************************************#
##########################################################################################################################################
#####################################   SELECTING ONLY GENES WITH L2FC > 2 EXPRESSION IN ALL SAMPLES   ###################################
##########################################################################################################################################
##########################################################################################################################################
l2fc1ly<-x.1[rownames(x.1) %in% rownames(geneshortList),]
dim(l2fc1ly)
#[1] 20   58

### Filterdata based on pvalue and TPM value
library(dplyr)
l2fc1ly_filterdf <- l2fc1ly %>% filter(pvalue <= 0.055 & Mean_Tumour_TPM >= 1) 
dim(l2fc1ly_filterdf)
#[1] 16 58
pv_TPM <- data.frame(l2fc1ly_filterdf)
write.csv(pv_TPM,"Results_SRA_CHOL_lncRNA\\11.pv_TPM__upr_CHOL_SRA_lncRNA.csv", row.names = TRUE)
pv_TPM %>% filter_at(vars(6), any_vars(. < 0.055)) 
pv_TPM.heatmap <- pv_TPM[7:40] 
##############***************&&&&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$#########################@@@@@@@@@@@@@@@@@@@@@@@
##############***************&&&&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$#########################@@@@@@@@@@@@@@@@@@@@@@@
##############***************&&&&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$#########################@@@@@@@@@@@@@@@@@@@@@@@
#impute missing data
library(impute) #bioconductor package  - https://www.bioconductor.org/packages/release/bioc/html/impute.html
pv.heatmap.imp<-impute.knn(as.matrix(pvalue_0.055_filtered_TPM))$data 
pv_l2fc.heatmap.imp<-impute.knn(as.matrix(pvalue_0.055_L2FC_2_filtered_TPM))$data 
pv_l2fc_TPM.heatmap.imp<-impute.knn(as.matrix(pvalue_0.055_l2fc2_TPM_1_filtered_TPM))$data 
pv_TPM.heatmap.imp<-impute.knn(as.matrix(pv_TPM.heatmap))$data            #### $data    add if see Error in data.imputed[1:features,]: incorrect number of dimensions
pv_l2fc_CnT_TPM.heatmap.imp<-impute.knn(as.matrix(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM))$data 
#
write.csv(pv_TPM.heatmap.imp,"Results_SRA_CHOL_lncRNA\\12.forheatmap_CHOL_SRA_lncRNA.csv", row.names = TRUE, col.names = TRUE)
#
##################################################
#Plot Up-regulated heatmap with column clustering#
source("arrayTools-noClanc.R")                   #
#####################################
# Basic Heatmap                     #
library(gplots)                     #
library(pheatmap)                   #
##################################################
gc()

pv_TPM <- tibble(pv.heatmap.imp)
pv_TPM <- data.frame(pv_TPM)
write.csv(pv_TPM,"Results_SRA_CHOL_lncRNA\\13.for_heatmap_lncRNA.csv", row.names = TRUE)

pv.heatmap<-read.table("Results_SRA_CHOL_lncRNA\\12.forheatmap_CHOL_SRA_lncRNA.csv", header = TRUE, sep = ",")
colnames(pv.heatmap)[1]  <- "lncRNASymbol"    # change column name for x column
###############################################################################################################################################
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\3_highexppv_TPMHeatmap.pdf")
pv_TPM_filtered <- data.frame(pv.heatmap$TPMC1,pv.heatmap$TPMC2,pv.heatmap$TPMC3,pv.heatmap$TPMC4,pv.heatmap$TPMC5,pv.heatmap$TPMC6,pv.heatmap$TPMC7,pv.heatmap$TPMC8,pv.heatmap$TPMC9,pv.heatmap$TPMC10,pv.heatmap$TPMC11,pv.heatmap$TPMC12,pv.heatmap$TPMC13,pv.heatmap$TPMC14,pv.heatmap$TPMC15,pv.heatmap$TPMC16,pv.heatmap$TPMC17,pv.heatmap$TPMC18,pv.heatmap$TPMT1,pv.heatmap$TPMT2,pv.heatmap$TPMT3,pv.heatmap$TPMT4,pv.heatmap$TPMT5,pv.heatmap$TPMT6,pv.heatmap$TPMT7,pv.heatmap$TPMT8,pv.heatmap$TPMT9,pv.heatmap$TPMT10,pv.heatmap$TPMT11,pv.heatmap$TPMT12,pv.heatmap$TPMT13,pv.heatmap$TPMT14,pv.heatmap$TPMT15,pv.heatmap$TPMT16)
colnames(pv_TPM_filtered) <- c("TPMC1","TPMC2","TPMC3","TPMC4","TPMC5","TPMC6","TPMC7","TPMC8","TPMC9","TPMC10","TPMC11","TPMC12","TPMC13","TPMC14","TPMC15","TPMC16","TPMC17","TPMC18","TPMT1","TPMT2","TPMT3","TPMT4","TPMT5","TPMT6","TPMT7","TPMT8","TPMT9","TPMT10","TPMT11","TPMT12","TPMT13","TPMT14","TPMT15","TPMT16")
rownames(pv_TPM_filtered) <- pv.heatmap$lncRNASymbol
pv_TPM_filtered <- data.matrix(pv_TPM_filtered)
pheatmap(pv_TPM_filtered ,scale = "row", cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("green","red"))(50),fontsize =8)
dev.off()
################################################################################################################################################
#####################################
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\4_pvalueL2FC2pv_l2fcheat.pdf")
myHeatmap((pv_l2fc.heatmap.imp))   ############ plot Heatmap ##################
dev.off()
##################################################
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\5_pvaluel2fc2TPM1heat7.pdf")
myHeatmap(( pv_l2fc_TPM.heatmap.imp))     ############ plot Heatmap ##################
dev.off()
##################################################
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\6_pvalue_TPM1_pv_TPM_heat.pdf")
myHeatmap(( pv_TPM.heatmap.imp))  ############ plot Heatmap ##################
dev.off() 
##################################################
#pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\7_pvaluel2fc2TPM1heat0.pdf")
#myHeatmap(( pv_l2fc_CnT_TPM.heatmap.imp))     ############ plot Heatmap ##################
#dev.off()
##################################################

pv_TPM_TPM <- tibble(pv_TPM.heatmap.imp)
pv_TPM_TPM <- data.frame(pv_TPM_TPM)
write.csv(pv_TPM_TPM,"Results_SRA_CHOL_lncRNA\\13.for_heatmap_lncRNA.csv", row.names = TRUE)

#Plot Up-regulated heatmap without column clustering
pv_TPM.heatmap<-read.table("Results_SRA_CHOL_lncRNA\\12.forheatmap_CHOL_SRA_lncRNA.csv", header = TRUE, sep = ",")
colnames(pv_TPM.heatmap)[1]  <- "lncRNASymbol"    # change column name for x column
#pv_TPM.heatmap <- data.matrix(pv_TPM.heatmap)
#####################################
# Basic Heatmap                     #
library(gplots)                     #
library(pheatmap)                   #
#####################################
################################################################################################################################################
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\8_highexppv_TPMHeatmap.pdf")
pv_TPM_TPM_filtered <- data.frame(pv_TPM.heatmap$TPMC1,pv_TPM.heatmap$TPMC2,pv_TPM.heatmap$TPMC3,pv_TPM.heatmap$TPMC4,pv_TPM.heatmap$TPMC5,pv_TPM.heatmap$TPMC6,pv_TPM.heatmap$TPMC7,pv_TPM.heatmap$TPMC8,pv_TPM.heatmap$TPMC9,pv_TPM.heatmap$TPMC10,pv_TPM.heatmap$TPMC11,pv_TPM.heatmap$TPMC12,pv_TPM.heatmap$TPMC13,pv_TPM.heatmap$TPMC14,pv_TPM.heatmap$TPMC15,pv_TPM.heatmap$TPMC16,pv_TPM.heatmap$TPMC17,pv_TPM.heatmap$TPMC18,pv_TPM.heatmap$TPMT1,pv_TPM.heatmap$TPMT2,pv_TPM.heatmap$TPMT3,pv_TPM.heatmap$TPMT4,pv_TPM.heatmap$TPMT5,pv_TPM.heatmap$TPMT6,pv_TPM.heatmap$TPMT7,pv_TPM.heatmap$TPMT8,pv_TPM.heatmap$TPMT9,pv_TPM.heatmap$TPMT10,pv_TPM.heatmap$TPMT11,pv_TPM.heatmap$TPMT12,pv_TPM.heatmap$TPMT13,pv_TPM.heatmap$TPMT14,pv_TPM.heatmap$TPMT15,pv_TPM.heatmap$TPMT16)
colnames(pv_TPM_TPM_filtered) <- c("TPMC1","TPMC2","TPMC3","TPMC4","TPMC5","TPMC6","TPMC7","TPMC8","TPMC9","TPMC10","TPMC11","TPMC12","TPMC13","TPMC14","TPMC15","TPMC16","TPMC17","TPMC18","TPMT1","TPMT2","TPMT3","TPMT4","TPMT5","TPMT6","TPMT7","TPMT8","TPMT9","TPMT10","TPMT11","TPMT12","TPMT13","TPMT14","TPMT15","TPMT16")
rownames(pv_TPM_TPM_filtered) <- pv_TPM.heatmap$lncRNASymbol
pv_TPM_TPM_filtered <- data.matrix(pv_TPM_TPM_filtered)
pheatmap(pv_TPM_TPM_filtered ,scale = "row", cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("green","red"))(50),fontsize =8)
dev.off()
################################################################################################################################################
#####################################
library(ComplexHeatmap) #3D Heatmap #
#####################################
################################################################################################################################################
max(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)
min(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)
#pvalue_0.055_l2fc2_TPM_1_filtered_TPM["lnc-TNFSF14-4:1", "TPMT21"] <- 50
max(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)
min(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)
View(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\9_highexp_sig7_3DHeat.pdf")
  pvalue_0.055_l2fc2_TPM_1_filtered_TPM <- data.matrix(pvalue_0.055_l2fc2_TPM_1_filtered_TPM)
  mat7 = pvalue_0.055_l2fc2_TPM_1_filtered_TPM
  Sample_Name = as.character(colnames(mat7))
  ha_column = HeatmapAnnotation(
    Sample = anno_text(Sample_Name, rot = 90, location = unit(0.5, "npc"), just = "top")
  )
  col_fun = circlize::colorRamp2(c(0, 1, 75, 150, 260), c( "#ffffff", "#00e1ff", "#ffdd00", "#fe9501", "#ff0000"))
  ht_opt$TITLE_PADDING = unit(c(15, 2), "mm") 
  Heatmap3D(mat7, name = "TPM", col = col_fun,
            cluster_columns = FALSE, show_row_dend = FALSE, 
            show_column_names = FALSE,
            row_names_side = "left", row_names_gp = gpar(fontsize = 4),
            column_title = 'pvalue<0.055 _ l2fc>2 _ TPM>1 _ CHOL SRA',
            bottom_annotation = ha_column,
            heatmap_legend_param = list(at = c(0, 1, 75, 150, 260), 
                                        labels = c("0TPM", "1TPM", "75TPM", "150TPM", "260TPM")),
            # new arguments for Heatmap3D()
            bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm")
  )
dev.off()
################################################################################################################################################
################################################################################################################################################
  max(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM)
  min(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM)
  pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\10_highexp_sig23_3DHeat.pdf")
  pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM <- data.matrix(pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM)
  mat23 = pvalue_0.055_l2fc2_C_0.1_Tumour_TPM1_filtered_TPM
  Sample_Name = as.character(colnames(mat23))
  ha_column = HeatmapAnnotation(
    Sample = anno_text(Sample_Name, rot = 90, location = unit(0.5, "npc"), just = "top")
  )
  col_fun = circlize::colorRamp2(c(0, 1, 100, 500, 1164), c("white", "cornflowerblue", "yellow", "orange", "red"))
  ht_opt$TITLE_PADDING = unit(c(15, 2), "mm")
  Heatmap3D(mat23, name = "TPM", col = col_fun,
            cluster_columns = FALSE, show_row_dend = FALSE, 
            show_column_names = FALSE,
            row_names_side = "left", row_names_gp = gpar(fontsize = 4),
            column_title = 'pvalue<0.055 _ L2FC>2_ controlTPM<0.1 & TumourTPM>1 _ CHOL SRA',(fontsize = 3),
            bottom_annotation = ha_column,
            heatmap_legend_param = list(at = c(0, 1, 100, 500, 1164), 
                                        labels = c("0TPM", "1TPM", "100TPM", "500TPM", "1164TPM")),
            # new arguments for Heatmap3D()
            bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm")
  )
  dev.off()
#####################################
library(ComplexHeatmap) #3D Heatmap #
#####################################
max(pv_TPM_TPM_filtered)
min(pv_TPM_TPM_filtered)
View(pv_TPM_TPM_filtered)
#pv_TPM_TPM_filtered["lnc-TNFSF14-4:1", "TPMT21"] <- 50
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\11_hexp_pv_TPM_3DHeat.pdf")
pv_TPM_TPM_filtered <- data.matrix(pv_TPM_TPM_filtered)
pv_TPM = pv_TPM_TPM_filtered
Sample_Name = as.character(colnames(pv_TPM))
ha_column = HeatmapAnnotation(
  Sample = anno_text(Sample_Name, rot = 90, location = unit(0.5, "npc"), just = "top")
)
col_fun = circlize::colorRamp2(c(0, 1, 75, 150, 260), c( "#ffffff", "#00e1ff", "#ffdd00", "#fe9501", "#ff0000"))
ht_opt$TITLE_PADDING = unit(c(15, 2), "mm") 
Heatmap3D(pv_TPM, name = "TPM", col = col_fun,
          cluster_columns = FALSE, show_row_dend = FALSE, 
          show_column_names = FALSE,
          row_names_side = "left", row_names_gp = gpar(fontsize = 4),
          column_title = 'pvalue<0.055 _ TPM>1 _ CHOL SRA',
          bottom_annotation = ha_column,
          heatmap_legend_param = list(at = c(0, 1, 75, 150, 260), 
                                      labels = c("0TPM", "1TPM", "75TPM", "150TPM", "260TPM")),
          # new arguments for Heatmap3D()
          bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm")
)
dev.off()
################################################################################################################################################################################################################################################################################################
max(pvalue_0.055_L2FC_2_filtered_TPM)
min(pvalue_0.055_L2FC_2_filtered_TPM)
pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\12_highexp_pv_l2fc_3DHeat.pdf")
pvalue_0.055_L2FC_2_filtered_TPM <- data.matrix(pvalue_0.055_L2FC_2_filtered_TPM)
pv_l2fc = pvalue_0.055_L2FC_2_filtered_TPM
Sample_Name = as.character(colnames(pv_l2fc))
ha_column = HeatmapAnnotation(
  Sample = anno_text(Sample_Name, rot = 90, location = unit(0.5, "npc"), just = "top")
)
col_fun = circlize::colorRamp2(c(0, 1, 10, 16, 75, 260), c("white", "cornflowerblue", "yellow", "orange", "red", "grey"))
ht_opt$TITLE_PADDING = unit(c(15, 2), "mm")
Heatmap3D(pv_l2fc, name = "TPM", col = col_fun,
          cluster_columns = FALSE, show_row_dend = FALSE, 
          show_column_names = FALSE,
          row_names_side = "left", row_names_gp = gpar(fontsize = 4),
          column_title = 'pvalue<0.055 _ L2FC>2 _ CHOL SRA',(fontsize = 3),
          bottom_annotation = ha_column,
          heatmap_legend_param = list(at = c(0, 1, 10, 16, 75, 260), 
                                      labels = c("0TPM", "1TPM", "15TPM", "16TPM", "75TPM", "260TPM")),
          # new arguments for Heatmap3D()
          bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm")
)
dev.off()
################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################

pdf(file="Results_SRA_CHOL_lncRNA\\Figures_SRA_CHOL_lncRNA\\13_highexp_pv_l2fc_3DHeat.pdf")
pvalue_0.055_L2FC_2_filtered_TPM <- data.matrix(pvalue_0.055_L2FC_2_filtered_TPM)
pv_l2fc = pvalue_0.055_L2FC_2_filtered_TPM
Sample_Name = as.character(colnames(pv_l2fc))
ha_column = HeatmapAnnotation(
  Sample = anno_text(Sample_Name, rot = 90, location = unit(0.5, "npc"), just = "top")
)
col_fun = circlize::colorRamp2(c(0, 1, 10, 50, 100, 500, 1000, 5000), c("violet", "#4b0082", "blue", "green", "yellow", "orange", "red", "white"))
ht_opt$TITLE_PADDING = unit(c(15, 2), "mm")
Heatmap3D(pv_l2fc, name = "TPM", col = col_fun,
          cluster_columns = FALSE, show_row_dend = FALSE, 
          show_column_names = FALSE,
          row_names_side = "left", row_names_gp = gpar(fontsize = 4),
          column_title = 'pvalue<0.055 _ L2FC>2 _ CHOL SRA',(fontsize = 3),
          bottom_annotation = ha_column,
          heatmap_legend_param = list(at = c(0, 1, 10, 50, 100, 500, 1000, 5000), 
                                      labels = c("0TPM", "1TPM", "10TPM", "50TPM", "100TPM", "500TPM", "1000TPM", "5000TPM")),
          # new arguments for Heatmap3D()
          bar_rel_width = 1, bar_rel_height = 1, bar_max_length = unit(1, "cm")
)
dev.off()
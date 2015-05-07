Imprinting analysis for Baran et al. 

Tuuli Lappalainen, New York Genome Center
May 6th 2015

Note that this file includes only the steps to produce results and figures of the manuscript - many additional analyses were performed to ensure their robustness. 
Due to limitations in the public availability of GTEx data, most users will not be able to rerun the analysis as is. 
Thus, the main motivation for sharing the code is transcaprency of the analysis. 


###########################################################################################################################################
#                                                    Raw data processing and formatting                                                   #
###########################################################################################################################################



##########  Matlab statistic files  ##########


cd ~/tuuli_lab/tlappalainen/imprint/data/imprint/gtex
perl ~/imprinting/utilities/CombineResList_filter5_model15_pertissue.pl GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt gtex_alltissues.txt impglr_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr
(cat impglr_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.ADPSBQ | head -n 1 ; cat impglr_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.* | grep -v tissue) >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr

cd ~/tuuli_lab/tlappalainen/imprint/data/imprint/gtex/brainsub_res_proc
perl ~/imprinting/utilities/CombineResList_filter5_model15_pertissue_brain.pl GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.IMPR47.txt ../brainsub/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.IMPR47.txt.tis impglr_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.txt.impglr
(cat impglr_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.txt.impglr.BRNACC | head -n 1 ; cat impglr_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.txt.impglr.* ) | grep -v tissue >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.txt.impglr

cd ~/tuuli_lab/tlappalainen/imprint/data/imprint/gencord
perl ~/imprinting/utilities/CombineResList_filter5_model15_pertissue_nochr.pl GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt gencord_alltissues.txt impglr_proc/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr
(cat impglr_proc/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.LCL | head -n 1 ; cat impglr_proc/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.* | grep -v tissue) >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr

cd ~/tuuli_lab/tlappalainen/imprint/data/imprint/geuvadis
perl ~/imprinting/utilities/CombineResList_filter5_model15_pertissue_nochr.pl GD462.ASE.COV8.ANNOTPLUS.SINFO.txt geuvadis_alltissues.txt GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr



##########  Gene lists  ##########

cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint/

#remList has been modified to add:
#remAnn  RP11-701H24.3
#remAnn  AL132709.5
#remAnn  RP11-395B7.2
#remFam  NLRP2
#remLCL  HM13

#ensg2nametab.txt modified to change INPP5F to INPP5F_V2


data <- read.table("known.novel.classify.txt", as.is=T, header=T)
data[,4] <- replace(data[,4], data[,4]=="INPP5F", "INPP5F_V2")
write.table(data, "known.novel.classify.txt", quote=F, row.names=F, col.names=T, sep='\t')


final <- read.table("known.novel.classify.txt", as.is=T, fill=T, header=T)
rem <- read.table("remList.txt", as.is=T)
keep <- vector("integer")
for(i in 1:nrow(final)) {
	if (!any(rem[,2]==final[,4][i])) {
		keep <- c(keep, i)
	}
}
final <- final[keep,]
final[final[,1]=="../gdresBi",2] <- "GD_LCL" 
final[final[,1]=="../gencordresBi" & final[,2]=="LCL",2] <- "GC_LCL" 
final[final[,1]=="../gencordresBi" & final[,2]=="TCELL",2] <- "GC_TCELL" 
final[final[,1]=="../gencordresBi" & final[,2]=="FIBRBLS",2] <- "GC_FIBRBLS" 
final_imp_gene <- unique(final[,4][final[,5]=="-"])
keep <- vector("integer")
for(i in 1:nrow(final)) {
	if (any(final_imp_gene==final[,4][i])) {
		keep <- c(keep, i)
	}
}
final_imp <- final[keep,]
write.table(final_imp[,c(4,3,2)], "final_impr_tab.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(final_imp[,2:ncol(final)], "final_impr_bi_imp.txt", row.names=F, col.names=T, sep="\t", quote=F)
final_known <- final[final[,8]=="+" & final[,9]!="M",]
write.table(final_known[,c(4,3,2)], "final_known_tab.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(final_known[,2:ncol(final)], "final_known_tab_bi_imp.txt", row.names=F, col.names=T, sep="\t", quote=F)


cut -f 1 final_impr_tab.txt | sort -u  >final_impr_genelist.txt
cut -f 1,2 final_impr_tab.txt | sort -u | awk '{print($2 "\t" $1 )}' > final_impr_genelist_ensg2name.txt

cut -f 2 final_known_tab.txt | sort -u  >known_otago.txt
cut -f 1,2 final_known_tab.txt | sort -u | awk '{print($2 "\t" $1 )}' > known_otago_ensg2name.txt



#########  Genes of intrest  ##########

#genes_of_interest_name.txt
perl ~/utilities/intersect.pl genes_of_interest_name.txt 0 ensg2nametab.txt 1 >extra_ensg2name.txt
cat putgenes.txt extra_ensg2name.txt | sort -u >extra_put_ensg2name.txt
perl ~/utilities/setdiff.pl final_impr_genelist_ensg2name.txt 0 extra_known_ensg2name.txt 0 >extra_ensg2name_1.txt
perl ~/utilities/setdiff.pl known_otago_ensg2name.txt 0 extra_ensg2name_1.txt 0 >extra_ensg2name.txt


#########  all known and provisional genes in human  ##########

ot <- read.table("otago.list.txt", as.is=T)
otf <- ot[ot[,5]!="M",]
otf <- otf[otf[,2]!="-",]
og <- unique(otf[,2]) #133 genes
final <- read.table("final_impr_tab.txt", as.is=T)
#take only those that are not already included
og_sd <- setdiff(og, unique(final[,1]))
write.table(og_sd, "otago.list.human.genename.txt", quote=F, row.names=F, col.names=F)

perl ~/utilities/intersect.pl otago.list.human.genename.txt 0 ensg2nametab.txt 1 >otago.human.all_ensg2name.txt
perl ~/utilities/setdiff.pl final_impr_genelist_ensg2name.txt 0 otago.human.all_ensg2name.txt 0 >otago.human.bi_ensg2name.txt


#########  Matlab haplotype files  ##########


## GTEx

setwd("~/tuuli_lab/tlappalainen/imprint/data/imprint/gtex")

##only genotyped SNPs, not imputed
tis <- read.table("GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.tis", as.is=T)[,1]
data <- vector("list", length(tis))
for(i in 1:length(tis)) {
	tmp <- vector("list", 22)
	for(j in 1:22) {
		tmp[[j]] <- read.table(paste("indxgen/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt", tis[i], j, "filter5.model15.indxgen.txt", sep="."), header=T, as.is=T, fill=T)
	}
	data[[i]] <- do.call(rbind, tmp)
}
for(i in 1:length(tis)) {
	data[[i]] <- cbind(rep(tis[i], nrow(data[[i]])), data[[i]])
}
tot <- do.call(rbind, data)
colnames(tot)[1] <- "tissue"
tot2 <- tot[!is.na(tot[,ncol(tot)]),] ##one row with a weird number - checked, OK to ignore
write.table(tot2, "GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.indxgen.txt", sep="\t", quote=F, row.names=F)

##including imputed SNPs
tis <- read.table("GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.tis", as.is=T)[,1]
data <- vector("list", length(tis))
for(i in 1:length(tis)) {
	tmp <- vector("list", 22)
	for(j in 1:22) {
		tmp[[j]] <- read.table(paste("indxgen_imp/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt", tis[i], j, "filter5.model15.imp.indxgen.txt", sep="."), header=T, as.is=T, fill=T)
	}
	data[[i]] <- do.call(rbind, tmp)
}
for(i in 1:length(tis)) {
	data[[i]] <- cbind(rep(tis[i], nrow(data[[i]])), data[[i]])
}
tot <- do.call(rbind, data)
colnames(tot)[1] <- "tissue"
tot2 <- tot[!is.na(tot[,ncol(tot)]),] ##one row with a weird number - checked, OK to ignore
write.table(tot2, "GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.txt", sep="\t", quote=F, row.names=F)

#brain subtissues
tis <- read.table("~/tuuli_lab/tlappalainen/imprint/data/imprint/gtex/brainsub/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.IMPR47.txt.tis", as.is=T)[,1]
data <- vector("list", length(tis))
for(i in 1:length(tis)) {
	tmp <- vector("list", 22)
	for(j in 1:22) {
		if (file.exists(paste("brainsub_impglr/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.IMPR47.txt", tis[i], j, "filter5.model15.imp.indxgen.txt", sep="."))) {
			tmp[[j]] <- read.table(paste("brainsub_impglr/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.IMPR47.txt", tis[i], j, "filter5.model15.imp.indxgen.txt", sep="."), header=T, as.is=T, fill=T)
		}
	}
	data[[i]] <- do.call(rbind, tmp)
}
for(i in 1:length(tis)) {
	data[[i]] <- cbind(rep(tis[i], nrow(data[[i]])), data[[i]])
}
tot <- do.call(rbind, data)
colnames(tot)[1] <- "tissue"
tot2 <- tot[!is.na(tot[,ncol(tot)]),] ##one row with a weird number. ignore
write.table(tot2, "GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.IMPR47.txt.filter5.model15.imp.indxgen.txt", sep="\t", quote=F, row.names=F)


##Gencord

setwd("~/tuuli_lab/tlappalainen/imprint/data/imprint/gencord")

tis <- c("TCELL", "LCL", "FIBRBLS")
data <- vector("list", length(tis))
for(i in 1:length(tis)) {
	data[[i]] <- read.table(paste("indxgen/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt", tis[i], "filter5.model15.indxgen.txt", sep="."), header=T, as.is=T, fill=T)
}
for(i in 1:length(tis)) {
	data[[i]] <- cbind(rep(tis[i], nrow(data[[i]])), data[[i]])
}
tot <- do.call(rbind, data)
colnames(tot)[1] <- "tissue"
write.table(tot2, "GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt", sep="\t", quote=F, row.names=F)


##Geuvadis

setwd("~/tuuli_lab/tlappalainen/imprint/data/imprint/geuvadis")

tis <- c("LCL")
data <- vector("list", length(tis))
for(i in 1:length(tis)) {
	data[[i]] <- read.table(paste("indxgen/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt", tis[i], "filter5.model15.indxgen.txt", sep="."), header=T, as.is=T, fill=T)
}
for(i in 1:length(tis)) {
	data[[i]] <- cbind(rep(tis[i], nrow(data[[i]])), data[[i]])
}
tot <- do.call(rbind, data)
colnames(tot)[1] <- "tissue"
write.table(tot, "GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt", sep="\t", quote=F, row.names=F)







##########  Data of imprinted genes  ##########


cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint/gtex
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt 21 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.indxgen.txt 1 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.indxgen.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.txt 1 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr 1 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../known_otago_ensg2name.txt 0 ../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr 1 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.known.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../extra_ensg2name.txt 0 ../../../data/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt 21 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.extra.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../extra_ensg2name.txt 0 ../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.txt 1 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.extra.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../otago.human.bi_ensg2name.txt 0 ../../../data/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt 21 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.otago.human.bi.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../otago.human.bi_ensg2name.txt 0 ../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.txt 1 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.otago.human.bi.txt' | qsub

d <- read.table("GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.final_impr.txt", header=T, as.is=T)
d[,3] <- replace(d[,3], d[,3]=="INPP5F", "INPP5F_V2")
write.table(d, "GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.final_impr.txt", sep="\t", quote=F, row.names=F)

d <- read.table("GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.known.txt", header=T, as.is=T)
d[,3] <- replace(d[,3], d[,3]=="INPP5F", "INPP5F_V2")
write.table(d, "GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.known.txt", sep="\t", quote=F, row.names=F)


#retrieve only 11 tissues
data <- read.table("../../../data/imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr", sep="\t", header=T, as.is=T, fill=T)
tis <- c("ADPSBQ", "ARTTBL", "BRAIN", "HRTLV", "LCL", "LUNG", "MSCLSK", "NERVET", "SKINS", "THYROID", "WHLBLD")
l <- vector("list", length(tis))
for(i in 1:length(tis)) {
	l[[i]] <- data[data[,1]==tis[i],]
}
dataf <- do.call(rbind, l)
dataf2 <- dataf[dataf[,20]>=5 & dataf[,21]>=2,]
dataf2[,3] <- replace(dataf2[,3], dataf2[,3]=="INPP5F", "INPP5F_V2")
write.table(dataf2, "GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.passed", row.names=F, sep="\t", quote=F)
cut -f 2 GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.passed | sort -u >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.passed.genes

#brain regions
cd ~/tuuli_lab/tlappalainen/imprint/data/gtex

sh brain_subreg_impr.sh | head -n1
perl ~/Geneva/utilities/intersect_keepheader.pl ~/tuuli_lab/tlappalainen/imprint/analysis/imprint/final_impr_genelist_ensg2name.txt 0 ~/Geneva/data/imprint/data/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN_BRNACC 21 >~/Geneva/data/imprint/data/gtex/brsub/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN_BRNACC.impr.txt
##same for all subregions

cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint/geuvadis/
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt 21 >GD462.ASE.COV8.ANNOTPLUS.SINFO.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt 1 >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.LCL 1 >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../known_otago_ensg2name.txt 0 ../../../data/imprint/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.LCL 1 >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.known.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../extra_ensg2name.txt 0 ../../../data/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt 21 >GD462.ASE.COV8.ANNOTPLUS.SINFO.extra.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../extra_ensg2name.txt 0 ../../../data/imprint/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt 1 >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.extra.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../otago.human.bi_ensg2name.txt 0 ../../../data/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt 21 >GD462.ASE.COV8.ANNOTPLUS.SINFO.otago.human.bi.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../otago.human.bi_ensg2name.txt 0 ../../../data/imprint/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt 1 >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.otago.human.bi.txt' | qsub

awk '$20>=5 && $21>=2 ' ../../../data/imprint/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.LCL >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.LCL.passed
cut -f 2 GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.LCL.passed | sort -u >GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.LCL.passed.genes

echo 'perl ~/utilities/intersect_keepheader.pl ../failed_examples_ensg2name.txt 0 ../../../data/geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt 21 >GD462.ASE.COV8.ANNOTPLUS.SINFO.failed_examples.txt' | qsub

d <- read.table("GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", header=T, as.is=T)
d[,3] <- replace(d[,3], d[,3]=="INPP5F", "INPP5F_V2")
write.table(d, "GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", sep="\t", quote=F, row.names=F)

d <- read.table("GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.known.txt", header=T, as.is=T)
d[,3] <- replace(d[,3], d[,3]=="INPP5F", "INPP5F_V2")
write.table(d, "GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.known.txt", sep="\t", quote=F, row.names=F)


cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint/gencord/
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt 21 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt 1 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../final_impr_genelist_ensg2name.txt 0 ../../../data/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr 1 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../known_otago_ensg2name.txt 0 ../../../data/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr 1 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.known.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../extra_ensg2name.txt 0 ../../../data/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt 21 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.extra.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../extra_ensg2name.txt 0 ../../../data/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt 1 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.extra.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../otago.human.bi_ensg2name.txt 0 ../../../data/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt 21 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.otago.human.bi.txt' | qsub
echo 'perl ~/utilities/intersect_keepheader.pl ../otago.human.bi_ensg2name.txt 0 ../../../data/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.txt 1 >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.otago.human.bi.txt' | qsub

awk '$20>=5 && $21>=2 ' ../../../data/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.passed
cut -f 2 GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.passed | sort -u >GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.passed.genes

d <- read.table("GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", header=T, as.is=T)
d[,3] <- replace(d[,3], d[,3]=="INPP5F", "INPP5F_V2")
write.table(d, "GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", sep="\t", quote=F, row.names=F)

d <- read.table("GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.known.txt", header=T, as.is=T)
d[,3] <- replace(d[,3], d[,3]=="INPP5F", "INPP5F_V2")
write.table(d, "GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.known.txt", sep="\t", quote=F, row.names=F)







##########  Biallelic/imprinted classifications  ##########


cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint

data <- read.table("final_impr_bi_imp.txt", header=T, as.is=T)
tis <- unique(data[,1])
gene <- unique(data[,3])
tab <- matrix(length(gene), length(tis), data=NA)
tab2 <- matrix(length(gene), length(tis), data=NA)
rownames(tab) <- gene
colnames(tab) <- tis
rownames(tab2) <- gene
colnames(tab2) <- tis
for(i in 1:length(tis)) {
	for(j in 1:length(gene)) {
		tmp <- data[data[,1]==tis[i] & data[,3]==gene[j],]
		if (nrow(tmp)>0) {
			tab[gene[j], tis[i]] <- tmp[1,4]
			tab2[gene[j], tis[i]] <- tmp[1,5]
		}
	}
}
write.table(tab, "biallelic.strong.final.tab.txt", sep="\t", quote=F)
write.table(tab2, "biallelic.weak.final.tab.txt", sep="\t", quote=F)


#known genes
data <- read.table("final_known_tab_bi_imp.txt", header=T, as.is=T)
tis <- unique(data[,1])
gene <- unique(data[,3])
tab <- matrix(length(gene), length(tis), data=NA)
tab2 <- matrix(length(gene), length(tis), data=NA)
rownames(tab) <- gene
colnames(tab) <- tis
rownames(tab2) <- gene
colnames(tab2) <- tis
for(i in 1:length(tis)) {
	for(j in 1:length(gene)) {
		tmp <- data[data[,1]==tis[i] & data[,3]==gene[j],]
		if (nrow(tmp)>0) {
			tab[gene[j], tis[i]] <- tmp[1,4]
			tab2[gene[j], tis[i]] <- tmp[1,5]
		}
	}
}
write.table(tab, "biallelic.strong.known.tab.txt", sep="\t", quote=F)
write.table(tab2, "biallelic.weak.known.tab.txt", sep="\t", quote=F)








###########################################################################################################################################
#                                                              Clonality                                                                  #
###########################################################################################################################################


#statistics of total and biallelic ASE site counts
gtex <- read.table("chrX/ase_stats_chrX_8_05_gtex.txt", as.is=T, header=T, sep="\t")
gd <- read.table("chrX/ase_stats_chrX_8_05_gd.txt", as.is=T, header=T, sep="\t")
gc <- read.table("chrX/ase_stats_chrX_16_gc.txt", as.is=T, header=T, sep="\t")
rownames(gtex) <- gtex[,1]
rownames(gd) <- gd[,1]
rownames(gc) <- gc[,1]

info <- read.table("../../data/gtex/GTEx.SubjectPhenotypesDS.v9.1.SUBJID_GENDER_AGE.txt", header=T, as.is=T, sep="\t")
fem <- info[info[,2]==2,]
gtexinfo <- read.table("~/tuuli_lab/tlappalainen/gtex/ug_analysis/ase/imputed/GTExDecReleaseFreezeAseImp1582.Sampleinfo.StandardNames.Colors.20130815.txt", as.is=T, header=T, sep="\t")
rownames(gtexinfo) <- gtexinfo[,1]
int <- intersect(gtexinfo[,1], gtex[,1])
gtex <- gtex[int,]
gtexinfo <- gtexinfo[int,]
tvec <- gtexinfo$TISSUE_ABBRV
tvec[grep("BRN", tvec)] <- "BRAIN"
sub <- unlist(lapply(strsplit(gtex[,1], "-"), function(x) {paste(x[1], x[2], sep="-")}))
log <- rep(0, length(sub))
for(i in 1:length(sub)) { 
	if (any(fem[,1]==sub[i])) {
		log[i] <- 1
	}
}
tvecf <- tvec[log==1]
gtp <- (gtex[,3]/gtex[,2])[log==1]
tab <- cbind(gtex[log==1,3], gtex[log==1,2], gtp)
colnames(tab) <- c("RNAseq_het_sites", "RNAseq_total_sites", "Proportion")
rownames(tab) <- rownames(gtex)[log==1]
write.table(tab, "GTEx_females_chrX_RNAseqhet_over_totalhetsites.txt", col.names=F, quote=F)

sub <- unlist(lapply(strsplit(gd[,1], ".", fixed=T), function(x) {x[1]}))
gdsex <- read.table("~/tuuli_lab/tlappalainen/geuvadis/data/qc_stats/GD462.sex.txt", as.is=T)
log <- rep(0, length(sub))
for(i in 1:length(sub)) { 
	if (any(gdsex[gdsex[,1]==sub[i],2]==2)) {
		log[i] <- 1
	}
}
gdp <- (gd[,3]/gd[,2])[log==1]
gdvec <- rep("LCL", length(gdp))

het <- read.table("GENCORD2.chrX.het", header=T, as.is=T)
sub <- sub("F", "", gc[,1])
sub <- sub("T", "", sub)
sub <- sub("L", "", sub)
sub <- sub("B", "", sub)
log <- rep(0, length(sub))
for(i in 1:length(sub)) { 
	if (het[het[,1]==sub[i],2]<260000) {
		log[i] <- 1
	}
}
gcp <- (gc[,3]/gc[,2])[log==1]
name <- gc[,1]
name <- sub("UC", "", name)
gcvec <- rep("GC-FIBRBLS", length(name))
gcvec[grep("T", name)] <- "GC-TCELL"
gcvec[grep("L", name)] <- "GC-LCL"
gcvec[grep("B", name)] <- "GC-LCL"
gcvec <- gcvec[log==1]

#proportion of biallelic sites per tissue as barplot
spl <- split(gtp, tvecf)
spl2 <- split(gcp, gcvec)
totlist <- c(spl, spl2, list(gdp))
names(totlist)[length(totlist)] <- "GD-LCL"
lim <- c(rep(0.9, length(spl)), rep(0.8, length(spl2)), 0.7)
cl <- vector("numeric", length(totlist))
for(i in 1:length(cl)) {
	cl[i] <- length(which(totlist[[i]]<lim[i]))
}

br <- seq(0, 1, 0.05)
X11(height=3, width=9)
#pdf("rplots/chrx_biall_histograms.pdf", height=3, width=9)
par(mfrow=c(1,3))
hist(gtp, breaks=br, main="GTEx", col="grey", xlab="chrX biallelic/total", ylab="Samples", cex.axis=1.5, cex.lab=1.5)
abline(v=0.9, lwd=2, col="royalblue2")
hist(gcp, breaks=br, main="Gencord", col="grey", xlab="chrX biallelic/total", ylab="Samples", cex.axis=1.5, cex.lab=1.5)
abline(v=0.8, lwd=2, col="royalblue2")
hist(gdp, breaks=br, main="Geuvadis", col="grey", xlab="chrX biallelic/total", ylab="Samples", cex.axis=1.5, cex.lab=1.5)
abline(v=0.7, lwd=2, col="royalblue2")
dev.off()

X11(height=5, width=7)
#pdf("rplots/chrx_monoclonal_barplot.pdf", height=5, width=7)
par(mar=c(7,5,5,3))
barplot(rbind(unlist(lapply(totlist, length))-cl, cl), las=3, ylab="Samples", col=c("royalblue2", "cyan3"), cex.axis=1.2, cex.lab=1.2)
legend("topleft", fill=c("royalblue2", "cyan3"), legend=c("Polyclonal", "Monoclonal"), bty="n", cex=1.2)
dev.off()

prop <-  cl/unlist(lapply(totlist, length))
tab <- cbind(cl, unlist(lapply(totlist, length)), prop)
colnames(tab) <- c("PUTATIVE_CLONAL_SAMPLES_FEMALE", "TOTAL_SAMPLES_FEMALE", "PROP")
write.table(tab, "putative_clonal_samples_stats.txt", sep="\t", quote=F)





###########################################################################################################################################
#                                                   Novelty and comparison to the mouse                                                             #
###########################################################################################################################################


ot <- read.table("otago.list.txt", as.is=T)
o <- unique(ot[,c(2,5,6)])
h <- read.table("biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
hn <- mat.or.vec(nrow(h),2)
for(i in 1:nrow(h)) {
	tmp <- o[o[,1]==rownames(h)[i],]
	if (nrow(tmp)==1) {
		hn[i,1] <- tmp[1,2]
		hn[i,2] <- tmp[1,3]
	}
	if (nrow(tmp)==0) {
		hn[i,1] <- NA
		hn[i,2] <- NA
	}
}
rownames(hn) <- rownames(h)

table(as.data.frame(hn))
#   - GI Prov
#  B 18  0    0
#  H  8  2    2

status <- hn[,1]
status[is.na(status)] <- "Novel"
status[hn[,2]=="GI"] <- "Provisional"
status[hn[,2]=="Prov"] <- "Provisional"
status[hn[,2]=="-"] <- "Known"

write.table(status, "final_knownstatus.txt", col.names=F, sep="\t", quote=F)

known_loci <- c("MEG9", "SNHG14", "LPAR6")
provisional <- rownames(hn)[!is.na(hn[,2]) & (hn[,2]=="GI" | hn[,2]=="Prov")]
human_mouse <- rownames(hn)[!is.na(hn[,1]) & hn[,1]=="B"]
human <- setdiff(rownames(hn)[!is.na(hn[,1]) & hn[,1]=="H"], provisional)
novel <- setdiff(rownames(hn)[is.na(hn[,1])], known_loci)
l <- list(human_mouse, human, provisional, known_loci, novel)
lapply(l, function(x) {paste(x, collapse=", ")})

#[1] "PEG10, SNRPN, ZDBF2, KCNQ1, MEG3, MEST, H19, IGF2, DLK1, PLAGL1, GRB10, SNURF, NAP1L5, PEG3, NDN, UBE3A, MEG8, MAGEL2"
#[1] "FAM50B, CPA4, ZNF597, INPP5F_V2, L3MBTL1, PWRN1, IGF2-AS, SGK2"
#[1] "MAGI2, ZNF331, DIRAS3, NTM"
#[1] "MEG9, SNHG14, LPAR6"
#[1] "THEGL, PPIEL, UTS2, SYCE1, RP11-7F17.7, PRSS50, KIF25, UGT2B4, CST1"


vec <- c(unlist(lapply(l, length)))
names(vec) <- c("HS+MM", "HSconf", "HSprov", "NovelKnownLoci", "Novel")
X11()
pdf("rplots/genes_classes_barplot.pdf")
par(mar=c(10,4,3,2))
barplot(vec, las=3, col=c("springgreen3", "blue", "cyan3", "purple", "magenta3"), cex.names=1.3, ylab="Genes", cex.axis=1.3, cex.lab=1.3)
dev.off()



ot <- read.table("otago.list.txt", as.is=T)
known_bitab2 <- read.table("biallelic.weak.known.tab.txt", as.is=T, row.names=1, header=T)
known_bitab <- read.table("biallelic.strong.known.tab.txt", as.is=T, row.names=1, header=T)
final_bitab2 <- read.table("biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
final_bitab <- read.table("biallelic.strong.final.tab.txt", as.is=T, row.names=1, header=T)
known_bitab2 <- unique(rbind(known_bitab2, final_bitab2))
known_bitab <- unique(rbind(known_bitab, final_bitab))

#prepare Table S6
Imprinting_status <- vector("character", nrow(ot))
for(i in 1:nrow(ot)) {
	if (ot[i,2]=="-") {
		Imprinting_status[i] <- "not evaluated"
	} 
	if (ot[i,5]=="M") {
		Imprinting_status[i] <- "not evaluated"
	} 
	if (ot[i,2]!="-" & ot[i,5]!="M") {
		if (any(rownames(known_bitab2)!=ot[i,2])) {
			Imprinting_status[i] <- "no data"
		} 
		if (any(rownames(known_bitab2)==ot[i,2])) {
			vec <- na.omit(as.character(known_bitab[ot[i, 2], ]))
			vec2 <- na.omit(as.character(known_bitab2[ot[i, 2], ]))
			Imprinting_status[i] <- "inconclusive data"
			if (any(vec2=="+")) { 
				Imprinting_status[i] <- "consistent with biallelic"
			}
			if (any(vec=="+")) { 
				Imprinting_status[i] <- "biallelic"
			}
			if (any(vec2=="-")) { 
				Imprinting_status[i] <- "consistent with imprinted"
			}
			if (any(vec=="-")) {
				Imprinting_status[i] <- "imprinted"
			}
		} 
	}
}
Imprinting_status[grep("MIR", ot[,2])] <- "small RNA gene"
Imprinting_status[grep("SNOR", ot[,2])] <- "small RNA gene"

table(Imprinting_status)
#                biallelic consistent with biallelic consistent with imprinted 
#                       41                         4                        28 
#                imprinted                   no data             not evaluated 
#                       32                        19                        36 
#           small RNA gene 
#                       37 
                       
colnames(ot) <- c("Human_GeneID", "Human_Genename", "Mouse_GeneID", "Mouse_Genename", "Imprinted_species", "Human_notes", "Mouse_notes")
ot[,5][ot[,5]=="H"] <- "human"
ot[,5][ot[,5]=="B"] <- "human, mouse"
ot[,5][ot[,5]=="M"] <- "mouse"
ot[,6][ot[,6]=="pdf"] <- "from otago pdf"
ot[,6][ot[,6]=="GI"] <- "from geneimprint"
ot[,6][ot[,6]=="Prov"] <- "provisional"
ot[,6][ot[,6]=="Conf"] <- "conflicting"
ot[,7][ot[,7]=="Ques"] <- "questionable"
ot[,7][ot[,7]=="Prov"] <- "provisional"
ot[,7][ot[,7]=="X"] <- "chrX"

Imprinting_status_Baran <- Imprinting_status

ot2 <- cbind(ot, Imprinting_status_Baran)
ot2[,8] <- replace(as.character(ot2[,8]), as.character(ot2[,2])=="NLRP2", "failed validation")
ot2[,8] <- replace(as.character(ot2[,8]), as.character(ot2[,2])=="HM13", "failed validation")

ot2f <- ot2[ot2[,5]!="mouse",]
write.table(ot2f, "table_S6_no_mouse.txt", row.names=F, sep="\t", quote=F)

ot2f <- ot2f[ot2f[,6]=="-" | ot2f[,6]=="from otago pdf",]
write.table(ot2f, "table_S6_no_mouse_otago_only.txt", row.names=F, sep="\t", quote=F)

length(unique(ot2f[,2]))
#107

length(unique(ot2f[,2][ot2f[,8]=="no data"]))
#16

length(unique(ot2f[,2][ot2f[,8]=="small RNA gene"]))
#miRNA/snoRNA 12

ot2f <- ot2f[ot2f[,8]!="small RNA gene",]
write.table(ot2f, "table_S6_no_mouse_otago_only_no_smallRNA.txt", row.names=F, sep="\t", quote=F)

length(unique(ot2f[,2]))
#95

length(unique(ot2f[,2][ot2f[,8]=="biallelic" | ot2f[,8]=="consistent with biallelic" | ot2f[,8]=="imprinted" | ot2f[,8]=="consistent with imprinted"]))
#79

length(unique(ot2f[,2][ot2f[,8]=="imprinted"]))
#26

length(unique(ot2f[,2][ot2f[,8]=="biallelic"]))
#31 biallelic, no signs of imprinting



###########################################################################################################################################
#                                                              Gene scatterplots                                                          #
###########################################################################################################################################




##########  SNP-based ref/alt counts  ##########


cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint

datal <- vector("list", 3)
datal[[1]] <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.final_impr.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.final_impr.txt", as.is=T, header=T, sep="\t", fill=T)
colnames(datal[[2]]) <- colnames(datal[[1]])
datal[[2]] <- datal[[2]][,1:28]
datal[[3]] <- read.table("gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.final_impr.txt", as.is=T, header=T, sep="\t")
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,1], datal[[i]][,2], sep="-")
}
datal[[2]][,"TISSUE"] <- "GD-LCL"
datal[[3]][,"TISSUE"][datal[[3]][,"TISSUE"]=="LCL"] <- "GC-LCL"
datal[[3]][,"TISSUE"][datal[[3]][,"TISSUE"]=="FIBRBLS"] <- "GC-FIBRBLS"

data <- do.call(rbind, datal)
spl <- split(data, data[,22])
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
g <- gen[,2]
ge <- gen[,1]
tis <- names(table(data$TISSUE))
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
spl2 <- spl2[ge]
spl <- spl[ge]
for(i in 1:length(ge)) {
	X11(height=5,width=5)
#	pdf(paste("rplots_snp_scatter/tissuescatters_", g[i], ".pdf", sep=""))
	for(j in 1:length(spl2[[i]]))  {
		ma <- max(c(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9]))
		plot(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9], main=paste(g[i], names(spl2[[i]])[j]), xlab="REF", ylab="NONREF", xlim=c(0,ma), ylim=c(0,ma))
	}
	dev.off()
}


#other genes of interest

datal <- vector("list", 3)
datal[[1]] <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.extra.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.extra.txt", as.is=T, header=T, sep="\t", fill=T)
colnames(datal[[2]]) <- colnames(datal[[1]])
datal[[2]] <- datal[[2]][,1:28]
datal[[3]] <- read.table("gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.extra.txt", as.is=T, header=T, sep="\t")
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,1], datal[[i]][,2], sep="-")
}
datal[[2]][,"TISSUE"] <- "GD-LCL"
datal[[3]][,"TISSUE"][datal[[3]][,"TISSUE"]=="LCL"] <- "GC-LCL"
datal[[3]][,"TISSUE"][datal[[3]][,"TISSUE"]=="FIBRBLS"] <- "GC-FIBRBLS"

data <- do.call(rbind, datal)
spl <- split(data, data[,22])
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gen <- read.table("extra_ensg2name.txt", as.is=T)
g <- gen[,2]
ge <- gen[,1]
tis <- names(table(data$TISSUE))
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
spl2 <- spl2[ge]
spl <- spl[ge]
for(i in 1:length(ge)) {
#	X11(height=5,width=5)
	pdf(paste("rplots_snp_scatter_extra/tissuescatters_", g[i], ".pdf", sep=""))
	for(j in 1:length(spl2[[i]]))  {
		ma <- max(c(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9]))
		plot(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9], main=paste(g[i], names(spl2[[i]])[j]), xlab="REF", ylab="NONREF", xlim=c(0,ma), ylim=c(0,ma))
	}
	dev.off()
}


g <- c("PAX8")
t1 <- "LUNG"
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
#X11(height=5,width=5)
pdf(paste("rplots/snpscatter_PAX8.pdf", sep=""), height=5,width=5)
j=1
i=1
td1 <- spl[[ge[i]]]
td1 <- td1[td1[,28]=="LUNG",]
plot(td1[,8], td1[,9], xlab=paste("REF count"), ylab=paste("ALT count"), col="black", lwd=2, cex.axis=1.5, cex.lab=1.5)
legend("topleft", legend=c("PAX8 LUNG"), cex=1.5, bty="n")
dev.off()




#known genes

datal <- vector("list", 3)
datal[[1]] <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.otago.human.bi.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.otago.human.bi.txt", as.is=T, header=T, sep="\t", fill=T)
colnames(datal[[2]]) <- colnames(datal[[1]])
datal[[2]] <- datal[[2]][,1:28]
datal[[3]] <- read.table("gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.otago.human.bi.txt", as.is=T, header=T, sep="\t")
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,1], datal[[i]][,2], sep="-")
}
datal[[2]][,"TISSUE"] <- "GD-LCL"
datal[[3]][,"TISSUE"][datal[[3]][,"TISSUE"]=="LCL"] <- "GC-LCL"
datal[[3]][,"TISSUE"][datal[[3]][,"TISSUE"]=="FIBRBLS"] <- "GC-FIBRBLS"

data <- do.call(rbind, datal)
spl <- split(data, data[,22])
#length(spl)
#[1] 77
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gen <- read.table("otago.human.bi_ensg2name.txt", as.is=T)
g <- gen[,2]
ge <- gen[,1]
tis <- names(table(data$TISSUE))
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
spl2 <- spl2[ge]
spl <- spl[ge]
#length(spl)
#[1] 77
for(i in 1:length(ge)) {
#	X11(height=5,width=5)
	pdf(paste("rplots_snp_scatter_known/tissuescatters_", g[i], ".pdf", sep=""))
	for(j in 1:length(spl2[[i]]))  {
		ma <- max(c(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9]))
		plot(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9], main=paste(g[i], names(spl2[[i]])[j]), xlab="REF", ylab="NONREF", xlim=c(0,ma), ylim=c(0,ma))
	}
	dev.off()
}




datal <- vector("list", 1)
datal[[2]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.failed_examples.txt", as.is=T, header=T, sep="\t", fill=T)
colnames(datal[[2]]) <- colnames(datal[[1]])
datal[[2]] <- datal[[2]][,1:28]
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,1], datal[[i]][,2], sep="-")
}
datal[[2]][,"TISSUE"] <- "GD-LCL"
data <- do.call(rbind, datal)
spl <- split(data, data[,22])
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)

g <- c("MAP2K3")
t1 <- "GD-LCL"
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
#X11(height=5,width=5)
pdf(paste("rplots/snpscatter_MAP2K3.pdf", sep=""), height=5,width=5)
j=1
i=1
td1 <- spl[[ge[i]]]
plot(td1[,8], td1[,9], xlab=paste("REF count"), ylab=paste("ALT count"), col="black", lwd=2, cex.axis=1.5, cex.lab=1.5)
legend("topright", legend=c("MAP2K3 GD-LCL"), cex=1.5, bty="n")
dev.off()

g <- c("UQCRFS1")
t1 <- "GD-LCL"
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
#X11(height=5,width=5)
pdf(paste("rplots/snpscatter_UQCRFS1.pdf", sep=""), height=5,width=5)
j=1
i=1
td1 <- spl[[ge[i]]]
plot(td1[,8], td1[,9], xlab=paste("REF count"), ylab=paste("ALT count"), col="black", lwd=2, cex.axis=1.5, cex.lab=1.5)
legend("topleft", legend=c("UQCRFS1 GD-LCL"), cex=1.5, bty="n")
dev.off()






##brain subtissues

br_subs <- read.table("../../data/gtex/braintissues.txt", as.is=T)[,1]
datal <- vector("list", length(br_subs))
for(i in 1:length(br_subs)) {
	datal[[i]] <- read.table(paste("~/tuuli_lab/tlappalainen/imprint/data/gtex/brsub/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN_", br_subs[i], ".impr.txt", sep=""), as.is=T, header=F, sep="\t")
	rownames(datal[[i]]) <- paste(datal[[i]][,1], datal[[i]][,2], sep="-")
}
data <- do.call(rbind, datal)
spl <- split(data, data[,22])
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("~/tuuli_lab/tlappalainen/imprint/data/imprint/gtex/ensg2name.txt", as.is=T, sep=",")
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
g <- gen[,2]
ge <- gen[,1]
tis <- names(table(data$TISSUE))
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
spl2 <- spl2[ge]
spl <- spl[ge]

for(i in 1:length(ge)) {
#	X11(height=5,width=5)
	pdf(paste("rplots_snp_scatter/tissuescatters_brsub_", g[i], "_pub.pdf", sep=""))
	for(j in 1:length(spl2[[i]]))  {
		ma <- max(c(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9]))
		plot(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9], main=paste(g[i], names(spl2[[i]])[j]), xlab="REF", ylab="NONREF", xlim=c(0,ma), ylim=c(0,ma))
	}
	dev.off()
}


g <- c("L3MBTL1")
ge <- gen[gen[,2]==g,1]
i <- which(names(spl2)==ge)
indx <- match(c("BRNCHB", "BRNPTM", "BRNACC", "BRNCTXA"), names(spl2[[i]]), nomatch=F)
pdf(paste("rplots/tissuescatters_brsub_", g, ".pdf", sep=""))
#X11()
par(mfrow=c(2,2))
for(j in indx)  {
ma <- max(c(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9]))
	plot(spl2[[i]][[j]][,8], spl2[[i]][[j]][,9], main=paste(names(spl2[[i]])[j]), xlab="REF count", ylab="NONREF count", xlim=c(0,ma), ylim=c(0,ma), cex.axis=1.3, cex.lab=1.3, lwd=2)
}
dev.off()






##########  Haplotype-based counts  ##########


datal <- vector("list", 3)
datal[[1]] <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.final_impr.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.final_impr.txt", as.is=T, header=T, sep="\t")
datal[[3]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.filter5.model15.indxgen.final_impr.txt", as.is=T, header=T, sep="\t", fill=T)
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,2], datal[[i]][,3], datal[[i]][,1], sep="-")
}
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="LCL"] <- "GC-LCL"
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="FIBRBLS"] <- "GC-FIBRBLS"
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="TCELL"] <- "GC-TCELL"
datal[[3]][,"tissue"] <- "GD-LCL"
data <- do.call(rbind, datal)
spl <- split(data, data[,2])
spl2 <- lapply(spl, function(x) {split(x, x[,1])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,3], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
g <- gen[,2]
ge <- gen[,1]
tis <- names(table(data$TISSUE))
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
spl2 <- spl2[ge]
spl <- spl[ge]
for(i in 1:length(ge)) {
#	X11(height=5,width=5)
	pdf(paste("rplots_haplo_scatter/tissuescatters_haplo_imp_", g[i], ".pdf", sep=""))
	for(j in 1:length(spl2[[i]]))  {
		ma <- max(c(spl2[[i]][[j]][,4], spl2[[i]][[j]][,5]))
		plot(spl2[[i]][[j]][,4], spl2[[i]][[j]][,5], main=paste(g[i], names(spl2[[i]])[j]), xlab="REF", ylab="NONREF", xlim=c(0,ma), ylim=c(0,ma))
	}
	dev.off()
}


data <- datal[[1]]
spl <- split(data, data[,2])
spl2 <- lapply(spl, function(x) {split(x, x[,1])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,3], sep="-")
	}
}
g <- c("MEST")
t1 <- "LUNG"
t2 <- "TESTIS"
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]

#X11(height=5,width=5)
pdf(paste("rplots/haploscatter_MEST.pdf", sep=""), height=5,width=5)
j=1
i=1
td1 <- spl2[[ge[i]]][[t1]]
plot(td1[,4], td1[,5], xlab=paste("Haplo 1 count"), ylab=paste("Haplo 2 count"), col="royalblue3", lwd=2, cex.axis=1.5, cex.lab=1.5, xlim=c(0,100), ylim=c(0,100))
td2 <- spl2[[ge[i]]][[t2[j]]]
points(td2[,4], td2[,5], col="black", lwd=2)
legend("topleft", fill=c("royalblue3", "black"), legend=c("MEST LUNG", "MEST TESTIS"), bty="n")
dev.off()
#excluding 7 data points >100 counts

data <- datal[[1]]
spl <- split(data, data[,2])
spl2 <- lapply(spl, function(x) {split(x, x[,1])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,3], sep="-")
	}
}
g <- c("DLK1")
t1 <- "HRTLV"
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]

#X11(height=5,width=5)
pdf(paste("rplots/haploscatter_DLK1.pdf", sep=""), height=5,width=5)
i=1
td1 <- spl2[[ge[i]]][[t1]]
plot(td1[,4], td1[,5], xlab=paste("Haplo 1 count"), ylab=paste("Haplo 2 count"), col="royalblue3", lwd=2, cex.axis=1.5, cex.lab=1.5, xlim=c(0,150), ylim=c(0,150))
legend("topleft", text.col=c("royalblue3"), legend=c("DLK1 HRTLV"), bty="n")
dev.off()







##########  SNP-based ref-ratio tissue direction  ##########


data <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.final_impr.txt", as.is=T, header=T)
rownames(data) <- paste(data[,1], data[,2], sep="-")
spl <- split(data, data[,22])
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]

g <- c("ZDBF2", "GRB10", "IGF2")
t1 <- "BRAIN"
t2 <- "MSCLSK"
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
col <- c("blue", "magenta", "red2")
lmlist <- vector("list", 3)
rhotab <- mat.or.vec(length(g), 2)
pdf(paste("rplots/tissueplots_alltis_IGF2_GRB10_ZDBF2.pdf", sep=""), height=5,width=5)
j=1
for (i in 1:3) {
	int <- intersect(rownames(spl2[[ge[i]]][[t1]]), rownames(spl2[[ge[i]]][[t2[j]]]))
	td1 <- spl2[[ge[i]]][[t1]][int,]
	td2 <- spl2[[ge[i]]][[t2[j]]][int,]
	r1 <- (td1[,8]/td1[,10])
	r2 <- (td2[,8]/td2[,10])
	lmlist[[i]] <- lm(r2 ~ r1)
	rhotab[i,1] <- cor.test(r1, r2, method="s")$estimate
	rhotab[i,2] <- cor.test(r1, r2, method="s")$p.value
	if (i==1) {
		plot(td1[,8]/td1[,10], td2[,8]/td2[,10], xlab=paste("Ref/Total Brain"), ylab=paste("Ref/Total Muscle"), xlim=c(0,1), ylim=c(0,1), col=col[i], lwd=2, cex.axis=1.5, cex.lab=1.5)
	} 
	if (i>1) {
		points(td1[,8]/td1[,10], td2[,8]/td2[,10], main=g[i], col=col[i], lwd=2)
	} 
	abline(lmlist[[i]]$coefficients, col=col[i])	
}
legend(x=0.5, y=1, fill=c("blue", "magenta3", "red2"), legend=c("ZDBF2", "GRB10", "IGF2"), cex=1.1, bg="white")
dev.off()

max(rhotab[,2])
#0.0004105211 , correlation p-value




#####  IGF2 & GRB10 brain vs all scatters  ######


data <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.final_impr.txt", as.is=T, header=T)
rownames(data) <- paste(data[,1], data[,2], sep="-")
spl <- split(data, data[,22])
spl2 <- lapply(spl, function(x) {split(x, x[,28])})
for(i in 1:length(spl2)) {
	for(j in 1:length(spl2[[i]])) {
		rownames(spl2[[i]][[j]]) <- paste(spl2[[i]][[j]][,2], spl2[[i]][[j]][,27], sep="-")
	}
}
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]

g <- c("GRB10", "IGF2")
t1 <- "BRAIN"
tis <- unique(data$TISSUE)
t2 <- tis[tis!=t1]
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]

for(i in 1:length(ge)) {
	pdf(paste("rplots/tissueplots_alltis_", g[i], ".pdf", sep=""), height=15,width=7.5)
	par(mfrow=c(8,4))
	for(j in 1:length(t2))	{
		if(any(names(spl2[[ge[i]]])==t2[j])) {	
			int <- intersect(rownames(spl2[[ge[i]]][[t1]]), rownames(spl2[[ge[i]]][[t2[j]]]))
			if(length(int)>0) {	
			td1 <- spl2[[ge[i]]][[t1]][int,]
			td2 <- spl2[[ge[i]]][[t2[j]]][int,]
			plot(td1[,8]/td1[,10], td2[,8]/td2[,10], xlab=paste("Ref ratio", t1), ylab=paste("Ref ratio", t2[j]), xlim=c(0,1), ylim=c(0,1))
			abline(h=0.5, col="blue")
			abline(v=0.5, col="blue")
			}
		}		
	}
	dev.off()
}








###########################################################################################################################################
#                                                   Imprinting across tissues                                                             #
###########################################################################################################################################


#########   Heatmaps    ###########

cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint

datal <- vector("list", 3)
datal[[1]] <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.final_impr.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", as.is=T, header=T, sep="\t")
datal[[3]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", as.is=T, header=T, sep="\t", fill=T)
known <- read.table("known_otago_ensg2name.txt", as.is=T)[,1]
status <- read.table("final_knownstatus.txt", row.names=1, sep="\t", as.is=T)
putative <- rownames(status)[status=="Provisional"]
k <- rownames(status)[status=="Known"]
ensg <- read.table("ensg2nametab.txt", as.is=T, sep="\t", row.names=1)
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="LCL"] <- "GC-LCL"
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="FIBRBLS"] <- "GC-FIBRBLS"
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="TCELL"] <- "GC-TCELL"
datal[[3]][,"tissue"] <- "GD-LCL"
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,2], datal[[i]][,1], sep="-")
}
data <- do.call(rbind, datal)
tis <- unique(unlist(lapply(datal, function(x) {x[,1]})))
g <- unique(unlist(lapply(datal, function(x) {x[,3]})))
gn <- g
for(i in 1:length(g)) {
	if (any(putative==g[i])) {
		gn[i] <- paste(g[i], "*", sep="")
	}
	if (any(k==g[i])) {
		gn[i] <- paste(g[i], "**", sep="")
	}
}

stat <- c("max_m") ## == tau
tabs <- vector("list", length(stat))
names(tabs) <- stat
for(i in 1:length(stat)) {
	tabs[[i]] <- mat.or.vec(length(g), length(tis))
	rownames(tabs[[i]]) <- g
	colnames(tabs[[i]]) <- tis
	for(x in 1:length(g)) {
		for(y in 1:length(tis)) {
			if (length(data[data[,3]==g[x] & data[,1]==tis[y],stat[i]])==0) {
				tabs[[i]][x,y] <- NA 
			} else {
				tabs[[i]][x,y] <- data[data[,3]==g[x] & data[,1]==tis[y],stat[i]]
			}
		}
	}
}
printtab <- tabs[["max_m"]]
rownames(printtab) <- g
write.table(printtab, "impr_tau_tab.txt", quote=F, sep="\t")

for(i in 1:length(stat)) {
	rownames(tabs[[i]]) <- gn
}



pos <- read.table("~/tuuli_lab/tlappalainen/resource/gencode.v12.annotation.tab.gene", as.is=T)
rownames(pos) <- unlist(lapply(strsplit(pos[,9], ".", fixed=T), function(x) {x[1]}))
ensgi <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
ensgi <- ensgi[!duplicated(ensgi[,2]),]
row.names(ensgi) <- ensgi[,2]
poss <- pos[ensgi[g,1],]
tss <- poss[,4]
tss[poss[,7]=="-"] <- poss[,5][poss[,7]=="-"]
p <- paste(poss[,1], round(tss/1000000, 1), sep=":")
c <- as.numeric(substr(poss[,1], 4, nchar(poss[,1])))
p2 <- paste(c, round(tss/1000000, 1), sep=":")
s <- order(c, tss, decreasing=T)


#save positions and data for landscape plots

pos_filt <- pos[pos[,14]=="lincRNA" | pos[,14]=="protein_coding",]
gene_bed <- pos_filt[,c(1,4,5,13)]
col <- rep("255,255,255", nrow(gene_bed))
gene_bed <- cbind(gene_bed, rep("1000", nrow(gene_bed)), pos_filt[,7], pos_filt[,4], pos_filt[,5], col)
write.table(gene_bed, "gene_bed", sep="\t", col.names=F, row.names=F, quote=F)
p_print <- paste(poss[,1], paste(tss-1500000, tss+1500000, sep="-"), sep=":")
write.table(sort(p_print), "positions_1.5MB.txt", quote=F, row.names=F, col.names=F)

i=1
#	X11(height=9, width=8)
pdf(paste("rplots/heatmap_impr_tau.pdf", sep=""), height=9, width=8)
par(mar=c(7,8,3,6))
ptab <- tabs[[i]][s,]
natab <- tabs[[i]][s,]
natab[is.na(tabs[[i]][s,])] <- 1
natab[!is.na(tabs[[i]][s,])] <- NA
image((max(ptab, na.rm=T))-(t(ptab)), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, col=heat.colors(20))
image(t(natab), add=T, xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, col="grey")
axis(at=0:(nrow(t(ptab))-1)/(nrow(t(ptab))-1), labels=rownames(t(ptab)), side=1, las=3)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=colnames(t(ptab)), side=2, las=1)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=p[s], side=4, las=1)		
dev.off()

pdf(paste("rplots/heatmap_impr_tau_legend.pdf", sep=""))
plot(0,0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
legend("topleft", fill=c(heat.colors(20), "grey"), legend=rep("", 21), bty="n", x.intersp=-1, cex=2, ncol=21)
text(x=seq(-0.925, 0.875, length.out=21), y=rep(0.75, 21), las=3, c(round(seq(max(ptab, na.rm=T), min(ptab, na.rm=T), length.out=20), 2), "NA"), srt=90, cex=1.5)
dev.off()




#########   Heatmap of brain subtissues    ###########


###brain subtissues
cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint

datal <- vector("list", 1)
datal[[1]] <- read.table("~/tuuli_lab/tlappalainen/imprint/data/imprint/gtex/brainsub_res_proc/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BRAIN.txt.impglr", as.is=T, header=F, sep="\t")
known <- read.table("known_otago_ensg2name.txt", as.is=T)[,1]
status <- read.table("final_knownstatus.txt", row.names=1, sep="\t", as.is=T)
putative <- rownames(status)[status=="Provisional"]
k <- rownames(status)[status=="Known"]
ensg <- read.table("ensg2nametab.txt", as.is=T, sep="\t", row.names=1)
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,2], datal[[i]][,1], sep="-")
}
data <- do.call(rbind, datal)
tis <- unique(unlist(lapply(datal, function(x) {x[,1]})))
g <- unique(unlist(lapply(datal, function(x) {x[,3]})))
gn <- g
for(i in 1:length(g)) {
	if (any(putative==g[i])) {
		gn[i] <- paste(g[i], "*", sep="")
	}
	if (any(k==g[i])) {
		gn[i] <- paste(g[i], "**", sep="")
	}
}
gn[gn=="INPP5F"] <- "INPP5F_V2**"
stat <- c("max_m")
tabs <- vector("list", length(stat))
names(tabs) <- stat
for(i in 1:length(stat)) {
	tabs[[i]] <- mat.or.vec(length(g), length(tis))
	rownames(tabs[[i]]) <- g
	colnames(tabs[[i]]) <- tis
	for(x in 1:length(g)) {
		for(y in 1:length(tis)) {
			if (length(data[data[,3]==g[x] & data[,1]==tis[y],32])==0) {
				tabs[[i]][x,y] <- NA 
			} else {
				tabs[[i]][x,y] <- data[data[,3]==g[x] & data[,1]==tis[y],32]
			}
		}
	}
}
printtab <- tabs[["max_m"]]
rownames(printtab) <- g
write.table(printtab, "impr_tau_brainsub_tab.txt", quote=F, sep="\t")
for(i in 1:length(stat)) {
	rownames(tabs[[i]]) <- gn
}

pos <- read.table("~/tuuli_lab/tlappalainen/resource/gencode.v12.annotation.tab.gene", as.is=T)
rownames(pos) <- unlist(lapply(strsplit(pos[,9], ".", fixed=T), function(x) {x[1]}))
ensgi <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
ensgi <- ensgi[!duplicated(ensgi[,2]),]
row.names(ensgi) <- ensgi[,2]
poss <- pos[ensgi[g,1],]
tss <- poss[,4]
tss[poss[,7]=="-"] <- poss[,5][poss[,7]=="-"]
p <- paste(poss[,1], round(tss/1000000, 1), sep=":")
c <- as.numeric(substr(poss[,1], 4, nchar(poss[,1])))
p2 <- paste(c, round(tss/1000000, 1), sep=":")
s <- order(c, tss, decreasing=T)
i=1
pdf(paste("rplots/heatmap_brainsub_tau.pdf", sep=""), height=9, width=5)
par(mar=c(7,8,3,6))
ptab <- tabs[[i]][s,]
natab <- tabs[[i]][s,]
natab[is.na(tabs[[i]][s,])] <- 1
natab[!is.na(tabs[[i]][s,])] <- NA
image(max(ptab, na.rm=T)-t(ptab), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, col=heat.colors(20))
image(t(natab), add=T, xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, main=paste("log2", stat[i]), col="grey")
axis(at=0:(nrow(t(ptab))-1)/(nrow(t(ptab))-1), labels=rownames(t(ptab)), side=1, las=3)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=colnames(t(ptab)), side=2, las=1)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=p[s], side=4, las=1)		
dev.off()




#########   Heatmap of known genes    ###########



datal <- vector("list", 3)
datal[[1]] <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.known.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.known.txt", as.is=T, header=T, sep="\t")
datal[[3]] <- read.table("geuvadis/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.impglr.known.txt", as.is=T, header=T, sep="\t", fill=T)
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="LCL"] <- "GC-LCL"
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="FIBRBLS"] <- "GC-FIBRBLS"
datal[[2]][,"tissue"][datal[[2]][,"tissue"]=="TCELL"] <- "GC-TCELL"
datal[[3]][,"tissue"] <- "GD-LCL"
for(i in 1:length(datal)) {
	rownames(datal[[i]]) <- paste(datal[[i]][,2], datal[[i]][,1], sep="-")
}
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
data <- do.call(rbind, datal)
tis <- unique(unlist(lapply(datal, function(x) {x[,1]})))
g <- unique(unlist(lapply(datal, function(x) {x[,3]})))
stat <- c("max_m")
tabs <- vector("list", length(stat))
names(tabs) <- stat
for(i in 1:length(stat)) {
	tabs[[i]] <- mat.or.vec(length(g), length(tis))
	rownames(tabs[[i]]) <- g
	colnames(tabs[[i]]) <- tis
	for(x in 1:length(g)) {
		for(y in 1:length(tis)) {
			if (length(data[data[,3]==g[x] & data[,1]==tis[y],stat[i]])==0) {
				tabs[[i]][x,y] <- NA 
			} else {
				tabs[[i]][x,y] <- median(data[data[,3]==g[x] & data[,1]==tis[y],stat[i]])
			}
		}
	}
}
##remove imprinted genes
for(i in 1:length(stat)) {
	sd <- setdiff(rownames(tabs[[i]]), gen[,2])
	tabs[[i]] <- tabs[[i]][sd,]
	g <- rownames(tabs[[i]])
	rownames(tabs[[i]]) <- paste(rownames(tabs[[i]]), "**", sep="")
}
pos <- read.table("~/tuuli_lab/tlappalainen/resource/gencode.v12.annotation.tab.gene", as.is=T)
rownames(pos) <- unlist(lapply(strsplit(pos[,9], ".", fixed=T), function(x) {x[1]}))
ensgi <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
ensgi <- ensgi[!duplicated(ensgi[,2]),]
row.names(ensgi) <- ensgi[,2]
poss <- pos[ensgi[g,1],]
tss <- poss[,4]
tss[poss[,7]=="-"] <- poss[,5][poss[,7]=="-"]
p <- paste(poss[,1], round(tss/1000000, 1), sep=":")
c <- as.numeric(substr(poss[,1], 4, nchar(poss[,1])))
p2 <- paste(c, round(tss/1000000, 1), sep=":")
s <- order(c, tss, decreasing=T)
i=1
pdf(paste("rplots/heatmap_known_tau.pdf", sep=""), height=12, width=8)
par(mar=c(7,8,3,6))
ptab <- tabs[[i]][s,]
natab <- tabs[[i]][s,]
natab[is.na(tabs[[i]][s,])] <- 1
natab[!is.na(tabs[[i]][s,])] <- NA
image((max(ptab, na.rm=T))-(t(ptab)), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, col=heat.colors(20))
image(t(natab), add=T, xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, col="grey")
axis(at=0:(nrow(t(ptab))-1)/(nrow(t(ptab))-1), labels=rownames(t(ptab)), side=1, las=3)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=colnames(t(ptab)), side=2, las=1)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=p[s], side=4, las=1)		
dev.off()








#########   Barplots    ###########


bitab <- read.table("biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
bitab2 <- read.table("biallelic.strong.final.tab.txt", as.is=T, row.names=1, header=T)
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
bitab <- bitab[intersect(gen[,2], rownames(bitab)),]
bitab2 <- bitab2[rownames(bitab),]

counts <- mat.or.vec(ncol(bitab), 6)
colnames(counts) <- c("BISTR", "BIW", "MOSTR", "MOW", "QU", "XNA") 
rownames(counts) <- colnames(bitab)
for(i in 1:ncol(bitab)) {
	counts[i,1] <- length(which(bitab2[,i]=="+" & !is.na(bitab2[,i]) & bitab[,i]!="?"))
	counts[i,2] <- length(which(bitab2[,i]=="?" & bitab[,i]=="+" & !is.na(bitab2[,i])  & !is.na(bitab[,i])))
	counts[i,3] <- length(which(bitab2[,i]=="-" & !is.na(bitab2[,i]) & bitab[,i]!="?"))
	counts[i,4] <- length(which(bitab2[,i]=="?" & bitab[,i]=="-" & !is.na(bitab2[,i]) & !is.na(bitab[,i])))
	counts[i,5] <- length(which(bitab[,i]=="?" & !is.na(bitab[,i])))
	counts[i,6] <- length(which(!is.na(bitab[,i])))
}

s3 <- order((counts[,3] + counts[,4])/counts[,6], decreasing=T)
#X11(width=7, height=5)
pdf("rplots/biallelic.impr.barplot.pdf", width=10, height=7)
par(mar=c(7,5,3,2))
barplot(t(apply(counts[,c(1,2,4,3,5)], 2, function(x) {x/counts[,6]}))[,s3], las=3, col=c("yellow", "goldenrod2", "red", "red3", "grey"), ylab="Proportion of genes", cex.axis=1.3, cex.lab=1.3)
dev.off()


pvec <- mat.or.vec(ncol(bitab), 3) 
for(i in 1:nrow(counts)) {
	others <- c(sum(counts[-i,1]+counts[-i,2]), sum(counts[-i,3]+counts[-i,4]))
	test <- c(sum(counts[i,1]+counts[i,2]), sum(counts[i,3]+counts[i,4]))
	pvec[i,1] <- fisher.test(cbind(test, others))$p.value
	pvec[i,2] <- test[1]/sum(test)
	pvec[i,3] <- others[1]/sum(others)
}
rownames(pvec) <- colnames(bitab)
colnames(pvec) <- c("FisherP", "MA_tissue", "MA_others")	
	
plim <- 0.05/nrow(pvec)
#pvec[pvec[,1]<plim,]
pvec["TESTIS",1]*nrow(pvec)
#0.0002392241



counts <- mat.or.vec(nrow(bitab), 6)
colnames(counts) <- c("BISTR", "BIW", "MOSTR", "MOW", "QU", "XNA") 
rownames(counts) <- rownames(bitab)
for(i in 1:nrow(bitab)) {
	counts[i,1] <- length(which(bitab2[i,]=="+" & !is.na(bitab2[i,]) & bitab[i,]!="?"))
	counts[i,2] <- length(which(bitab2[i,]=="?" & bitab[i,]=="+" & !is.na(bitab2[i,])  & !is.na(bitab[i,])))
	counts[i,3] <- length(which(bitab2[i,]=="-" & !is.na(bitab2[i,]) & bitab[i,]!="?"))
	counts[i,4] <- length(which(bitab2[i,]=="?" & bitab[i,]=="-" & !is.na(bitab2[i,]) & !is.na(bitab[i,])))
	counts[i,5] <- length(which(bitab[i,]=="?" & !is.na(bitab[i,])))
	counts[i,6] <- length(which(!is.na(bitab[i,])))
}

length(which(counts[,2]>0))
#34

s3 <- order((counts[,3] + counts[,4])/counts[,6], decreasing=T)
#X11(width=7, height=5)
pdf("rplots/biallelic.impr.barplot.genes.pdf", width=10, height=7)
par(mar=c(7.5,5,3,2))
barplot(t(counts[,c(1,2,4,3,5)])[,s3], las=3, col=c("yellow", "goldenrod2", "red", "red3", "grey"), ylab="Tissues", ylim=c(0,55), cex.axis=1.3, cex.lab=1.3)
legend("topright", fill=c("red3", "red", "goldenrod2", "yellow", "grey"), legend=c("Strong MAE", "MAE", "BAE", "Strong BAE", "?"), ncol=3, bty="n")
dev.off()

bipr <- (apply(counts[,3:4], 1, sum)/counts[,6])[counts[,6]>=15]
#X11(height=4, width=5)
pdf("rplots/monoa_prop_genehist.pdf", height=4, width=5)
hist(bipr, col="red", breaks=10, xlab="Proportion of monoallelic tissues", ylab="Genes", cex.axis=1.5, cex.lab=1.5, main="") 
dev.off()


known <- read.table("biallelic.weak.known.tab.txt", as.is=T, row.names=1, header=T)
int <- intersect(rownames(known), rownames(counts))
sd <- setdiff(rownames(counts), rownames(known))
ke <- counts[int,]
ne <- counts[sd,]
knum <- apply(ke, 1, function(x) {sum(x[3:4])/sum(x[1:4])})
nnum <- apply(ne, 1, function(x) {sum(x[3:4])/sum(x[1:4])})
#X11(width=8, height=4)
pdf("rplots/imprprop_known_novel.pdf", width=8, height=4)
par(mfrow=c(1,2))
hist(knum, breaks=seq(0, 1, length.out=10), col="grey", main="Known imprinted", xlab="Proportion of IMP+cIPM tissues", ylab="Genes", cex.axis=1.5, cex.lab=1.5)
hist(nnum, breaks=seq(0, 1, length.out=10), col="grey", main="Novel imprinted", xlab="Proportion of IMP+cIPM tissues", ylab="Genes", cex.axis=1.5, cex.lab=1.5)
dev.off()



### Known imprinted genes

bitab <- read.table("biallelic.weak.known.tab.txt", as.is=T, row.names=1, header=T)
bitab2 <- read.table("biallelic.strong.known.tab.txt", as.is=T, row.names=1, header=T)
gen <- read.table("ensg2nametab.txt", as.is=T)
bitab <- bitab[intersect(gen[,2], rownames(bitab)),]
bitab2 <- bitab2[rownames(bitab),]

counts <- mat.or.vec(ncol(bitab), 6)
colnames(counts) <- c("BISTR", "BIW", "MOSTR", "MOW", "QU", "XNA") 
rownames(counts) <- colnames(bitab)
for(i in 1:ncol(bitab)) {
	counts[i,1] <- length(which(bitab2[,i]=="+" & !is.na(bitab2[,i]) & bitab[,i]!="?"))
	counts[i,2] <- length(which(bitab2[,i]=="?" & bitab[,i]=="+" & !is.na(bitab2[,i])  & !is.na(bitab[,i])))
	counts[i,3] <- length(which(bitab2[,i]=="-" & !is.na(bitab2[,i]) & bitab[,i]!="?"))
	counts[i,4] <- length(which(bitab2[,i]=="?" & bitab[,i]=="-" & !is.na(bitab2[,i]) & !is.na(bitab[,i])))
	counts[i,5] <- length(which(bitab[,i]=="?" & !is.na(bitab[,i])))
	counts[i,6] <- length(which(!is.na(bitab[,i])))
}

s3 <- order((counts[,3] + counts[,4])/counts[,6], decreasing=T)
#X11(width=7, height=5)
pdf("rplots/biallelic.known.barplot.pdf", width=10, height=7)
par(mar=c(7,5,3,2))
barplot(t(apply(counts[,c(1,2,4,3,5)], 2, function(x) {x/counts[,6]}))[,s3], las=3, col=c("yellow", "goldenrod2", "red", "red3", "grey"), ylab="Proportion of genes", cex.axis=1.3, cex.lab=1.3)
#legend("topright", fill=c("red", "yellow"), legend=c("Monoall", "Biall"), ncol=2, bty="n", cex=1.2)
dev.off()



##Geneplots

final <- read.table("biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
sd <- setdiff(rownames(bitab), rownames(final))
bitab <- bitab[sd,]
bitab2 <- bitab2[sd,]

counts <- mat.or.vec(nrow(bitab), 6)
colnames(counts) <- c("BISTR", "BIW", "MOSTR", "MOW", "QU", "XNA") 
rownames(counts) <- rownames(bitab)
for(i in 1:nrow(bitab)) {
	counts[i,1] <- length(which(bitab2[i,]=="+" & !is.na(bitab2[i,]) & bitab[i,]!="?"))
	counts[i,2] <- length(which(bitab2[i,]=="?" & bitab[i,]=="+" & !is.na(bitab2[i,])  & !is.na(bitab[i,])))
	counts[i,3] <- length(which(bitab2[i,]=="-" & !is.na(bitab2[i,]) & bitab[i,]!="?"))
	counts[i,4] <- length(which(bitab2[i,]=="?" & bitab[i,]=="-" & !is.na(bitab2[i,]) & !is.na(bitab[i,])))
	counts[i,5] <- length(which(bitab[i,]=="?" & !is.na(bitab[i,])))
	counts[i,6] <- length(which(!is.na(bitab[i,])))
}


s3 <- order(counts[,6], decreasing=F)
#X11(width=10, height=5)
pdf("rplots/biallelic.known.barplot.genes.pdf", width=12, height=7)
par(mar=c(7.5,5,3,2))
barplot(t(counts[,c(1,2,4,3,5)])[,s3], las=3, col=c("yellow", "goldenrod2", "red", "red3", "grey"), ylab="Tissues", ylim=c(0,45), cex.axis=1.3, cex.lab=1.3)
legend("topright", fill=c("red3", "red", "goldenrod2", "yellow", "grey"), legend=c("Strong MAE", "MAE", "BAE", "Strong BAE", "?"), ncol=3, bty="n")
dev.off()





########   Maternal/paternal differences  #########

matpat <- read.table("impr_known_matpat.txt", as.is=T, header=T, row.names=1)
matpat <- matpat[rownames(matpat)!="IGF2" & rownames(matpat)!="GRB10",]
mae <- rownames(matpat)[matpat[,3]=="Maternal"]
pae <- rownames(matpat)[matpat[,3]=="Paternal"]

bitab <- read.table("biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
bitab2 <- read.table("biallelic.strong.final.tab.txt", as.is=T, row.names=1, header=T)
gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
bitab <- bitab[intersect(gen[,2], rownames(bitab)),]
bitab2 <- bitab2[rownames(bitab),]

#bitab <- bitab[,colnames(bitab)!="TESTIS"]
#bitab2 <- bitab2[,colnames(bitab2)!="TESTIS"]

counts <- mat.or.vec(ncol(bitab), 6)
colnames(counts) <- c("MATM", "MATB", "PATM", "PATB", "MATE", "PATE") 
rownames(counts) <- colnames(bitab)
for(i in 1:ncol(bitab)) {
	counts[i,1] <- length(intersect(rownames(bitab2)[bitab[,i]=="-" & !is.na(bitab[,i]) & bitab[,i]!="?"], mae))
	counts[i,2] <- length(intersect(rownames(bitab2)[bitab[,i]=="+" & !is.na(bitab[,i]) & bitab[,i]!="?"], mae))
	counts[i,3] <- length(intersect(rownames(bitab2)[bitab[,i]=="-" & !is.na(bitab[,i]) & bitab[,i]!="?"], pae))
	counts[i,4] <- length(intersect(rownames(bitab2)[bitab[,i]=="+" & !is.na(bitab[,i]) & bitab[,i]!="?"], pae))
	counts[i,5] <- length(intersect(rownames(bitab2)[!is.na(bitab2[,i])], mae))
	counts[i,6] <- length(intersect(rownames(bitab2)[!is.na(bitab2[,i])], pae))
}



#X11(width=7, height=5)
pdf("rplots/biallelic.impr.biall.barplot1_matpat.pdf", width=10, height=7)
par(mar=c(8,5,3,2))
s3 <- order((counts[,1]+counts[,2])/apply(counts[,1:4], 1, sum), decreasing=T)
barplot(t(counts[,c(1:4)]/apply(counts[,1:4], 1, sum))[,s3], las=3, col=c( "red", "red4", "royalblue2", "royalblue4"), ylab="Proportion of genes", cex.axis=1.5, cex.lab=1.5, cex.names=1.3, ylim=c(0,1.2))
legend("topright", fill=c("red", "red4", "royalblue3", "royalblue4"), legend=c("MAE Mat Expr", "BAE Mat Expr", "MAE Pat Expr", "BAE Pat Expr"), ncol=4, bty="n", cex=1.3)
dev.off()




counts <- mat.or.vec(nrow(bitab), 3)
colnames(counts) <- c("MAE", "BAE", "TOT") 
rownames(counts) <- rownames(bitab)
for(i in 1:nrow(bitab)) {
	counts[i,1] <- length(which(bitab[i,]=="-" & !is.na(bitab[i,])))
	counts[i,2] <- length(which(bitab[i,]=="+" & !is.na(bitab[i,])))
	counts[i,3] <- length(which(!is.na(bitab[i,])))
}
##1 = maternal expression, 0 = paternal expression
dirvec <- c(rep(1, length(mae)), rep(0, length(pae)))
names(dirvec) <- c(mae, pae)
int <- intersect(names(dirvec), rownames(counts))
dirvec <- dirvec[int]
countsf <- counts[int,]

tmp1 <- countsf[dirvec==1,1:2]
tmp1 <- tmp1[order(tmp1[,2]/(tmp1[,2]+tmp1[,1]), decreasing=T),]
tmp2 <- countsf[dirvec==0,1:2]
tmp2 <- tmp2[order(tmp2[,2]/(tmp2[,2]+tmp2[,1]), decreasing=T),]
posvec <- c(rep(0.2, nrow(tmp1)), 2.5, rep(0.2, nrow(tmp2)-1))
pdf("rplots/biallelic.patmat.pergene.barplot.merged.pdf", height=5, width=7)
#X11(height=5, width=7)
par(mar=c(7,4.5,2,2))
barplot(t(rbind(tmp1, tmp2)), col=c("red", "yellow"), las=3, ylab="Tissues", cex.axis=1.5, cex.lab=1.5, space=posvec, ylim=c(0,45))
text(x=5.5, y=37, labels="Maternal", cex=1.4)
text(x=25, y=37, labels="Paternal", cex=1.4)
legend("topright", fill=c("red", "yellow"), legend=c("Imprinted", "Biallelic"), cex=1.1, bty="n", ncol=2)
dev.off()

##permutation test
np <- 100000
pres <- mat.or.vec(np, 2)
colnames(pres) <- c("pi_perm", "mi_perm")
for(i in 1:np) {
	s <- sample(c(rep(1, length(which(dirvec==1))), rep(0, length(which(dirvec==0)))))
	pres[i,1] <- sum(countsf[s==1,1])/sum(countsf[s==1,3])
	pres[i,2] <- sum(countsf[s==0,1])/sum(countsf[s==0,3])
}
pi_emp <- sum(countsf[dirvec==1,1])/sum(countsf[dirvec==1,3])
mi_emp <- sum(countsf[dirvec==0,1])/sum(countsf[dirvec==0,3])
length(which(pi_emp>pres[,1]))/np
#p=0.06252







tau <- read.table("impr_tau_tab.txt", as.is=T, sep="\t", header=T, row.names=1)

matpat <- read.table("impr_known_matpat.txt", as.is=T, header=T, row.names=1)
matpat <- matpat[rownames(matpat)!="IGF2" & rownames(matpat)!="GRB10",]
mae <- intersect(rownames(tau), rownames(matpat)[matpat[,3]=="Maternal"])
pae <- intersect(rownames(tau), rownames(matpat)[matpat[,3]=="Paternal"])


tau_mp <- tau[c(mae, pae),]
tau_m <- tau[mae,]
tau_p <- tau[pae,]


X11(height=4, width=8)
pdf("rplots/mat_pat_tau.pdf", height=4, width=8)
par(mfrow=c(1,2))
hist(unlist(tau_m), main="Maternally expr", xlab="tau", col="grey")
hist(unlist(tau_p), main="Paternally expr", xlab="tau", col="grey")
dev.off()


lm <- length(mae)
lp <- length(pae)
np <- 1000
pres <- mat.or.vec(np, 2)
colnames(pres) <- c("pi_perm", "mi_perm")
for(i in 1:np) {
	s <- sample(c(mae, pae))
	tau_m_p <- tau[s[1:lm],]
	tau_p_p <- tau[s[-c(1:lm)],]
	pres[i,1] <- wilcox.test(unlist(tau_m_p), unlist(tau_p_p))$p.value
}
p_obs <- wilcox.test(unlist(tau_m), unlist(tau_p))$p.value










###########################################################################################################################################
#                                                          Gene expression                                                                #
###########################################################################################################################################



#from median expression level of each gene per tissue
cd ~/tuuli_lab/tlappalainen/imprint/analysis/expr_revision_2
cp ~/tuuli_lab/tlappalainen/gtex/data/quantification/gene_tissue/RPKM_GeneLevel_December_analysis_freeze_tissuemedian.txt .


#########  Heatmap  #########  

datal <- vector("list", 5)
datal[[1]] <- read.table("expr/RPKM_GeneLevel_December_analysis_freeze_tissuemedian.txt", as.is=T, header=T, sep="\t")
datal[[2]] <- read.table("expr/GENCORD_GENE_RPKM_F.median", as.is=T, header=T, sep="\t")
datal[[3]] <- read.table("expr/GENCORD_GENE_RPKM_L.median", as.is=T, header=T, sep="\t")
datal[[4]] <- read.table("expr/GENCORD_GENE_RPKM_T.median", as.is=T, header=T, sep="\t")
datal[[5]] <- read.table("expr/GD462.GeneQuantRPKM.50FN.median", as.is=T, header=T, sep="\t", fill=T)
#combine brain
annot <- datal[[1]][,1:2]
nobr <- datal[[1]][,grep("BRN", colnames(datal[[1]]), invert=T)]
nobr <- nobr[,3:ncol(nobr)]
br <- datal[[1]][,grep("BRN", colnames(datal[[1]]))]
BRAIN <- apply(br, 1, median)
all <- cbind(nobr, BRAIN)
all <- all[,order(colnames(all))]
datal[[1]] <- cbind(annot, all)
datal[[5]] <- datal[[5]][,-c(3,4)]
tn <- c("GC-FIBRBLS", "GC-LCL", "GC-TELL", "GD-LCL")
for(i in 1:length(datal)) {
	tmp <- as.matrix(datal[[i]][,3:ncol(datal[[i]])])
	n <- unlist(lapply(strsplit(datal[[i]][,1], ".", fixed=T), function(x) {x[1]}))
	rownames(tmp) <- n
	datal[[i]] <- tmp
}
for(i in 2:5) {
	colnames(datal[[i]]) <- tn[i-1]
}
genecount <- c(apply(datal[[1]], 2, function(x) {length(x[x!=0 & !is.na(x)])}), unlist(lapply(datal[2:5], function(x) {length(x[,1][x[,1]!=0 & !is.na(x)])})))
names(genecount)[(ncol(datal[[1]])+1):(ncol(datal[[1]])+4)] <- tn
datalo <- datal

gen <- read.table("final_impr_genelist_ensg2name.txt", as.is=T)
g <- gen[,2]
ge <- gen[,1]
known <- read.table("known_otago_ensg2name.txt", as.is=T)[,1]
status <- read.table("final_knownstatus.txt", row.names=1, sep="\t", as.is=T)
putative <- rownames(status)[status=="Provisional"]
kn <- rownames(status)[status=="Known"]
ensg <- read.table("ensg2nametab.txt", as.is=T, sep="\t", row.names=1)

for(i in 1:length(datal)) {
	int <- intersect(ge, rownames(datal[[i]]))
	datal[[i]] <- as.data.frame(datal[[i]][int,])
}
for(i in 2:5) {
	colnames(datal[[i]]) <- tn[i-1]
}

data <- matrix(length(ge), sum(unlist(lapply(datal, ncol))), data=NA)
rownames(data) <- ge
colnames(data) <- unlist(lapply(datal,  colnames))
for(i in 1:length(datal)) {
	for(k in 1:ncol(datal[[i]])) {
		int <- intersect(rownames(datal[[i]]), ge)
		data[int, colnames(datal[[i]])[k]] <- datal[[i]][int,k]
	}
}
data <- data[apply(data, 1, function(x) {!all(is.na(x))}),]

tis <- colnames(data)
gn <- g
for(i in 1:length(g)) {
	if (any(putative==g[i])) {
		gn[i] <- paste(g[i], "*", sep="")
	}
	if (any(kn==g[i])) {
		gn[i] <- paste(g[i], "**", sep="")
	}
}

tabs <- data
rownames(tabs) <- g
write.table(tabs, "impr_expr_tab.txt", sep="\t", quote=F)
gn[gn=="INPP5F_V2**"] <- "INPP5F**"
rownames(tabs) <- gn

pos <- read.table("/data/research/tuuli_lab/tlappalainen/resource/gencode.v12.annotation.tab.gene", as.is=T)
rownames(pos) <- unlist(lapply(strsplit(pos[,9], ".", fixed=T), function(x) {x[1]}))
ensgi <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
ensgi <- ensgi[!duplicated(ensgi[,2]),]
row.names(ensgi) <- ensgi[,2]
poss <- pos[ensgi[g,1],]
tss <- poss[,4]
tss[poss[,7]=="-"] <- poss[,5][poss[,7]=="-"]
p <- paste(poss[,1], round(tss/1000000, 1), sep=":")
c <- as.numeric(substr(poss[,1], 4, nchar(poss[,1])))
p2 <- paste(c, round(tss/1000000, 1), sep=":")
s <- order(c, tss, decreasing=T)
#	X11(height=9, width=8)
pdf(paste("rplots/heatmap_impr_expr_logrpkm.pdf", sep=""), height=9, width=8)
par(mar=c(7,8,3,6))
ptab <- tabs[s,]
natab <- tabs[s,]
natab[is.na(tabs[s,])] <- 1
natab[!is.na(tabs[s,])] <- NA
image(log2(max(ptab, na.rm=T))-log2(t(ptab)), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, main="log2 RPKM", col=heat.colors(20))
#image(t(natab), add=T, xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.5, cex.lab=1.5, col="grey")
axis(at=0:(nrow(t(ptab))-1)/(nrow(t(ptab))-1), labels=rownames(t(ptab)), side=1, las=3)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=colnames(t(ptab)), side=2, las=1)
axis(at=0:(ncol(t(ptab))-1)/(ncol(t(ptab))-1), labels=p[s], side=4, las=1)		
dev.off()

pdf(paste("rplots/heatmap_impr_expr_legend.pdf", sep=""))
plot(0,0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
mi <- log2(ptab)
mi <- min(mi[mi!="-Inf" & !is.na(mi)])
ma <- log2(ptab)
ma <- max(ma[ma!="Inf" & !is.na(ma)])
legend("topleft", fill=c(heat.colors(20), "white"), legend=rep("", 21), bty="n", x.intersp=-1, cex=2, ncol=21)
text(x=seq(-0.9, 0.9, length.out=21), y=rep(0.75, 21), las=3, pos=2, c(round(seq(ma, mi, length.out=20), 2), "Not Expr"), srt=90, cex=1.5)
dev.off()




#########  Boxplot  #########  

spl <- vector("list", ncol(tabs))
names(spl) <- colnames(tabs)
for(i in 1:ncol(tabs)) {
	tmp <- tabs[,i]
	tmp <- replace(tmp, is.na(tmp), 0)
	tmp <- tmp+0.001
	spl[[i]] <- tmp
}
med <- unlist(lapply(spl, function(x) {median(log2(x)[log2(x)!="-Inf"], na.rm=T)}))
s <- order(med, decreasing=T)

l <- vector("list", length(spl)*3)
n <- vector("character", length(spl)*3)
pos <- vector("numeric", length(spl)*3)
spls <- spl[s]
spls <- spls[grep("GC.", names(spls), invert=T)]
spls <- spls[grep("GD.", names(spls), invert=T)]
spls <- lapply(spls, function(x) {x+0.001})
spls <- lapply(spls, function(x) {replace(x, is.na(x), 0.001)})

#X11(height=4, width=8)
pdf("rplots/tissue_rpkm_boxplot.pdf", height=4, width=8)
par(mar=c(7,4,4,4))
boxplot(lapply(spls, function(x) {log2(x)[log2(x)!="-Inf"]}), las=3, col="grey", ylab="log2 RPKM", cex.axis=1.2, cex.lab=1.3)
dev.off()


##novel genes more tissue-specific?
rownames(tabs) <- gn
ke <- tabs[grep("**", rownames(tabs), fixed=T),]
ne <- tabs[grep("**", rownames(tabs), fixed=T, invert=T),]
lim <- 0.1
knum <- apply(ke, 1, function(x) {length(which(x>lim))})
nnum <- apply(ne, 1, function(x) {length(which(x>lim))})
pdf("rplots/expr_known_novel.pdf", width=8, height=4)
par(mfrow=c(1,2))
hist(knum, breaks=seq(0, 37, length.out=10), col="grey", main="Known imprinted", xlab="Tissues with RPKM>0.1", ylab="Genes", cex.axis=1.5, cex.lab=1.5)
hist(nnum, breaks=seq(0, 37, length.out=10), col="grey", main="Novel imprinted", xlab="Tissues with RPKM>0.1", ylab="Genes", cex.axis=1.5, cex.lab=1.5)
dev.off()





#########  imprinted status vs expression level  #########  


expr <- read.table("impr_expr_tab.txt", header=T, row.names=1, check.names=F)
colnames(expr)[colnames(expr)=="GC-FIBRBLS"] <- "GC_FIBRBLS"
colnames(expr)[colnames(expr)=="GC-TELL"] <- "GC_TCELL"
colnames(expr)[colnames(expr)=="GC-LCL"] <- "GC_LCL"
colnames(expr)[colnames(expr)=="GD-LCL"] <- "GD_LCL"
bitab <- read.table("biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
expr <- expr[rownames(bitab),]
expr <- expr[,colnames(bitab)]
elist <- vector("list", nrow(expr))
names(elist) <- rownames(expr)
for(i in 1:nrow(expr)) {
	elist[[i]] <- vector("list", 2)
	elist[[i]][[1]] <- expr[i,colnames(bitab)[bitab[i,]=="+" & !is.na(bitab[i,])]]
	elist[[i]][[2]] <- expr[i,colnames(bitab)[bitab[i,]=="-" & !is.na(bitab[i,])]]
	elist[[i]][[1]] <- elist[[i]][[1]][!is.na(elist[[i]][[1]])]
	elist[[i]][[2]] <- elist[[i]][[2]][!is.na(elist[[i]][[2]])]
}
#1 = biallelic, 2 is imprinted
tab <- mat.or.vec(nrow(expr), 6)
rownames(tab) <- rownames(expr)
colnames(tab) <- c("N_BAE", "N_MAE", "MEDIAN_BAE", "MEDIAN_MAE", "P_MW", "P_MW_DOUBLE")
for(i in 1:nrow(expr)) {
	tab[i,1] <- length(elist[[i]][[1]])	
	tab[i,2] <- length(elist[[i]][[2]])	
	tab[i,3] <- median(as.numeric(elist[[i]][[1]]), na.rm=T)	
	tab[i,4] <- median(as.numeric(elist[[i]][[2]]), na.rm=T)
	if (tab[i,1]>0 & tab[i,2]>0) {
		tab[i,5] <- wilcox.test(as.numeric(elist[[i]][[1]]), as.numeric(elist[[i]][[2]]))$p.value
		tab[i,6] <- wilcox.test(as.numeric(elist[[i]][[1]])*2, as.numeric(elist[[i]][[2]]))$p.value
	} else {
		tab[i,5] <- NA
		tab[i,6] <- NA
	}
}

elist <- elist[!is.na(tab[,3]) & !is.na(tab[,4])]
posvec <- vector("numeric", length(elist))
posvec[1] <- 1
for(i in 2:length(elist)) {
	posvec[i] <- posvec[i-1]+5
}

val <- c(unlist(lapply(elist, function(x) {log2(x[[1]]+0.001)})), unlist(lapply(elist, function(x) {log2(x[[2]]+0.001)})))
val <- val[val!="-Inf" & val!="Inf"]
s <- order(unlist(lapply(elist, function(x) {median(log2(x[[1]]+0.001))})) - unlist(lapply(elist, function(x) {median(log2(x[[2]]+0.001))})))
mi <- min(val)
ma <- max(val)
#X11(height=5, width=15)
pdf("rplots/expr_bae_mae_pergene.pdf", height=5, width=15)
par(mar=c(10,4,2,2))
boxplot(lapply(elist, function(x) {log2(x[[1]]+0.001)})[s], at=posvec, col="yellow", las=3, border="goldenrod3", ylim=c(mi, ma), ylab="log2 RPKM", names=names(elist)[s], cex.axis=1.3, cex.lab=1.3, cex=0.5, lwd=1.2)
boxplot(lapply(elist, function(x) {log2(x[[2]]+0.001)})[s], at=posvec+1, col="red", add=T, names=rep("", length(elist)), border="red3", cex.axis=1.3, cex.lab=1.3, cex=0.5, lwd=1.2)
points(y=unlist(lapply(elist, function(x) {median(log2((x[[2]]+0.001)*2))}))[s], x=posvec, col="black", pch=3, lwd=1)
legend("topright", ncol=2, fill=c("yellow", "red"), legend=c("Biallelic", "Imprinted"), bty="n", cex=1.3)
dev.off()


tab <- tab[!is.na(tab[,3]) & !is.na(tab[,4]),]
tabf <- tab[tab[,1]>1 & tab[,2]>1,]

val <- c(tabf[,3], tabf[,4])
val <- val[!is.na(val)]
val <- log2(val)
val <- val[val!="-Inf" & val!="Inf"]
mi <- min(val)-1.5
ma <- max(val)+1
#X11(height=5, width=5)
pdf("rplots/expr_bae_mae_pergene_scatter_multitis.pdf", height=5, width=5)
plot(log2(tabf[,4]+0.001), log2(tabf[,3]+0.001), lwd=2, col="white", xlab="Median log2 RPKM, IMP tissues", ylab="Median log2 RPKM, BI tissues", cex.axis=1.3, cex.lab=1.3, xlim=c(mi, ma), ylim=c(mi, ma))
text(log2(tabf[,4]+0.001), log2(tabf[,3]+0.001), labels=rownames(tabf))
abline(0,1)
dev.off()


val <- c(tabf[,3], tabf[,4])
val <- val[!is.na(val)]
val <- log2(val)
val <- val[val!="-Inf" & val!="Inf"]
mi <- min(val)-1.5
ma <- max(val)+1
#X11(height=5, width=5)
pdf("rplots/expr_bae_mae_pergene_scatter_noname_multitis.pdf", height=5, width=5)
plot(log2(tabf[,4]+0.001), log2(tabf[,3]+0.001), lwd=2, col="white", xlab="Median log2 RPKM, IMP tissues", ylab="Median log2 RPKM, BI tissues", cex.axis=1.3, cex.lab=1.3, xlim=c(mi, ma), ylim=c(mi, ma))
points(log2(tabf[,4]+0.001)[log2(tabf[,4]+0.001)>log2(tabf[,3]+0.001)], log2(tabf[,3]+0.001)[log2(tabf[,4]+0.001)>log2(tabf[,3]+0.001)], col="red", lwd=2)
points(log2(tabf[,4]+0.001)[log2(tabf[,3]+0.001)>=log2(tabf[,4]+0.001)], log2(tabf[,3]+0.001)[log2(tabf[,3]+0.001)>=log2(tabf[,4]+0.001)], col="Goldenrod", lwd=2)
abline(0,1)
dev.off()


wilcox.test(tab[,3], tab[,4], paired=T)
#p-value = 0.6978
wilcox.test(tabf[,3], tabf[,4], paired=T)
#p=0.6507






###########################################################################################################################################
#                                                  Interindividual variation                                                              #
###########################################################################################################################################




datas <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.final_impr.txt", as.is=T, header=T, sep="\t")
data <- datas[,c("TISSUE", "SEVERE_GENE", "SUBJECT", "REF_COUNT", "NONREF_COUNT")]
tmp <- apply(data[,c(1,2,3)], 1, function(x) {paste(x, sep="_", collapse="_")})
spl <- split(data, tmp)
l <- vector("list", length(spl))
for(i in 1:length(spl)) {
	tmp <- apply(spl[[i]][,4:5], 1, function(x) {abs(0.5-(x[1]/(x[1]+x[2])))})
	vec <- c(spl[[i]][1,1], spl[[i]][1,2], spl[[i]][1,3], median(tmp))
	l[[i]] <- matrix(ncol=4, nrow=1, data=vec)
}
data <- as.data.frame(do.call(rbind, l))
ens2g <- read.table("ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
spl <- split(data, as.character(data[,2]))

##heterogeneity across genes
genes <- unique(data[,2])
genei <- vector("list", length(genes))
names(genei) <- genes
for(k in 1:length(genes)) {
	ge <- genes[k]
	td1 <- spl[[ge]]
	spli <- split(td1, td1[,3])
	splif <- spli[unlist(lapply(spli, nrow))>5]
	td2 <- do.call(rbind, splif)
	if (length(dim(td2))>0) {
		spli <- split(td2, td2[,1])
		splif <- spli[unlist(lapply(spli, nrow))>5]
		td2 <- do.call(rbind, splif)
		if (length(dim(td2))>0) {
			spli <- split(td2, td2[,3])
			genei[[k]] <- unlist(lapply(spli, function(x) {median(as.numeric(x[,4], na.rm=T))}))
			genei[[k]] <- genei[[k]][!is.na(genei[[k]])]
		}
	}
}
genei <- genei[unlist(lapply(genei, length))>0]
alli <- unlist(genei)
spli <- split(alli, names(alli))
spli <- spli[unlist(lapply(spli, length))>=10]

tmpl <- vector("list", length(genei))
for(i in 1:length(genei)) {
	tmpl[[i]] <- cbind(rep(names(genei)[i], length(genei[[i]])), names(genei[[i]]), genei[[i]])
}
td1 <- as.data.frame(do.call(rbind, tmpl))
spli <- split(td1, td1[,2])
splif <- spli[unlist(lapply(spli, nrow))>=10]
td2 <- do.call(rbind, splif)
spli <- split(td2, td2[,1])
splif <- spli[unlist(lapply(spli, nrow))>=10]
td2 <- do.call(rbind, splif)


##heatmap of tissues x samples
gen <- as.character(unique(td2[,1]))
sam <- as.character(unique(td2[,2]))
tab <- matrix(nrow=length(gen), ncol=length(sam), data=NA)
rownames(tab) <- gen
sam2 <- unlist(lapply(strsplit(sam, "-"), function(x) {x[2]}))
colnames(tab) <- sam
for(i in 1:nrow(td2)) {
	tab[as.character(td2[i,1]),as.character(td2[i,2])] <- as.numeric(as.vector(td2[i,3]))
}
colnames(tab) <- sam2
rownames(tab) <- ens2ginv[gen,2]
log <- apply(tab, 1, function(x) {length(which(is.na(x)))})
tab <- tab[log<=40,]
natab <- tab
natab[is.na(tab)] <- 1
natab[!is.na(tab)] <- NA
s <- order(apply(tab, 2, function(x) {median(x, na.rm=T)}))
tab <- tab[,s]
natab <- natab[,s]
s <- order(apply(tab, 1, function(x) {median(x, na.rm=T)}))
tab <- tab[s,]
natab <- natab[s,]
pdf(paste("rplots/het_ind_gene_snp_heatmap.pdf", sep=""), height=7, width=14)
par(mar=c(6,10,3,3))
image((max(tab, na.rm=T))-(t(tab)), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.7, cex.lab=1.7, col=heat.colors(20))
image(t(natab), add=T, xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.7, cex.lab=1.7, col="grey")
axis(at=0:(nrow(t(tab))-1)/(nrow(t(tab))-1), labels=rownames(t(tab)), side=1, las=3, cex.axis=1.3, cex.lab=1.3)
axis(at=0:(ncol(t(tab))-1)/(ncol(t(tab))-1), labels=colnames(t(tab)), side=2, las=1, cex.axis=1.3, cex.lab=1.3)
dev.off()



##from the analysis by gene:
hetg <- c("PPIEL", "KCNQ1", "ZNF331") 

tabsave <- vector("list", nrow(tab))
##all genes with lots of data
##look if loss of imprinting is shared between tissues
genes <- hetg
for(k in 1:length(genes)) {
g <- genes[k]
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
i=1
td1 <- spl[[ge[i]]]
spli <- split(td1, td1[,3])
splif <- spli[unlist(lapply(spli, nrow))>7]
td2 <- do.call(rbind, splif)
spli <- split(td2, td2[,1])
splif <- spli[unlist(lapply(spli, nrow))>5]
td2 <- do.call(rbind, splif)
if (length(dim(td2))>0) {

##heatmap of tissues x samples
tis <- unique(as.character(td2[,1]))
sam <- unique(as.character(td2[,3]))
tab <- matrix(nrow=length(tis), ncol=length(sam), data=NA)
rownames(tab) <- tis
sam2 <- unlist(lapply(strsplit(sam, "-"), function(x) {x[2]}))
colnames(tab) <- sam
for(i in 1:nrow(td2)) {
	tab[as.character(td2[i,1]),as.character(td2[i,3])] <- as.numeric(as.vector(td2[i,4]))
}
colnames(tab) <- sam2
natab <- tab
natab[is.na(tab)] <- 1
natab[!is.na(tab)] <- NA
s <- order(apply(tab, 2, function(x) {median(x, na.rm=T)}))
tab <- tab[,s]
natab <- natab[,s]
s <- order(apply(tab, 1, function(x) {median(x, na.rm=T)}))
tab <- tab[s,]
natab <- natab[s,]
tabsave[[k]] <- tab
pdf(paste("rplots/het_", g, "_ind_tis_snp_heatmap.pdf", sep=""), height=7, width=10)
par(mar=c(7,7,3,3))
image((max(tab, na.rm=T))-(t(tab)), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.7, cex.lab=1.7, col=heat.colors(20))
image(t(natab), add=T, xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.7, cex.lab=1.7, col="grey")
axis(at=0:(nrow(t(tab))-1)/(nrow(t(tab))-1), labels=rownames(t(tab)), side=1, las=3, cex.axis=1.3, cex.lab=1.3)
axis(at=0:(ncol(t(tab))-1)/(ncol(t(tab))-1), labels=colnames(t(tab)), side=2, las=1, cex.axis=1.3, cex.lab=1.3)
dev.off()

}}


##for main figure

g <- c("ZNF331")
ge <- ens2g[g,1]
ens2ginv <- ens2g
rownames(ens2ginv) <- ens2ginv[,1]
ge <- intersect(ge, names(spl))
g <- ens2ginv[ge,2]
i=1
td1 <- spl[[ge[i]]]
spli <- split(td1, td1[,3])
splif <- spli[unlist(lapply(spli, nrow))>=10]
td2 <- do.call(rbind, splif)
spli <- split(td2, td2[,1])
splif <- spli[unlist(lapply(spli, nrow))>5]
td2 <- do.call(rbind, splif)

##heatmap of tissues x samples
tis <- unique(as.character(td2[,1]))
sam <- unique(as.character(td2[,3]))
tab <- matrix(nrow=length(tis), ncol=length(sam), data=NA)
rownames(tab) <- tis
sam2 <- unlist(lapply(strsplit(sam, "-"), function(x) {x[2]}))
colnames(tab) <- sam
for(i in 1:nrow(td2)) {
	tab[as.character(td2[i,1]),as.character(td2[i,3])] <- as.numeric(td2[i,4])
}
colnames(tab) <- sam2
natab <- tab
natab[is.na(tab)] <- 1
natab[!is.na(tab)] <- NA
s <- order(apply(tab, 2, function(x) {median(x, na.rm=T)}))
tab <- tab[,s]
natab <- natab[,s]
s <- order(apply(tab, 1, function(x) {median(x, na.rm=T)}))
tab <- tab[s,]
natab <- natab[s,]
pdf(paste("rplots/ZNF331_ind_tis_snp_heatmap_main.pdf", sep=""), height=7, width=7)
par(mar=c(7,7,3,3))
image((max(tab, na.rm=T))-(t(tab)), xaxt="n", xlab="",  yaxt="n", ylab="", cex.axis=1.7, cex.lab=1.7, col=heat.colors(20))
image(t(natab), add=T, xaxt="", xlab="",  yaxt="n", ylab="", cex.axis=1.7, cex.lab=1.7, col="grey")
axis(at=0:(nrow(t(tab))-1)/(nrow(t(tab))-1), labels=rownames(t(tab)), side=1, las=3, cex.axis=1.5, cex.lab=1.5)
axis(at=0:(ncol(t(tab))-1)/(ncol(t(tab))-1), labels=colnames(t(tab)), side=2, las=1, cex.axis=1.5, cex.lab=1.5)
dev.off()






###########################################################################################################################################
#                                                   Age versus sex effects                                                                #
###########################################################################################################################################





cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint

#########  correlation with overall level of imprinting  #########  

 
#Consistency with SNP results (to rule out imputation errors) has been checked and is OK

gtab <- read.table("final_impr_tab.txt", as.is=T)
gspl <- split(gtab[,3], gtab[,1])
data <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.final_impr.txt", as.is=T, header=T, sep="\t")
ensg <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
rownames(ensg) <- ensg[,1]
n <- ensg[data[,2],2]
data <- cbind(data, n)
spl <- split(data, n)
spl2 <- lapply(spl, function(x) {split(x, x[,1])})
for(i in 1:length(spl)) {
	tis <- gspl[[names(spl2)[i]]]
	int <- intersect(tis, names(spl2[[i]]))
	if (length(int)==0) {
		spl[[i]] <- NA
	} else {
		tmp <- spl2[[i]][int]
		tmp2 <- do.call(rbind, tmp)
		ispl <- split(tmp2, tmp2[,3])
		tmp3 <- do.call(rbind, lapply(ispl, function(x) {apply(x[,4:5], 2, sum)}))
		tmp4 <- do.call(rbind, lapply(ispl, function(x) {x[1,1:3]}))
		spl[[i]] <- cbind(tmp4, tmp3)
	}
}
spl <- spl[!unlist(lapply(spl, function(x) {any(is.na(x))}))]

data <- do.call(rbind, spl)
spl <- split(data, data[,3])
info <- read.table("../../data/gtex/GTEx.SubjectPhenotypesDS.v9.1.SUBJID_GENDER_AGE.txt", header=T, as.is=T, sep="\t")
info <- info[!duplicated(info[,1]),]
rownames(info) <- info[,1]
info <- rbind(info)
int <- intersect(names(spl), rownames(info))
spl <- spl[int]
info <- info[int,]
med <- unlist(lapply(spl, function(x) {median(0.5+abs(0.5-(x[,4]/(x[,4]+x[,5]))))}))

cor.test(info[,3], med, method="s")
#rho = -0.0547, p=0.5158

wilcox.test(med[info[,2]==1], med[info[,2]==2])
#p=0.802





#########  separately for tissues  #########  



gtab <- read.table("final_impr_tab.txt", as.is=T)
gspl <- split(gtab[,3], gtab[,1])
data <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.final_impr.txt", as.is=T, header=T, sep="\t")
ensg <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
rownames(ensg) <- ensg[,1]
n <- ensg[data[,2],2]
data <- cbind(data, n)
spl <- split(data, n)
spl2 <- lapply(spl, function(x) {split(x, x[,1])})
#get only imprinted tissues
for(i in 1:length(spl)) {
	tis <- gspl[[names(spl2)[i]]]
	int <- intersect(tis, names(spl2[[i]]))
	if (length(int)==0) {
		spl[[i]] <- NA
	} else {
		spl2[[i]] <- spl2[[i]][int]
	}
}

spl <- lapply(spl2, function(x) {do.call(rbind, x)})
data <- do.call(rbind, spl)
spl <- split(data, data[,1])
spl2 <- lapply(spl, function(x) {split(x, x[,3])})
spl3 <- lapply(spl2, function(x) {lapply(x, function(y) {split(y, y[,2])})})
med <- lapply(spl3, function(x) {lapply(x, function(y) {unlist(lapply(y, function(z) {
	median(0.5+abs(0.5-(z[,4]/(z[,4]+z[,5]))))}))})})
med2 <- lapply(med, function(x) {unlist(lapply(x, function(y) {median(y)}))})

info <- read.table("../../data/gtex/GTEx.SubjectPhenotypesDS.v9.1.SUBJID_GENDER_AGE.txt", header=T, as.is=T, sep="\t")
info <- info[!duplicated(info[,1]),]
rownames(info) <- info[,1]

tab <- matrix(data=NA, nrow=length(med2), ncol=5)
rownames(tab) <- names(med2)
colnames(tab) <- c("AGE_RHO", "AGE_P", "SEX_M-F", "SEX_P", "N")
for(i in 1:length(med2)) {
	int <- intersect(names(med2[[i]]), rownames(info))
	ratio <- med2[[i]][int]
	names(ratio) <- int
	if (length(ratio)>=20) {
		c <- cor.test(ratio, info[int,3])
		tab[i,1] <- c$estimate
		tab[i,2] <- c$p.value
		tab[i,3] <- median(ratio[info[int,2]==1]) - median(ratio[info[int,2]==2])
		tab[i,4] <- wilcox.test(ratio[info[int,2]==1], ratio[info[int,2]==2])$p.value
		tab[i,5] <- length(ratio)
	}
}
tab <- tab[!is.na(tab[,1]),]
plim <- 0.05/(nrow(tab))
tab["MSCLSK",4]*nrow(tab)
#[1] 0.01345125

write.table(tab, "age_sex_tissues.txt", quote=F, sep="\t")

i <- which(names(med2)=="MSCLSK")
int <- intersect(names(med2[[i]]), rownames(info))
ratio <- med2[[i]][int]
names(ratio) <- int
#X11(height=5, width=3.5)
pdf("rplots/sexdiff_muscle.pdf", height=5, width=3.5)
boxplot(list(ratio[info[int,2]==1], ratio[info[int,2]==2]), names=c("Male", "Female"), cex.axis=1.3, cex.lab=1.3, lwd=2, col=c("blue", "red"), ylab="Allelic imbalance in muscle")
dev.off()


#per gene
gtab <- read.table("final_impr_tab.txt", as.is=T)
gspl <- split(gtab[,3], gtab[,1])
data <- read.table("gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.filter5.model15.imp.indxgen.final_impr.txt", as.is=T, header=T, sep="\t")
data <- data[data[,1]=="MSCLSK",]
ensg <- read.table("ensg2nametab.txt", as.is=T, sep="\t")
rownames(ensg) <- ensg[,1]
n <- ensg[data[,2],2]
data <- cbind(data, n)
spl <- split(data, n)
med <- lapply(spl, function(z) {(0.5+abs(0.5-(z[,4]/(z[,4]+z[,5]))))})
for(i in 1:length(med)) {
	names(med[[i]]) <- spl[[i]][,3]
}
l <- unlist(lapply(med, length))
med <- med[l>=10] #at least 10 individuals

count <- lapply(spl, function(z) {z[,4]+z[,5]})
for(i in 1:length(count)) {
	names(count[[i]]) <- spl[[i]][,3]
}
l <- unlist(lapply(count, length))
count <- count[l>=10] #at least 10 individuals

info <- read.table("../../data/gtex/GTEx.SubjectPhenotypesDS.v9.1.SUBJID_GENDER_AGE.txt", header=T, as.is=T, sep="\t")
info <- info[!duplicated(info[,1]),]
rownames(info) <- info[,1]

tab <- matrix(data=NA, nrow=length(med), ncol=7)
rownames(tab) <- names(med)
colnames(tab) <- c("SEX_M-F", "SEX_P", "N", "VAR", "READS_F", "READS_M", "READS_P")
for(i in 1:length(med)) {
	int <- intersect(names(med[[i]]), rownames(info))
	ratio <- med[[i]][int]
	c <- count[[i]][int]
		tab[i,1] <- median(ratio[info[int,2]==1]) - median(ratio[info[int,2]==2])
		tab[i,4] <- var(ratio[info[int,2]==1]) - var(ratio[info[int,2]==2])
		tab[i,2] <- wilcox.test(ratio[info[int,2]==1], ratio[info[int,2]==2])$p.value
		tab[i,3] <- length(ratio)
		tab[i,5] <- median(c[info[int,2]==2])
		tab[i,6] <- median(c[info[int,2]==1])
		tab[i,7] <- wilcox.test(c[info[int,2]==1], c[info[int,2]==2])$p.value
}

write.table(tab, "sex_muscle_pergene.txt", quote=F, sep="\t")

colvec <- rep("blue", nrow(tab))  
colvec[tab[,1]<0] <- "red"

#X11(height=5, width=5)
pdf("rplots/sexdiff_musle_pergene.pdf", height=5, width=5)
plot(tab[,2], tab[,1], xlim=c(-0.1,1.1), ylim=c(-0.1, 0.1), xlab="Male/female p-value", ylab="Median MAE M-F", cex.axis=1.5, cex.lab=1.5, col="white")
text(tab[,2], tab[,1], rownames(tab), col=colvec)
dev.off()

#X11(height=5, width=5)
pdf("rplots/sexdiff_musle_pergene_points.pdf", height=5, width=5)
plot(tab[,2], tab[,1], xlim=c(-0.1,1.1), ylim=c(-0.1, 0.1), xlab="Male/female p-value", ylab="Median MAE M-F", cex.axis=1.5, cex.lab=1.5)
dev.off()

#adjust label positions manually




#######   expression differences male/female : imprinted vs nonimprinted genes    #######


#Diff expr data from GTEx transcriptome paper

cd ~/tuuli_lab/tlappalainen/imprint/analysis/expr_revision_2

tot <- read.table("../imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.passed.genes", as.is=T)
imp <- read.table("../imprint/final_impr_tab.txt", as.is=T)
imp <- imp[grep("GD_", imp[,3], invert=T),]
imp <- imp[grep("GC_", imp[,3], invert=T),]
impg <- unique(imp[,2])
esex <- read.table("LMM_Sex_log_pvalues_adjp.txt", as.is=T, header=T)
rownames(esex) <- unlist(lapply(strsplit(esex[,1], ".", fixed=T), function(x) {x[1]}))
int <- intersect(rownames(esex), tot[,1])
esex <- esex[int,]
int <- intersect(impg, rownames(esex))
esexi <- esex[int,]

#per tissue
tisa <- c("ADPSBQ", "WHLBLD", "ARTTBL", "BRAIN", "HRTLV", "LUNG", "MSCLSK", "NERVET", "SKINS", "THYROID")
tis <- c("Adiposetissue", "Blood", "Bloodvessel", "Caudate", "Heart", "Lung", "Muscle", "Nerve", "Skin", "Thyroid")
imp <- read.table("../imprint/final_impr_tab.txt", as.is=T)
imp <- imp[grep("GD_", imp[,3], invert=T),]
imp <- imp[grep("GC_", imp[,3], invert=T),]
impg <- unique(imp[,2])

tot <- read.table("../imprint/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt.impglr.passed", as.is=T, header=T)
tesexl <- vector("list", length(tisa))
res <- mat.or.vec(length(tisa), 2)
rownames(res) <- tisa
for(i in 1:length(tisa)) {
	t1 <- tisa[i]
	t2 <- tis[i]
	totf <- tot[tot[,1]==t1,]
	totf <- totf[,2:5]
	tesexl[[i]] <- read.table(paste("~/tuuli_lab/tlappalainen/imprint/analysis/expr/", t2, "_DESeq_sex.txt", sep=""), as.is=T, header=T)
	tesex <- tesexl[[i]]
	rownames(tesex) <- unlist(lapply(strsplit(tesex[,1], ".", fixed=T), function(x) {x[1]}))
	int <- intersect(rownames(tesex), totf[,1])
	tesex <- tesex[int,]
	timp <- imp[imp[,3]==t1,]
	int <- intersect(timp[,2], rownames(tesex))
	tesexi <- tesex[int,]
	res[i,1] <- wilcox.test(abs(tesex[,6]), abs(tesexi[,6]))$p.value
	res[i,2] <- wilcox.test((tesex[,7]), (tesexi[,7]))$p.value
}

#a=male, b=female. positive fold change = higher in females
#X11(height=5, width=12)
pdf("rplots/exprchange_tissue_volcano.pdf", height=5, width=12)
par(mfrow=c(2,5))
for(i in 1:length(tisa)) {
	tesex <- tesexl[[i]]
	rownames(tesex) <- unlist(lapply(strsplit(tesex[,1], ".", fixed=T), function(x) {x[1]}))
	int <- intersect(rownames(tesex), totf[,1])
	tesex <- tesex[int,]
	timp <- imp[imp[,3]==t1,]
	int <- intersect(timp[,2], rownames(tesex))
	tesexi <- tesex[int,]
	mi <- (-1)*max(abs(tesex[,6])[abs(tesex[,6])!="Inf"], na.rm=T)
	ma <- max(abs(tesex[,6])[abs(tesex[,6])!="Inf"], na.rm=T)
	plot(tesex[,6], (-1)*log10(tesex[,7]), main=tisa[i], xlab="Log2 fold change", ylab="-log10 p-value", xlim=c(mi, ma), cex.axis=1.3, cex.lab=1.3, col="red")
	points(tesex[,6][tesex[,6]<0], ((-1)*log10(tesex[,7]))[tesex[,6]<0], main=tisa[i], xlab="Log2 fold change", ylab="-log10 p-value", xlim=c(mi, ma), cex.axis=1.3, cex.lab=1.3, col="royalblue2")
	points(tesexi[,6], (-1)*log10(tesexi[,7]), col="black")
	legend("topright", legend=paste("p =", round(res[i,1], 4)), bty="n", cex=1.3)
}	
dev.off()






#######   expression differences male/female : maternally vs paternally imprinted    #######

gen <- read.table("../imprint/final_impr_genelist_ensg2name.txt", as.is=T)
int <- intersect(gen[,1], rownames(tesex))
tesexi <- tesex[int,]
rownames(gen) <- gen[,2]
matpat <- read.table("../imprint/impr_known_matpat.txt", as.is=T, header=T, row.names=2)
matpat <- matpat[matpat[,1]!="IGF2" & matpat[,1]!="GRB10",]
int <- intersect(gen[,1], rownames(matpat))
matpat <- matpat[int,]
mae <- rownames(matpat)[matpat[,3]=="Maternal" & !is.na(matpat[,3])]
pae <- rownames(matpat)[matpat[,3]=="Paternal" & !is.na(matpat[,3])]

###collect fold change of each gene across tissues
fdl <- vector("list", length(unique(imp[,2])))
names(fdl) <- unique(imp[,2])
for(j in 1:length(fdl)) {
	for(i in 1:length(tisa)) {
		tesex <- tesexl[[i]]
		rownames(tesex) <- unlist(lapply(strsplit(tesex[,1], ".", fixed=T), function(x) {x[1]}))
		int <- intersect(rownames(tesex), names(fdl)[j])
		if (length(int)>0) {
			tesexi <- tesex[int,]
			fdl[[j]] <- c(fdl[[j]], tesexi[,6])
		}
	}
}

colvec <- rep("grey", length(fdl))
names(colvec) <- names(fdl)
colvec[mae] <- "red"
colvec[pae] <- "royalblue2"

ensg <- read.table("../imprint/ensg2nametab.txt", as.is=T, sep="\t")
rownames(ensg) <- ensg[,1]
fdls <- fdl[order(colvec)]
colvecs <- colvec[order(colvec)]
n <- ensg[names(fdls),2]
n[n=="INPP5F_V2"] <- "INPP5F"
pdf("rplots/pergene_fc_boxplot.pdf", height=6, width=10)
par(mar=c(9,5,2,2))
boxplot(fdls, las=3, names=n, col=colvecs, ylab="log2 fold change M/F", cex.axis=1.3, cex.names=1.3, cex.lab=1.3)
abline(h=0)
legend("topright", legend=c("Maternal expression", "Paternal expression", "Unknown"), fill=c("red", "royalblue2","grey"), cex=1.3, bty="n")
legend("topleft", legend=c("Expression F>M"), cex=1.3, bty="n")
legend("bottomleft", legend=c("Expression M>F"), cex=1.3, bty="n")
dev.off()


#difference between M>F expression for mat/pat expressed genes
mat <- fdls[colvecs=="red"]
pat <- fdls[colvecs=="royalblue2"]
mat <- mat[unlist(lapply(mat, length))>0]
pat <- pat[unlist(lapply(pat, length))>0]
empm <- length(which(unlist(mat)>0))/length(unlist(mat))
#0.53
empp <- length(which(unlist(pat)>0))/length(unlist(pat))
#0.5
perm <- 100000
ratio_perm <- vector("numeric", perm)
fdls_matpat <- fdls[colvecs!="grey"]
fdls_matpat <- fdls_matpat[unlist(lapply(fdls_matpat, length))>0]
vec <- c(rep(1, length(mat)), rep(0, length(pat)))
for(i in 1:perm) {
	vecs <- sample(vec)
	ratio_perm[i] <- length(which(unlist(fdls_matpat[vecs==1])>0))/length(unlist(fdls_matpat[vecs==1]))
}
length(which(ratio_perm>empm))/perm
#0.265



###########################################################################################################################################
#                                                   GRB10 transcript specificity                                                          #
###########################################################################################################################################



#########  GRB10 transcript expression levels  #########  



awk '$2=="ENSG00000106070.12"' /data/NYGC/Resources/GTEx/GTEx_Analysis_RNA-seq_Flux1.2.3_transcript_rpkm__Pilot_2013_01_31.txt >GTEx_Analysis_RNA-seq_Flux1.2.3_transcript_rpkm__Pilot_2013_01_31.GRB10.txt


data <- read.table("GTEx_Analysis_RNA-seq_Flux1.2.3_transcript_rpkm__Pilot_2013_01_31.GRB10.txt", as.is=T)
header <- read.table("/data/NYGC/Resources/GTEx/GTEx_Analysis_RNA-seq_Flux1.2.3_transcript_rpkm__Pilot_2013_01_31.txt", nrow=1, as.is=T)
colnames(data) <- header
colnames(data) <- unlist(lapply(strsplit(colnames(data), ".", fixed=T), function(x) {x[1]}))

data_annot <- data[,1:4]
data <- data[,5:ncol(data)]
rownames(data) <- data_annot[,1]
samples <- read.table("~/tuuli_lab/tlappalainen/gtex/ug_analysis/qc_notes/GTExDecReleaseFreeze1641.Sampleinfo.StandardNames.Colors.20130815.txt", as.is=T, sep="\t", header=T, row.names=1)
int <- intersect(rownames(samples), colnames(data))
sf <- samples[int,] 
data <- data[,int]
datat <- t(data)

tis <- unique(sf$TISSUE_ABBRV)
tmed <- mat.or.vec(length(tis), nrow(data))
rownames(tmed) <- tis
colnames(tmed) <- rownames(data)
for(i in 1:length(tis)) {
	tmp <- as.data.frame(data[,sf$TISSUE_ABBRV==tis[i]])
	tmed[tis[i],] <- apply(tmp, 1, median)
}
log <- apply(tmed, 2, function(x) {length(which(x>1))>1})
tmed <- tmed[,log]
colnames(tmed) <- unlist(lapply(strsplit(colnames(tmed), ".", fixed=T), function(x) {x[1]}))
tmedxbr <- tmed[grep("BRN", rownames(tmed), invert=T),] 
tmedbr <- tmed[grep("BRN", rownames(tmed)),] 
BRAIN <- apply(tmedbr, 2, median)
tmed <- rbind(BRAIN, tmedxbr)

tmedl <- log2(tmed)
colvec <- unique(sf[,c("TISSUE_ABBRV", "TISSUE_RCOL")])
colvec <- rbind(colvec, c("BRAIN", "yellow2"))
rownames(colvec) <- colvec[,1]

br <- tmed["BRAIN",]/sum(tmed["BRAIN",])
ms <- tmed["MSCLSK",]/sum(tmed["MSCLSK",])
bl <- tmed["WHLBLD",]/sum(tmed["WHLBLD",])
brms <- rbind(br, ms, bl)
brms <- brms[,apply(brms, 2, sum)>0.02]
#X11(height=5, width=5)
pdf("rplots/GRB10_brain_muscle_blood_rel_expr.pdf", height=5, width=5)
par(mar=c(10,5,4,2))
barplot(brms, beside=T, las=3, col=c(colvec[c("BRAIN", "MSCLSK","WHLBLD"),2]), ylab="Relative expression", cex.lab=1.3, cex.axis=1.3, cex.names=1.1, ylim=c(0,0.7))
legend("topright", ncol=3, fill=c(colvec[c("BRAIN", "MSCLSK", "WHLBLD"),2]), legend=c("Brain", "Muscle", "Blood"), bty="n", cex=1.1)
dev.off()







###########################################################################################################################################
#                                                   Methylation data                                                                      #
###########################################################################################################################################



cd /data/research/tuuli_lab/tlappalainen/imprint/data/methylation/gencord

#########  annotations  #########  


R
#gene bodies and promoter regions
data <- read.table("~/tuuli_lab/tlappalainen/resource/gencode.v12.annotation.tab.gene", as.is=T)
data <- data[data[,14]=="protein_coding" | data[,14]=="lincRNA",]
pos <- data[data[,7]=='+',]
neg <- data[data[,7]=='-',]
prom <- rbind(cbind(pos[,1], pos[,4]-1000, pos[,4]+2000, pos[,9], pos[,9], pos[,7], pos[,14]), cbind(neg[,1], neg[,5]-2000, neg[,5]+1000, neg[,9], neg[,9], neg[,7], neg[,14]))
genebody <- rbind(cbind(pos[,1], pos[,4]+2000, pos[,5], pos[,9], pos[,9], pos[,7], pos[,14]), cbind(neg[,1], neg[,4], neg[,5]-2000, neg[,9], neg[,9], neg[,7], neg[,14]))
genebody <- genebody[as.numeric(genebody[,2])<as.numeric(genebody[,3]),]
prom <- prom[as.numeric(prom[,2])<as.numeric(prom[,3]),]
write.table(prom, "gencodeV12_promoter_-1kbto2kb.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(genebody, "gencodeV12_geneBodies_2kbEndGene.bed", quote=F, row.names=F, col.names=F, sep="\t")

#court et al DMRs: Suppl table 1, sheet: intersection 145 VS 356
d <- read.table("court_dmrs.txt", header=T, as.is=T, sep="\t")
d$gene[d$gene=="h19"] <- "H19"
d <- d[d$gene!="",]
ens2g <- read.table("../../../analysis/imprint/ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
gene <- ens2g[d$gene,1]
bed <- cbind(d[,1:3], gene)
write.table(bed, "court_dmrs.bed", row.names=F, col.names=F, quote=F, sep="\t")


#########  extract data  #########  

/data/NYGC/software/bedtools/bedtools-2.17.0/bin/intersectBed -wao -a methylation_sites_FILT4.bed -b gencodeV12_geneBodies_2kbEndGene.bed >methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.bed
/data/NYGC/software/bedtools/bedtools-2.17.0/bin/intersectBed -wao -a methylation_sites_FILT4.bed -b gencodeV12_promoter_-1kbto2kb.bed >methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.bed
/data/NYGC/Software/bedtools/bedtools-2.17.0/bin/intersectBed -wao -a methylation_sites_FILT4.bed -b court_dmrs.bed >methylation_sites_FILT4.court_dmr.bed
awk '$5!="."' methylation_sites_FILT4.court_dmr.bed >methylation_sites_FILT4.court_dmr_filt.bed

echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.bed >Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.bed >Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.bed >Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.bed >Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.bed >Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.bed >Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt methylation_sites_FILT4.court_dmr_filt.bed >Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.court_dmr.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt methylation_sites_FILT4.court_dmr_filt.bed >Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.court_dmr.txt' | qsub
echo 'perl ~/imprinting/utilities/OurfSummaryStatsBySegment.pl Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt methylation_sites_FILT4.court_dmr_filt.bed >Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.court_dmr.txt' | qsub


#for landscape plots: whole gene region +- 10kb
perl ~/Geneva/utilities/GencodeGeneChrTssTes.pl ~/tuuli_lab/data/resources/gencode.v12.annotation.gtf >gencode.v12.tss.tes.txt
awk '{print("chr" $2 "\t" $3-10000 "\t" $4+10000 "\t" $1)}' gencode.v12.tss.tes.txt | awk '$2>0 && $3>0' >gencode.v12.tss.tes.10kb.bed
intersectBed -wao -a methylation_sites_FILT4.bed -b gencode.v12.tss.tes.10kb.bed >methylation_sites_FILT4.gencode.v12.tss.tes.10kb.bed
awk '$8!="."' methylation_sites_FILT4.gencode.v12.tss.tes.10kb.bed >methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.bed
cut -f 8 methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.bed | cut -f1 -d"." >tmp
paste methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.bed tmp >methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.bed
perl ~/utilities/intersect.pl ../../../analysis/imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt 1 methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.bed 9 >methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.final_impr.bed
echo 'perl ~/utilities/intersect.pl methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.final_impr.bed 3 Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt 0 >Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.final_impr' | qsub
echo 'perl ~/utilities/intersect.pl methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.final_impr.bed 3 Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt 0 >Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.final_impr' | qsub
echo 'perl ~/utilities/intersect.pl methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.final_impr.bed 3 Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt 0 >Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.final_impr' | qsub






#####   statistics of genes, promoters, DMRs    ######

cd /data/research/tuuli_lab/tlappalainen/imprint/analysis/methylation


metfiles <- c("Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.court_dmr.txt", 
"Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.court_dmr.txt", 
"Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.court_dmr.txt")
impname <- c("GC_FIBRBLS", "GC_LCL", "GC_TCELL") 
name <- c("F-DMR", "L-DMR", "T-DMR")
ens2g <- read.table("../imprint/ensg2nametab.txt", as.is=T)
ens2g <- ens2g[!duplicated(ens2g[,2]),]
rownames(ens2g) <- ens2g[,2]
ens2g2 <- ens2g
rownames(ens2g2) <- ens2g2[,1]

imp <- read.table("../imprint/biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
gslist <- vector("list", 3)
for(f in 1:length(metfiles)) {
    met <- read.table(paste("../../data/methylation/gencord/", metfiles[f], sep=""), as.is=T, header=T, sep="\t", row.names=NULL, fill=T)
	met <- met[!is.na(met[,1]),]
	rownames(met) <- met[,1]
	met <- met[,-1]
    rownames(met) <- unlist(lapply(strsplit(rownames(met), ".", fixed=T), function(x) {x[1]}))
    #int <- intersect(rownames(met), rownames(imp))
	

    spl <- apply(met, 1, function(x) {strsplit(x, ";")})
    metl <- vector("list", 6)
    for(i in 1:6) {
        metl[[i]] <- mat.or.vec(nrow(met), ncol(met))
        rownames(metl[[i]]) <- rownames(met)
        colnames(metl[[i]]) <- colnames(met)
    }
    names(metl) <- c("TOT_SITES", "SEMI_PROP", "MEAN", "MEDIAN", "VAR", "SEMI_PROP_XIND")
    for(i in 1:length(spl)) {
        metl[[1]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[1])}))
        metl[[2]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[2])/as.numeric(x[1])}))
        metl[[3]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[3])}))
        metl[[4]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[4])}))
        metl[[5]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[5])}))
        metl[[6]][i,] <- sum(unlist(lapply(spl[[i]], function(x) {as.numeric(x[2])})))/sum(unlist(lapply(spl[[i]], function(x) {as.numeric(x[1])})))
    }

    gs <- vector("list", 8)
    for(i in 1:length(metl)) {
        gs[[i]] <- apply(metl[[i]], 1, median)
    }
    gs[[7]] <- apply(metl[[2]], 1, function(x) {length(which(x>=0.5))/length(x)})
    gs[[8]] <- apply(metl[[3]], 1, var)
    names(gs) <- c("Meth sites", "Median ind semiprop", "Median ind mean", "Median ind median", "Median ind var", "Tot semiprop", "Prop semimeth ind", "Var ind median")
	gslist[[f]] <- gs

	tab <- do.call(cbind, gs)

	imp <- read.table("../imprint/biallelic.weak.final.tab.txt", as.is=T, row.names=1, header=T)
	imps <- imp[,impname[f]]
    names(imps) <- ens2g[rownames(imp),1]
	imps <- imps[rownames(tab)]
	impvec <- imps
	impvec <- replace(impvec, is.na(impvec), "NA")
	impvec <- replace(impvec, impvec=="+", "BAE")
	impvec <- replace(impvec, impvec=="-", "MAE")
	IMPR <- impvec
	tab2 <- cbind(tab, IMPR)	
	tab2f <- tab2[tab2[,9]!="NA",]
	GENE <- ens2g2[rownames(tab2f),2]
	write.table(cbind(GENE, tab2f), paste("DMR_", name[f], "_stats.txt", sep=""), sep="\t", quote=F)
}	



metfiles <- vector("list", 3)
metfiles[[1]] <- c("Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.txt", 
    "Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.txt")     
metfiles[[2]] <- c("Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.txt", 
    "Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.txt")
metfiles[[3]] <- c("Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_geneBodies_2kbEndGene.txt", 
    "Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.methylation_sites_FILT4-gencodeV12_promoter_-1kbto2kb.txt")
gen <- read.table("../imprint/final_impr_genelist_ensg2name.txt", as.is=T)
name <- c("F", "L", "T")

for(f in 1:length(metfiles)) {
    met <- read.table(paste("../../data/methylation/gencord/", metfiles[[f]][1], sep=""), as.is=T, header=T, sep="\t", row.names=1, fill=T)
    rownames(met) <- unlist(lapply(strsplit(rownames(met), ".", fixed=T), function(x) {x[1]}))
    int <- intersect(rownames(met), gen[,1])
    sd <- setdiff(rownames(met), gen[,1])
    meti <- met[int,]
    metni <- met[sd,]

    spl <- apply(met, 1, function(x) {strsplit(x, ";")})
    metl <- vector("list", 6)
    for(i in 1:6) {
        metl[[i]] <- mat.or.vec(nrow(met), ncol(met))
        rownames(metl[[i]]) <- rownames(met)
        colnames(metl[[i]]) <- colnames(met)
    }
    names(metl) <- c("TOT_SITES", "SEMI_PROP", "MEAN", "MEDIAN", "VAR", "SEMI_PROP_XIND")
    for(i in 1:length(spl)) {
        metl[[1]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[1])}))
        metl[[2]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[2])/as.numeric(x[1])}))
        metl[[3]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[3])}))
        metl[[4]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[4])}))
        metl[[5]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[5])}))
        metl[[6]][i,] <- sum(unlist(lapply(spl[[i]], function(x) {as.numeric(x[2])})))/sum(unlist(lapply(spl[[i]], function(x) {as.numeric(x[1])})))
    }

    gs <- vector("list", 8)
    for(i in 1:length(metl)) {
        gs[[i]] <- apply(metl[[i]], 1, median)
    }
    gs[[7]] <- apply(metl[[2]], 1, function(x) {length(which(x>=0.5))/length(x)})
    gs[[8]] <- apply(metl[[3]], 1, var)
    names(gs) <- c("Meth sites", "Median ind semiprop", "Median ind mean", "Median ind median", "Median ind var", "Tot semiprop", "Prop semimeth ind", "Var ind median")
	gs1 <- gs
	int1 <- int
	sd1 <- sd

    met <- read.table(paste("../../data/methylation/gencord/", metfiles[[f]][2], sep=""), as.is=T, header=T, sep="\t", row.names=1, fill=T)
    rownames(met) <- unlist(lapply(strsplit(rownames(met), ".", fixed=T), function(x) {x[1]}))
    int <- intersect(rownames(met), gen[,1])
    sd <- setdiff(rownames(met), gen[,1])
    meti <- met[int,]
    metni <- met[sd,]

    spl <- apply(met, 1, function(x) {strsplit(x, ";")})
    metl <- vector("list", 6)
    for(i in 1:6) {
        metl[[i]] <- mat.or.vec(nrow(met), ncol(met))
        rownames(metl[[i]]) <- rownames(met)
        colnames(metl[[i]]) <- colnames(met)
    }
    names(metl) <- c("TOT_SITES", "SEMI_PROP", "MEAN", "MEDIAN", "VAR", "SEMI_PROP_XIND")
    for(i in 1:length(spl)) {
        metl[[1]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[1])}))
        metl[[2]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[2])/as.numeric(x[1])}))
        metl[[3]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[3])}))
        metl[[4]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[4])}))
        metl[[5]][i,] <- unlist(lapply(spl[[i]], function(x) {as.numeric(x[5])}))
        metl[[6]][i,] <- sum(unlist(lapply(spl[[i]], function(x) {as.numeric(x[2])})))/sum(unlist(lapply(spl[[i]], function(x) {as.numeric(x[1])})))
    }

    gs <- vector("list", 8)
    for(i in 1:length(metl)) {
        gs[[i]] <- apply(metl[[i]], 1, median)
    }
    gs[[7]] <- apply(metl[[2]], 1, function(x) {length(which(x>=0.5))/length(x)})
    gs[[8]] <- apply(metl[[3]], 1, var)
    names(gs) <- c("Meth sites", "Median ind semiprop", "Median ind mean", "Median ind median", "Median ind var", "Tot semiprop", "Prop semimeth ind", "Var ind median")
	gs2 <- gs
	int2 <- int
	sd2 <- sd

	dmr <- read.table(paste("DMR_", name[f], "-DMR_stats.txt", sep=""), as.is=T, sep="\t")
	colvec <- dmr$IMPR
	colvec <- replace(colvec, colvec=="MAE", "red")
	colvec <- replace(colvec, colvec=="BAE", "black")

	ptab <- mat.or.vec(9, 2)
	rownames(ptab) <- c(names(gs), "0.5-Median ind median")
	colnames(ptab) <- c("GeneBody", "Promoter")
#    X11(height=9, width=7)
    pdf(paste("rplots/imp_all_", name[f], "_boxplot.pdf", sep=""), height=9, width=9)
    par(mfrow=c(3,3))
	for(i in 1:8) {
    boxplot(list(gs1[[i]][int1], gs1[[i]][sd1], gs2[[i]][int2], gs2[[i]][sd2]), lwd=2, cex=1.5, ylab=names(gs)[i], col=c("orange", "grey"), ylim=c(quantile(gs[[i]], 0.05),1.1*max(c(quantile(gs[[i]], 0.98), dmr[,i+1]))),  names=c("Gene Impr", "Gene Others", "Prom Impr", "Prom Others"), las=3, cex.axis=1.5, cex.lab=1.5, xlim=c(0,7))
	points(x=rep(6, nrow(dmr)), y=dmr[,i+1], pch=3, lwd=2, col=colvec, cex=1.2)
#	text(x=rep(6, nrow(dmr)), y=dmr[,i+1], labels=dmr[,1], col=colvec, cex=1.2)
	ptab[i,1] <- wilcox.test(gs1[[i]][int1], gs1[[i]][sd1])$p.value
	ptab[i,2] <- wilcox.test(gs2[[i]][int2], gs2[[i]][sd2])$p.value
	}
    boxplot(list(abs(0.5-gs1[[4]])[int1], abs(0.5-gs1[[4]])[sd1], abs(0.5-gs2[[4]])[int2], abs(0.5-gs2[[4]])[sd2]), col=c("orange", "grey"), lwd=2, cex=1.5, 
    	ylab="|0.5-Median ind median|", ylim=c(min(c(abs(0.5-dmr[,5]), quantile(abs(0.5-gs[[4]]), 0.01))),quantile(abs(0.5-gs[[4]]), 0.99)), las=3, 
    		names=c("Gene Impr", "Gene Others", "Prom Impr", "Prom Others"), cex.axis=1.5, cex.lab=1.5, xlim=c(0,7))
#	text(x=rep(6, nrow(dmr)), y=abs(0.5-dmr[,5]), labels=dmr[,1], col=colvec, cex=1.2)
	points(x=rep(6, nrow(dmr)), y=abs(0.5-dmr[,5]), pch=3, lwd=2, col=colvec, cex=1.2)

	ptab[9,1] <- wilcox.test(abs(0.5-gs1[[4]])[int1], abs(0.5-gs1[[4]])[sd1])$p.value
	ptab[9,2] <- wilcox.test(abs(0.5-gs2[[4]])[int2], abs(0.5-gs2[[4]])[sd2])$p.value
	dev.off()

	write.table(ptab, paste("methstats_ptab_", name[f], ".txt", sep=""), sep="\t", quote=F)

}


#####    landscape plots    #####


metfiles <- c("Normalized_Beta_values_Fs_Genchord2_114_Y_float_OURF_filt.final_impr", 
    "Normalized_Beta_values_LCLs_Genchord2_119_Y_float_OURF_filt.final_impr", 
    "Normalized_Beta_values_Ts_Genchord2_69_Y_float_OURF_filt.final_impr")
dmr <- read.table("../../data/methylation/gencord/methylation_sites_FILT4.court_dmr_filt.bed", as.is=T)
gene <- read.table("../../data/methylation/gencord/gencode.v12.tss.tes.txt", as.is=T)
rownames(gene) <- unlist(lapply(strsplit(gene[,1], ".", fixed=T), function(x) {x[1]}))
name <- c("F", "L", "T")
bed <- read.table("../../data/methylation/gencord/methylation_sites_FILT4.gencode.v12.tss.tes.10kb.ensg.nov.final_impr.bed", as.is=T)
spl <- split(bed, bed[,10])
example <- read.table("../imprint/gencord/GENCORD.536.2MIMP.ASE.COV16.ANNOTPLUS.SINFO.txt.impglr.final_impr.txt", as.is=T, header=T)
genames <- read.table("../imprint/ensg2nametab.txt", as.is=T)
rownames(genames) <- genames[,1]

geneinfo <- read.table("~/tuuli_lab/tlappalainen/resource/gencode.v12.annotation.tab.gene", as.is=T, sep="\t")
rownames(geneinfo) <- unlist(lapply(strsplit(geneinfo[,9], ".", fixed=T), function(x) {x[1]}))
geneinfo <- geneinfo[example[,2],]

for(f in 1:length(metfiles)) {
    met <- read.table(paste("../../data/methylation/gencord/", metfiles[f], sep=""), as.is=T, header=F, sep="\t", row.names=1, fill=T)
    for(i in 1:length(spl)) {
        m <- met[spl[[i]][,4],]
        gst <- gene[names(spl)[i],3]
        gen <- gene[names(spl)[i],4]
        if (nrow(m)>0) {
            annot <- m[,1:3]
            m <- m[,4:ncol(m)]
            s <- order(annot[,3], decreasing=F)
            annot <- annot[s,]
            m <- m[s,]
			log <- !is.na(m[,1])
			m <- m[log,]
			annot <- annot[log,]
            l <- as.list(as.data.frame(t(m)))
			dmrf <- dmr[match(names(l), dmr[,4], nomatch=F),]
            colv <- match(names(l), dmr[,4], nomatch=F)
            colvec <- replace(colv, colv==0, "grey")
            colvec <- replace(colvec, colv!=0, "cyan3")
            colvec2 <- rep("black", length(l))
            colvec2 <- replace(colvec2, (annot[,3]>gst & annot[,3]<gen), "red2")
			if (any(rownames(genames)==names(spl)[i]) & !is.na(geneinfo[names(spl[i]),1])) {
			fiven <- lapply(l, function(x) {c(quantile(x, p=0.1), quantile(x, p=0.5), quantile(x, p=0.9))})

#            X11(height=5, width=min(c(15, max(c(length(l)/3, 3)))))
            pdf(paste("rplots/methylation_landscape_lineplot_", genames[names(spl[i]),2], "_", name[f], ".pdf", sep=""), height=5, width=min(c(15, max(c(length(l)/3, 3)))))
			plot(x=annot[,3], y=unlist(lapply(fiven, function(x) {x[2]})), ylab="Methylation", xlab="Methylation site", ylim=c(0,1), cex=0.5, main=paste(genames[names(spl[i]),2], name[f]), cex.axis=1.3, cex.lab=1.3, type="l")    
			lines(x=annot[,3], y=unlist(lapply(fiven, function(x) {x[1]})), col="grey80")    
			lines(x=annot[,3], y=unlist(lapply(fiven, function(x) {x[3]})), col="grey80")    
			if(geneinfo[names(spl[i]), 7]=="+") {			
				arrows(x1=geneinfo[names(spl[i]), 4], y1=0, x0=geneinfo[names(spl[i]), 5], y0=0, lwd=2)
			}
			if(geneinfo[names(spl[i]), 7]=="-") {			
				arrows(x1=geneinfo[names(spl[i]), 5], y1=0, x0=geneinfo[names(spl[i]), 4], y0=0, lwd=2)
			}
			for(j in 1:nrow(dmrf)) {
				lines(x=dmrf[j,6:7], y=c(1,1), col="cyan3", lwd=3)
			}
            abline(h=0.5, col="blue", lwd=2)
            dev.off()
			            
#            X11(height=5, width=min(c(15, max(c(length(l)/3, 3)))))
            pdf(paste("rplots/methylation_landscape_boxplot_", genames[names(spl[i]),2], "_", name[f], ".pdf", sep=""), height=5, width=min(c(15, max(c(length(l)/3, 3)))))
			genen <- paste(genames[names(spl[i]),2], name[f], "strand", geneinfo[names(spl[i]), 7])
            boxplot(l, names=rep("", length(l)), ylab="Methylation", border=colvec2, col=colvec, xlab="Methylation site", ylim=c(0,1), cex=0.5, main=genen, cex.axis=1.3, cex.lab=1.3)    
            abline(h=0.5, col="blue", lwd=2)

            dev.off()
			}
        }
    }
}






###########################################################################################################################################
#                                                   mmPCR-seq validation                                                                  #
###########################################################################################################################################




cd ~/tuuli_lab/tlappalainen/imprint/analysis/imprint/mmpcr
cp ~/tuuli_lab/tlappalainen/gtex/data/validation/ase_star/GTEX.121.MMPCR.STAR_IMPR.ASE.COV8.ANNOTPLUS.ID.SINFO.txt .
perl ~/Geneva/utilities/intersect2columns_keepheader.pl GTEX.121.MMPCR.STAR_IMPR.ASE.COV8.ANNOTPLUS.ID.SINFO.txt 1 25 ../../../data/gtex/GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.txt 1 0 >GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.MMPCR.txt


mm <- read.table("GTEX.121.MMPCR.STAR_IMPR.ASE.COV8.ANNOTPLUS.ID.SINFO.txt", as.is=T, header=T, sep="\t")
rna <- read.table("GTEX.1582.5MIMP.ASE.COV8.ANNOTPLUSGT.SINFO.BR.MMPCR.txt", as.is=T, header=T, sep="\t")
mmall <- mm
filt <- read.table("../../../data/gtex/imputed_exome_filter_all.txt", as.is=T)

rnaid <- unique(paste(rna[,27], rna[,2], sep="_"))
filtid <- paste(filt[,1], filt[,2], sep="_")
int <- (intersect(rnaid, filtid))
ma <- match(int, rnaid)
rna <- rna[-c(ma),]

mm <- mm[!is.na(mm[,26]),]
rownames(mm) <- paste(mm[,26], mm[,2], sep="_")
rownames(rna) <- paste(rna[,1], rna[,2], sep="_")
int <- intersect(rownames(mm), rownames(rna))
mm <- mm[int,]
rna <- rna[int,]

ma <- max(c(mmall[,10], rna[,10]))
#X11(height=3, width=9)
pdf("rplots/total_count_hist.pdf",height=3, width=9)
par(mfrow=c(1,3))
hist(rna[,10], breaks=20, xlab="Coverage", main="RNA-seq", col="grey", xlim=c(8,ma))
hist(mmall[,10], breaks=40, xlab="Coverage", main="mmPCR all sites", col="grey", xlim=c(8,ma))
hist(mm[,10], breaks=20, xlab="Coverage", main="mmPCR overlapping RNA-seq", col="grey", xlim=c(8,ma))
dev.off()


lim <- c(8, 20, 30, 50)
#X11()
pdf("rplots/ratio_scatters.pdf")
par(mfrow=c(2,2))
for(i in 1:length(lim)) {
	c <- round(cor.test((rna[,8]/rna[,10])[rna[,10]>=lim[i] & mm[,10]>=lim[i]], (mm[,8]/mm[,10])[rna[,10]>=lim[i] & mm[,10]>=lim[i]], methdos="s")$estimate, 2)
	plot((rna[,8]/rna[,10])[rna[,10]>=lim[i] & mm[,10]>=lim[i]], (mm[,8]/mm[,10])[rna[,10]>=lim[i] & mm[,10]>=lim[i]], xlab="Ref/Total RNA-seq", ylab="Ref/Total mmPCR", main=paste("Total counts >=", lim[i]), cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
	legend("topleft", paste("rho =", c), bty="n", cex=1.5, text.col="blue")
}
dev.off()

imp <- read.table("../final_impr_genelist_ensg2name.txt", as.is=T, sep="\t")
ensg <- read.table("../ensg2nametab.txt", as.is=T, sep="\t", row.names=1)

spl <- split(mm, mm[,22])
gene <- names(spl)
genen <- ensg[gene,1]
imprinted <- vector("character", length(genen))
for(i in 1:length(genen)) {
	if (any(genen[i]==imp[,2])) {
		imprinted[i] <- 1
	} else {
		imprinted[i] <- 0
	}		
}
snp_sample_total <- unlist(lapply(spl, function(x) {nrow(x)}))
ind <- unlist(lapply(spl, function(x) {length(unique(x[,28]))}))
snp <- unlist(lapply(spl, function(x) {length(unique(x[,2]))}))
tis <- unlist(lapply(spl, function(x) {length(unique(x[,29]))}))
cov_median <- unlist(lapply(spl, function(x) {median(x[,10])}))
stats <- cbind(gene, genen, imprinted, snp_sample_total, ind, snp, tis, cov_median)
write.table(stats, "mmPCR_star_asedata.txt", row.names=F, quote=F, sep="\t")

median(mm[,10])
median(rna[,10])

splmm <- split(mm, mm[,22])
splrna <- split(rna, rna[,22])
splmm2 <- lapply(splmm, function(x) {split(x, x[,29])})
splrna2 <- lapply(splrna, function(x) {split(x, x[,28])})
ensg <- read.table("../../imprint/ensg2nametab.txt", as.is=T, sep="\t")
rownames(ensg) <- ensg[,1]

for(i in 1:length(splmm2)) {
	pdf(paste("rplots/scatter_snp_", ensg[names(splmm2)[i],1], "_mmpcr_star.pdf", sep=""), width=7, height=3.5)
	par(mfrow=c(1,2))
	n <- names(splmm2[[i]])
	for(j in 1:length(splmm2[[i]])) {
		if (nrow(splrna2[[i]][[n[j]]])>0) {
#			X11(width=7, height=3.5)
			ma <- max(c(splmm2[[i]][[j]][,8], splmm2[[i]][[j]][,9], splrna2[[i]][[j]][,8], splrna2[[i]][[j]][,9]))
			plot(splmm2[[i]][[j]][,8], splmm2[[i]][[j]][,9], ylim=c(0,ma), xlim=c(0,ma), xlab="REF READS", ylab="ALT READS", main=paste(ensg[names(splmm2)[i],1], n[j]))
			points(splrna2[[i]][[j]][,8], splrna2[[i]][[j]][,9], col="red")
			ma <- log(max(c(splmm2[[i]][[j]][,8], splmm2[[i]][[j]][,9], splrna2[[i]][[j]][,8], splrna2[[i]][[j]][,9]))+1)
			plot(log(splmm2[[i]][[j]][,8]+1), log(splmm2[[i]][[j]][,9]+1), ylim=c(0,ma), xlim=c(0,ma), xlab="LOG REF READS", ylab="LOG ALT READS", main=paste(ensg[names(splmm2)[i],1], n[j]))
			points(log(splrna2[[i]][[j]][,8]+1), log(splrna2[[i]][[j]][,9]+1), col="red")
		}
	}
	dev.off()
}












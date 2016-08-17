#!/usr/bin/env Rscript
#Tychele N. Turner 
#Laboratory of Aravinda Chakravarti, Ph.D.
#Johns Hopkins University School of Medicine
#Functionality extended by Ryan Morin, Ph. D., Simon Fraser University
#Protein Plotting Script
#Programming Language: R
#Updated 06/15/2013
#modified 08/05/2014 by R. Morin

#Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. 

#NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

#Required files:
##Mutation file: tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE
##Protein architecture file: tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information.
##Post-translational modification file: This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written.

#Usage:
## R --slave --vanilla < plotProtein.R mutationFile proteinArchitectureFile proteinLength nameOfYourQuery tickSize showLabels zoomIn zoomStart zoomEnd
# x
#without zoom
## R --slave --vanilla < plotProtein.R psen1_mutation_file.txt psen1_architecture_file.txt  463 Test 25 no no

#with zoom
## R --slave --vanilla < plotProtein.R psen1_mutation_file.txt psen1_architecture_file.txt  463 Test 25 no yes 25 50

#Arguments:
argv <- function(x){
    args <- commandArgs()
    return(args[x])
}


mutationFile <- argv(6) #This is the mutation file
print(mutationFile)
proteinArchitectureFile <- argv(7) #This is the protein architecture file
print(proteinArchitectureFile)
exonFile<-argv(8)
print(exonFile)
spliceFile<-argv(9)
print(spliceFile)
proteinLength <- argv(10) #Length of the protein isoform your looking at
nameOfYourQuery <- argv(11) #Here you can put whatever name you want to show up in the plot
tickSize <- as.numeric(argv(12)) #Specify the tick spacing for x-axis
showLabels <- argv(13) #yes/no
outFile<-argv(14)
zoomIn <- argv(15) #yes/no
if(zoomIn == "yes"){
	zoomStart <- as.numeric(argv(16))
	zoomEnd <- as.numeric(argv(17))
}

####################ANALYSIS####################
#Read in the files
var <- read.table(mutationFile, sep="\t")
splice <- try(read.table(spliceFile,sep="\t",header=1))
if(class(splice)=='try-error') {
splice = matrix(nrow=0,ncol=5,dimnames=list(c(),c("ensg","gene_symbol","position","colour","mutation_type")))
}
pa <- read.table(proteinArchitectureFile, sep="\t",header=TRUE,as.is=TRUE)
exon<-read.table(exonFile,sep="\t")

############PLOTTING#############
#x is the input data, y is rpt, z is rpa from HPRD
pdf(outFile, height=7.5, width=10)
#par(oma=c(8, 1.2, 8, 1.2))
layout(matrix(c(1,2),nrow=1), widths=c(1,3))
par(oma=c(4, 0, 4, 0), mar=c(5, 0, 4, 0) + 0.4)

#stable legend
plot((-30:-15), rep(-1, 16), col="white", type="l", ann=FALSE, bty="n", xaxt="n", yaxt="n", xlim=c(-160, -15), ylim=c(1,-5.5))
	
#query text
text(-100,-2,nameOfYourQuery, col="blue", cex=0.9, font=2)

xlimRegion <- c(0, proteinLength)
	if(zoomIn == "yes") {
          xlimRegion <- c(as.numeric(zoomStart), as.numeric(zoomEnd))
	}




xlimRegion <- c(0, as.numeric(proteinLength))
	if(zoomIn == "yes") {
          xlimRegion <- c(as.numeric(zoomStart), as.numeric(zoomEnd))
	}
	
plot((1:as.numeric(proteinLength)), rep(-2, as.numeric(proteinLength)), type="l", lwd=5, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", ylim=c(-1,-4), cex.lab=0.9, cex.main=1, yaxt="n", xlim=xlimRegion, xaxt="n", ann=FALSE, bty="n")

#Plot mutations
for(i in 1:nrow(var)){
	lines(c(var[i,3],var[i,3]),c(-2.2,-2))
}
symbols = c()
symbols['non-synonymous']<-19
symbols['synonymous']<-20
symbols['nonsense']<-18
symbols['indel']<-17

points(var[,3], rep(-2.2, length(var[,3])), pch=symbols[var[,7]], col=var[,6], cex=0.7)

if(showLabels == "yes"){
	#Label mutations
	for(i in 1:nrow(var)){
		text(var[i,3], rep(-2.4, length(var[i,3])), paste(as.character(var[i,4]), as.character(var[i,3]), as.character(var[i,5]), sep=""), col="blue", cex=0.7, srt=90, adj = 0)
        #add lines to figure here?
	}
}

ticks=seq(0,as.numeric(proteinLength), by=tickSize) 
axis(side = 1, at = ticks, las=3)

#labels
legend_vals<-c()
legend_cols<-c()
print(pa$colour)
print(pa$architecture_name)
for(i in 1:length(pa$start_site)){
        rect(as.numeric(pa$start_site[i]), -2.05, as.numeric(pa$end_site[i]), -1.95, col=pa$colour[i])
	legend_cols[pa$architecture_name[i]] = pa$colour[i]
	legend_vals<-c(legend_vals,pa$architecture_name[i])
}
#print(pa$architecture_name[i])
#print(pa$architecture_name[1])
#for(i in 1:length(pa$architecture_name)){
#        text(median(c(as.numeric(pa$start_site[i]), as.numeric(pa$end_site[i]))), -1.80, pa$architecture_name[i], cex=1)
#}
print(exon[1,])
for( i in 1:length(exon[,1])){
	print(paste(exon[i,1],exon[i,2]))
	rect(as.numeric(exon[i,1]), -1.75, as.numeric(exon[i,2]), -1.86)
}

#now plot any splice site mutations
points(splice[,3], rep(-1.7, length(splice[,3])), pch="*", col="blue", cex=1)

#updated this so it uses a different colour for each domain type!
#legend("topright", c("Protein Domain", "Post-Translational Modification"), fill=c("lightseagreen", "deeppink"),  box.col="white", bg="white", cex=1)
print(legend_vals)
print(legend_cols)
legend("topleft",legend=names(legend_cols),fill=legend_cols,box.col='white',bg='white',cex=0.6)
dev.off()



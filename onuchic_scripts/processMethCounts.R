#!/usr/bin/env RScript

# required library
library (GenomicRanges)

# arguments for input file, name of desired output file, and optionally a desired FDR value threshold
args = commandArgs(trailingOnly = TRUE)

inputFile = args[1]

outFile = args[2]

fdrCutoff = 0.1
if(length(args)>2){
   fdrCutoff = as.numeric(args[3])
}

chrOrder = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")


# necessary functions that will be called
export_meth_ai_table = function(methTableGRanges,outFileName){
    outDf = as.data.frame(methTableGRanges)
    outDf = outDf[,-c(4,5)]
    outDf$start = outDf$start-1
    colnames(outDf)[1] = "chromosome"
    write.table(outDf,outFileName,sep="\t",row.names = FALSE,quote = FALSE)
}

import_meth_ai_table = function(fileName){
    inputDf = read.table(fileName,header=TRUE,sep="\t")
    badRows = grep("chr[123456789XY][[:digit:]]?",inputDf[,1],invert=TRUE)
    badRows = c(badRows,grep("[[:digit:]]+",inputDf[,2],invert=TRUE))
    badRows = c(badRows,grep("[[:digit:]]+",inputDf[,3],invert=TRUE))
    if(length(badRows)>0){
        warning("One or more rows did not match the standard Chr\tStart\tEnd format")
        inputDf = inputDf[-badRows,]
    }
    order = c();
    for(chr in chrOrder){
        order = c(order,which(inputDf[,1]==chr))
    }
    inputDf = inputDf[order,]
    gRangesInput = GRanges(seqnames = Rle(as.character(inputDf[,1])),
    ranges = IRanges(start = as.numeric(inputDf[,3]),width = rep(1,nrow(inputDf)),names = paste(inputDf[,1],inputDf[,3],sep=":")),
    Allele1.Methylated = inputDf[,4],
    Allele1.Unmethylated = inputDf[,5],
    Allele2.Methylated = inputDf[,6],
    Allele2.Unmethylated = inputDf[,7],
    Number.of.good.reads = inputDf[,8],
    Is.on.heterozygous.CpG = inputDf[,9],
    Allele1 = inputDf[,10],
    Allele2 = inputDf[,11])
    gRangesInput = sort(gRangesInput)
    return(gRangesInput)
}

compute_p_values_methylation = function(aiTable,fdrCutoff=0.1){
    compute_p = function(x){
        p = try(fisher.test(matrix(x,ncol=2))$p.value,silent=TRUE)
        if(is.character(p)){
            warning(p)
            NA
        } else{
            p
        }
    }
    goodSnps = which(seqnames(aiTable)!="chrY")
    counts = as.matrix(mcols(aiTable)[,1:4])
    coverage = rep(NA,nrow(counts))
    coverage[goodSnps] = apply(cbind(aiTable$Allele1.Methylated+aiTable$Allele1.Unmethylated,
    aiTable$Allele2.Methylated+aiTable$Allele2.Unmethylated),
    1,min)[goodSnps]
    P.Values = rep(NA,nrow(counts))
    P.Values[goodSnps] = apply(counts[goodSnps,],1,compute_p)
    nAIs = c()
    for(i in 0:20){
        nAIs = c(nAIs,sum(p.adjust(P.Values[which(coverage>=i)],"fdr")<=fdrCutoff))
    }
    bestMinCount = seq(0,20)[which.max(nAIs)]
    goodSnps = which(coverage>=bestMinCount)
    FDR = rep(NA,nrow(counts))
    FDR[goodSnps] = p.adjust(P.Values[goodSnps],"fdr")
    mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),P.Values,FDR))
    return(aiTable)
}

compute_methylation_differences = function(aiTable,minCounts=1){
    compute_meth_diff = function(x){
        methLevel1 = x[1]/(x[1]+x[2])
        methLevel2 = x[3]/(x[3]+x[4])
        methLevel1-methLevel2
    }
    counts = as.matrix(mcols(aiTable)[,1:4])
    goodSnps = which((counts[,1]+counts[,2]>=minCounts) & (counts[,3]+counts[,4]>=minCounts) & seqnames(aiTable)!="chrY")
    Methylation.Difference = rep(NA,nrow(counts))
    Methylation.Difference[goodSnps] = apply(counts[goodSnps,],1,compute_meth_diff)
    mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),Methylation.Difference))
    return(aiTable)
}

compute_methylation_al1 = function(aiTable,minCounts=1){
    compute_meth = function(x){
        methLevel = (x[1])/(x[1]+x[2])
        methLevel
    }
    counts = as.matrix(mcols(aiTable)[,1:4])
    goodSnps = which((counts[,1]+counts[,2]>=minCounts) & seqnames(aiTable)!="chrY")
    Methylation.Allele1 = rep(NA,nrow(counts))
    Methylation.Allele1[goodSnps] = apply(counts[goodSnps,],1,compute_meth)
    mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),Methylation.Allele1))
    return(aiTable)
}

compute_methylation_al2 = function(aiTable,minCounts=1){
    compute_meth = function(x){
        methLevel = (x[3])/(x[3]+x[4])
        methLevel
    }
    counts = as.matrix(mcols(aiTable)[,1:4])
    goodSnps = which((counts[,3]+counts[,4]>=minCounts) & seqnames(aiTable)!="chrY")
    Methylation.Allele2 = rep(NA,nrow(counts))
    Methylation.Allele2[goodSnps] = apply(counts[goodSnps,],1,compute_meth)
    mcols(aiTable) = as.data.frame(cbind(as.data.frame(mcols(aiTable)),Methylation.Allele2))
    return(aiTable)
}





# loading and processing meth count file
aiTable = import_meth_ai_table(inputFile)

aiTable = compute_p_values_methylation(aiTable = aiTable,fdrCutoff)

aiTable = compute_methylation_al1(aiTable = aiTable)

aiTable = compute_methylation_al2(aiTable = aiTable)

aiTable = compute_methylation_differences(aiTable = aiTable)

export_meth_ai_table(methTableGRanges = aiTable, outFileName = outFile)

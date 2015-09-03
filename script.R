library(CNVPanelizer)

########################################################################################################
message("Load Input parameters")
########################################################################################################
source(file = "properties.R")

# TODO Find better way to log the imported properties
message(paste0("referenceDirectory                     :  ", referenceDirectory))
message(paste0("sampleDirectory                        :  ", sampleDirectory))
message(paste0("bedFilepath                            :  ", bedFilepath))
message(paste0("ampliconColumnNumber                   :  ", ampliconColumnNumber))
message(paste0("removePcrDuplicates  		       :  ", removePcrDuplicates))
message(paste0("outputDirectory      		       :  ", outputDirectory))
message(paste0("numberOfBootstrapReplicates            :  ", numberOfBootstrapReplicates))

########################################################################################################
message("Retrieve the bam filepaths")
########################################################################################################
# Reference data set
referenceFilenames <- list.files(path = referenceDirectory,
                                 pattern = ".bam$",
                                 full.names = TRUE)

# Sample data set
sampleFilenames <- list.files(path = sampleDirectory,
                              pattern = ".bam$",
                              full.names = TRUE)

########################################################################################################
message("Retrieve the bedfile defined regions")
########################################################################################################
#extract the information from a bed file
genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
                                           ampliconColumn = ampliconColumnNumber,
                                           split = "_")
metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
geneNames = metadataFromGenomicRanges["geneNames"][,1]
ampliconNames = metadataFromGenomicRanges["ampliconNames"][,1]

########################################################################################################
message("Retrieve the data and normalize it")
########################################################################################################
#count the reads in the files of interest
#for the reference and for the samples
referenceReadCounts <- ReadCountsFromBam(referenceFilenames,
                                         sampleNames = referenceFilenames,
                                         genomicRangesFromBed,
                                         ampliconNames = ampliconNames,
                                         removeDup = removePcrDuplicates)

sampleReadCounts <- ReadCountsFromBam(sampleFilenames,
                                      sampleNames = sampleFilenames,
                                      genomicRangesFromBed,
                                      ampliconNames = ampliconNames,
                                      removeDup = removePcrDuplicates)

# Normalize references and samples together
normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts,
                                                 ampliconNames = ampliconNames)

########################################################################################################
message("perform the bootstrap based analysis")
########################################################################################################
# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

# Perform the bootstrap based analysis
bootList <- BootList(geneNames,
                     samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = numberOfBootstrapReplicates)

backgroundNoise <- Background(geneNames,
                              samplesNormalizedReadCounts,
                              referenceNormalizedReadCounts,
                              bootList,
                              replicates = numberOfBootstrapReplicates)

reportTables <- ReportTables(geneNames,
                             samplesNormalizedReadCounts,
                             referenceNormalizedReadCounts,
                             bootList,
                             backgroundNoise)

########################################################################################################
message("Plotting bootstrap distributions")
########################################################################################################

PlotBootstrapDistributions(bootList,
                           reportTables,
                           outputFolder = plotsOutputDirectory,
                           save = TRUE)

########################################################################################################
message("Generate compilation of results")
########################################################################################################

# TODO Issue #10 Rename Report tables entries
#names(reportTables) <- as.vector((sapply((sapply(names(bootList), strsplit, split="_IonXpress")), "[[", 1)))
newNames <- c()
for( i in names(bootList)) {
  newName <- unlist(strsplit(i, "-"))
  newNames <- c(newNames, paste(newName[c(1,2)], collapse = "_"))
}
names(reportTables) <- newNames

dir.create(tablesOutputDirectory)
for(report in names(reportTables)) {
  print(paste("Writing report table for: ", report))
  write.table(reportTables[[report]], file = file.path(tablesOutputDirectory, paste0(report,".tsv")) , sep ='\t', row.names = TRUE, col.names = NA)
}

# TODO Issue #9 Reorder all genes in report tables to match the BED file..
geneNamesUnique <- unique(geneNames)
#reportTables[[1]][geneNamesUnique, ]

reportTablesOrderedAsBED <- list()
for(report in names(reportTables)) {
  reportTablesOrderedAsBED[[report]] <- reportTables[[report]][geneNamesUnique, ]
}

EvaluateReportTable <- function(reportTable, sampleName) {
  kk <- apply(reportTable, 1, GeneAberration)
  result <- do.call(rbind, Map(function(x,y) cbind(x), kk))
  rownames(result) <- rownames(reportTable)
  colnames(result) <- sampleName
  return(result)
}

GeneAberration <- function(stats) {
  geneAberration <- "neutral"
  if (stats["Signif."] & stats["AboveNoise"]) {
    if (stats["UBBtpRatio"] < 1) {
      geneAberration <- "deletion"
    } else {
      geneAberration <- "amplification"
    }
  }
  return(geneAberration)
}

# TODO Improve this..
resultsTable <- EvaluateReportTable(reportTablesOrderedAsBED [[1]], names(reportTablesOrderedAsBED)[1])
for (reports in names(reportTablesOrderedAsBED)[-1]) {
  resultsTable <- cbind(resultsTable, EvaluateReportTable(reportTablesOrderedAsBED[[reports]], reports))
}

WriteListToXLSX(reportTables, filepath = file.path(outputDirectory, "reportTables.xlsx"))

save(reportTables, file = file.path(outputDirectory, "reportTables.RData"))

write.table(resultsTable, file = file.path(outputDirectory, paste0("IonTorrentGeneEffects", ".tsv")) , sep ='\t', row.names = TRUE, col.names = NA)

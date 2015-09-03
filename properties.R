# Use full paths only..'~' not working..
dataRoot <- "/data/CNVPanelizer/ION_TORRENT_DATA"
outputDirectory <- "/wip/scripts/CNVPanelExomComparison/CNVPanelizer/output"
referenceDirectory <- file.path(dataRoot, "HIPOReference")
sampleDirectory <- file.path(dataRoot, "HIPO21")
bedFilepath <- file.path(dataRoot, "ACPv2_WITHOUT_CHR.bed")
ampliconColumnNumber <- 6
removePcrDuplicates <- TRUE
numberOfBootstrapReplicates <- 10000
plotsOutputDirectory <- file.path(outputDirectory, "plots")
tablesOutputDirectory <- file.path(outputDirectory, "tables")


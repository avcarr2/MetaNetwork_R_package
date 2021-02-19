#' Write WGCNA Excel Output
#' @param fullDataSet Object including all data for each accession.
#' @param modulesList List, with each element in the list corresponding to a single module.
#' @param filePath The file to ultimately write to
#' @name createResultsWGCNAExcelWorkbook
#' @export
createResultsWGCNAExcelWorkbook <- function(fullDataSet, modulesList,
                                            filePath = "Results/ResultsWGCNA.xlsx"){
  wb = openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb,
                         paste('WGCNA_workbook'))
  openxlsx::writeData(wb, sheet = 1, fullDataSet)
  for(ii in 1:length(modulesList)) {
    openxlsx::addWorksheet(wb, names(modulesList)[ii])
    openxlsx::writeData(wb, names(modulesList)[ii], modulesList[[ii]])
  }

  openxlsx::saveWorkbook(wb, file = filePath, overwrite = TRUE)
}

#' Merges modules broken up by list with the complete dataset sheet.
#' @param selectedDatabase User uploaded database of UniProt Accessions and Genes
#' @param allData The result of the WGCNA
#' @param dataList List containing module and expression data
#' @param ... passes arguments to workbook generation function
#' @return This function combines allData, containing the complete set of module data,
#' with the dataList object, including a list with each index corresponding to a
#' different module.
#' @name mergeAndWriteWGCNAWorkbook
#' @export
mergeAndWriteWGCNAWorkbook <- function(selectedDatabase,
                                       allData,
                                       dataList, ...){
  allData <- tibble::as_tibble(allData)
  dataList <- lapply(dataList, tibble::as_tibble)
  message("Pulling gene symbols from Uniprot Database...")
  userInputDatabaseSelectedColumns <- tibble::tibble(selectedDatabase$Entry,
                                             selectedDatabase$`Gene names`,
                                             selectedDatabase$`Protein names`)
  colnames(userInputDatabaseSelectedColumns) <- c(names(allData)[1],
                                                  "gene name",
                                                  "protein name")
  dat.resMerged <- dplyr::left_join(allData,
                             userInputDatabaseSelectedColumns,
                             by = names(allData)[1])

  list.cluster.datMerged <- list()
  for(i in seq_along(dataList)){
    list.cluster.datMerged[[i]] <- dplyr::left_join(
                                          dataList[[i]],
                                          userInputDatabaseSelectedColumns,
                                          by = names(allData[1])
                                          )
  }
  message("Saving workbook...")
  names(list.cluster.datMerged) <- names(dataList)
  createResultsWGCNAExcelWorkbook(dat.resMerged, list.cluster.datMerged, ...)
  message("WGCNA workbook saved")
}
#' Merges correlatd module eigenproteins
#' @name mergeModuleEigenproteins
#' @param exprData Matrix of expression data
#' @param moduleColors Module membership for each protein
#' @param cutHeight User input value corresponding to the desired merge height.
#' @export
#' @return Merged module membership and hclust objects
mergeModuleEigenproteins <- function(exprData, moduleColors, cutHeight){

  allDataModuleEigenProteins <- WGCNA::moduleEigengenes(expr = exprData,
                                                 colors = moduleColors)
  preMergeModuleEigenProteins <- allDataModuleEigenProteins$eigengenes
  MEDiss <- 1 - WGCNA::cor(preMergeModuleEigenProteins)
  METree <- fastcluster::hclust(as.dist(MEDiss),
                   method = "average")
  mergedClust <- WGCNA::mergeCloseModules(exprData,
                                   moduleColors,
                                   cutHeight = cutHeight,
                                   verbose = 3)
  mergedClust$METree <- METree
  mergedClust$MEDiss <- MEDiss
  return(mergedClust)
}


#' Plots module adjacency as a heatmap.
#' @param moduleEigenproteins Datatable of module eigenproteins and their corresponding values for each sample.
#' @param fileName Desired filename.
#' @param widthInches Set the width of the output. Deafult parameter is
#' 10 inches.
#' @return Creates a plot of the adjacency heatmap, which relates the adjacency between the
#'  module eigenproteins. Hierarchical clustering is performed using 1-adjacency,
#'  where adjacency is calculated by transforming the correlation. Can show
#'  interactions between two modules.
#'  @name plotAdjacencyHeatmap
#'  @export
plotAdjacencyHeatmap <- function(moduleEigenproteins,
                                 fileName = "Results/Eigenprotein_adjacency heatmap.pdf",
                                 widthInches = 10){
  MET <- WGCNA::orderMEs(MEs = moduleEigenproteins)


  pdf(file = fileName, width = widthInches)
  par(cex = 1.0)
  WGCNA::plotEigengeneNetworks(MET, "Eigenprotein Dendrogram",
                        marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)

  WGCNA::plotEigengeneNetworks(MET,
                        "Eigenprotein adjacency heatmap",
                        marDendro = c(3,4,2,2),
                        xLabelsAngle = 90)
  dev.off()
  message("Adjacency heatmap created")
}

# plotDendrograms <- function(x){
#   extracted_dendro_data <- dendro_data(x)
#
#   ggplot()+
#     geom_segment(data = extracted_dendro_data$segments,
#                  aes(x = x, y = y, xend = xend, yend = yend))+
#     geom_text(data = extracted_dendro_data$labels,
#               aes(x = x, y = y, label = label),
#               size = 3, vjust = 0, angle = 90)
#
# }

#' Plots the correlation-based eigenprotein dendrogram
#' @param moduleEigenproteins Datatable of module eigenproteins and their corresponding values for each sample.
#' @param fileName Desired filename.
#' @param widthInches Set the width of the output. Deafult parameter is 10 inches.
#' @return Creates a plot of the eigenproteins based on correlation dissimilarity (1-correlation(module_i, module_j)).
#' @name plotEigenproteinsNetwork
#' @export
plotEigenproteinsNetwork <- function(moduleEigenproteins,
                                     fileName = "Results/DendrogramEigenproteins.pdf",
                                     widthInches = 10){
  pdf(file = fileName, width = widthInches)
  WGCNA::plotEigengeneNetworks(moduleEigenproteins,
                        "EigenproteinNetwork",
                        marHeatmap = c(3,4,2,2),
                        marDendro = c(3,4,2,5),
                        plotDendrograms = TRUE,
                        xLabelsAngle = 90,
                        heatmapColors = WGCNA::blueWhiteRed(50))
  dev.off()
  message("Eigenproteins dendrogram exported")
}

#' Plots samples versus module eigenproteins in heatmap form.
#' @param moduleEigenproteins Datatable of module eigenproteins and their
#' corresponding values for each sample.
#' @param fileName Desired filename.
#' @param widthInches Set the width of the output. Deafult parameter is 10 inches.
#' @return Creates a heatmap that shows samples vs eigenprotein, with coloration
#' dependent on the value of the eigenprotein.
#' @importFrom grDevices pdf
#' @importFrom graphics par
#' @importFrom grDevices dev.off
#' @importFrom stats dist
#' @name plotmoduleEigenproteinsHeatmap
#' @export
plotmoduleEigenproteinsHeatmap <- function(moduleEigenproteins,
                                           fileName = "Results/ModuleEigenproteinsHeatmap.pdf",
                                           widthInches = 15){
  pdf(file = fileName, width = widthInches)
    heatmap3::heatmap3(moduleEigenproteins,
             distfun = function(x) dist(x, method = "euclidean"),
             main = "Module Eigenproteins",
             cexRow = 0.6, cexCol = 0.6)
  dev.off()
  message("Heatmap of eigenproteins successfully created")
}

#' Plots protein dendrogram with color bar identifying module membership.
#' @param proteinDendrogram Object of type hclust containing the clustering results.
#' @param mergedColors Vector of the module membership of each protein in the
#' proteinDendrogram.
#' @return Creates a hierarchical clustering dendrogram with the module
#' membership of each protein beneath.
#' #' @importFrom grDevices pdf
#' @importFrom graphics par
#' @importFrom grDevices dev.off
#' @name plotProteinDendrogram
#' @export
plotProteinDendrogram <- function(proteinDendrogram,
                                  mergedColors){

  pdf("Results/DendroColorMergedClust.pdf")
    WGCNA::plotDendroAndColors(proteinDendrogram,
                               mergedColors,
                               "Module membership",
                               rowText = mergedColors,
                               dendroLabels = FALSE,
                               hang = 0.03,
                               addGuide = TRUE,
                               guideHang = 0.05,
                               rowTextAlignment = "left",
                               addTextGuide = TRUE)
  dev.off()
  message("Protein dendrogram generated")
}

#' Plots correlation-based sample clustering.
#' @param dendro Dendrogram of samples based on dissimlarity correlation.
#' @param fontSize Default = 15.
#' @param fileName Preferred file name.
#' @return Creates a plot of the hierarchical clustering of the samples,
#' which can identify sample outliers.
#' @name plotSampleClusteringDenro
#' @export
plotSampleClusteringDendro <- function(dendro, fontSize = 15,
                                       fileName = "SampleClustering.pdf"){
  ggdendro::ggdendrogram(dendro, rotate = TRUE) +
    ggplot2::ggtitle("Sample Clustering to Detect Outliers") +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Height") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = fontSize)

    ggplot2::ggsave(filename = fileName,
                    path = "Results",
                    plot = ggplot2::last_plot())
  message("Sample dendrogram created")
}

#' Plots scale free topology diagnostic plots.
#' @param fileName preferred file name.
#' @param powersOutput Result of WGCNA::picksfotthreshhold.
#' @param powers Vector of powers from 1:user input value
#' @param cex1 Resizes the final graphics device
#' @param scaleFreeThreshold User input value corresponding to the R^2 topology threshold
#' @param widthInches Default size of the graphical output
#' @return Creates the diagnostic plot showing the satisfaction of the scale-free
#' topology and the effect of power on the mean connectivity of the proteins.
#' #' @importFrom grDevices pdf
#' @importFrom graphics par
#' @importFrom grDevices dev.off
#' @name plotScaleFreeTopology
#' @export
plotScaleFreeTopology <- function(fileName =
                                    "Results/ScaleFreeTopology.pdf",
                                  powersOutput,
                                  powers,
                                  cex1 = 0.9,
                                  scaleFreeThreshold,
                                  widthInches = 10){
  pdf(fileName, width = widthInches)
  par(mfrow = c(1,2))
  # Scale-free topology fit index as a
  # function of the soft-thresholding power
  ## note, plot was moved from the graphics package to the base package in R 4.0,
  ## so all functions reliant on plotting will break if not run in R 4.0.
  plot(powersOutput$fitIndices[,1],
       -sign(powersOutput$fitIndices[,3])*powersOutput$fitIndices[,2],
       xlab = "Power (beta)",
       ylab = "Scale Free Topology Model Fit,signed R^2",
       type = "n",
       main = paste("Scale independence"))
  graphics::text(powersOutput$fitIndices[,1],
       -sign(powersOutput$fitIndices[,3])*powersOutput$fitIndices[,2],
       labels = powers,
       cex = cex1,
       col = "red")
  # this adds the red line corresponding to R^2. Default = 0.85
  graphics::abline(h = scaleFreeThreshold,
         col = "red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(powersOutput$fitIndices[,1],
       powersOutput$fitIndices[,5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n",
       main = paste("Mean connectivity"))
  graphics::text(powersOutput$fitIndices[,1],
       powersOutput$fitIndices[,5],
       labels = powers,
       cex = cex1,
       col = "red")
  dev.off()
  message("Scale free topology plot successfully created.")
}

#' Plots the topological overlap matrix as a heatmap.
#' @param TOMData TOM object output by the WGCNA workflow
#' @param proteinDendro TOM-based hclust object from workflow
#' @param moduleColors Module colors post-merging
#' @param fileName desired file name
#' @param ... passes arguments to TOMplot function
#' @return Outputs the network heatmap.
#' @name plotTOM
#' @importFrom grDevices hcl.colors
#' @export
plotTOM <- function(TOMData,
                    proteinDendro,
                    moduleColors,
                    fileName = "Results/NetworkHeatmap.png", ...){
  message("Begin generating TOM plot.")
  message("Please be patient. This may take a while")

  grDevices::png(filename = fileName)
  WGCNA::TOMplot(TOMData,
                 proteinDendro,
                 moduleColors,
                 main = "Network heatmap plot, all proteins",
                 col = hcl.colors(n = 55,
                                  palette = "viridis",
                                  alpha = 1,
                                  rev = FALSE), ...) #rev(heat.colors()))
  grDevices::dev.off()
  message("TOM plot successfully generated")
}

#' Plots the topological overlap matrix as a heatmap. Saves as a pdf.
#' @param TOMData TOM object output by the WGCNA workflow
#' @param proteinDendro TOM-based hclust object from workflow
#' @param moduleColors Module colors post-merging
#' @param fileName desired file name
#' @param height height of pdf output
#' @param width width of pdf output
#' @return Outputs the network heatmap as a pdf. Use only for small TOMs
#' @name plotTOMpdf
#' @export
plotTOMpdf <- function(TOMData,
                       proteinDendro,
                       moduleColors,
                       fileName = "Results/NetworkHeatmap.pdf",
                       height=12,
                       width=12){
  message("Use only for small TOMs!")
  message("Begin generating TOM plot.")
  message("Please be patient. This may take a while")

  grDevices::pdf(file = fileName, height, width)
  WGCNA::TOMplot(TOMData,
                 proteinDendro,
                 moduleColors,
                 main = "Network heatmap plot, all proteins",
                 col = hcl.colors(n = 55,
                                  palette = "viridis",
                                  alpha = 1, rev = FALSE),
                 ... = list(labRow = TRUE,
                            colRow = TRUE)) #rev(heat.colors()))
  grDevices::dev.off()
  message("TOM plot successfully generated")
}
#' Tidies the module eigenprotein data for easy access by ggplot2.
#' @param dataInput untidy data as output by the module eigenproteins workflow.
#' @param moduleColors merged module colors
#' @param groupsFile user-uploaded file with experimental groups and sampleIDs
#' @return outputs "tidy" module data for easy plotting using ggplot2.
#' @name tidyModuleDataForOutput
#' @export
tidyModuleDataForOutput <- function(dataInput, moduleColors, groupsFile){

  col.length <- length(colnames(dataInput)) - 1

  dfAddedColumn <- tibble::add_column(dataInput,
                                      "Experiment" = rownames(dataInput))
  dataFile_tidy <- tidyr::gather(dfAddedColumn,
                          key = "Gene",
                          value = "Expression",
                          1:col.length)
  # Create intitial data frame. Using data frame because we get column
  # names that don't require the use of ``, which makes it easier to
  # maintain.
  genes_colors.df <- data.frame(colnames(dataInput), moduleColors)
  colnames(genes_colors.df) <- c("Gene", "moduleColor")
  dataColTidy <- dplyr::left_join(dataFile_tidy, genes_colors.df, "Gene")
  dataColTidy$moduleColor <- as.factor(dataColTidy$moduleColor)
  colnames(groupsFile) <- c("Experiment", "SampleID")
  groupsFile$Experiment <- stringr::str_replace_all(groupsFile$Experiment,
                                           pattern = " ",
                                           replacement = ".")

  allDataFinal <- dplyr::left_join(dataColTidy, groupsFile, "Experiment")
  message("Module data is now tidy")
  return(allDataFinal)
}

#' Writes a .csv file of the variance explained by the module eigenproteins
#' @param datExpr Samples vs proteins data frame
#' @param colors module colors
#' @param MEs module eigenproteins
#' @importFrom utils write.csv
#' @return Write a .csv file that shows how well the module
#' eigenproteins explain the variance in the data.
#' @name writeVarianceExplained
#' @export
writeVarianceExplained <- function(datExpr,
                                   colors,
                                   MEs){
  varianceExplained <- WGCNA::propVarExplained(datExpr,
                                        colors,
                                        MEs,
                                        corFnc = "cor",
                                        corOptions = "use = 'p'")
  write.csv(x = varianceExplained, file = "Results/varianceExplained.csv")
  message("Variance explained file successfully written")
}

#' Plots module eigenproteins by sample
#' @param MEPs Module eigenproteins object with
#' @param module Color-based module identification
#' @param filename The desired filename to save each plot to
#' @return Plots the module eigenproteins by sample
#' @name plotModuleEigenproteinsBySample
#' @export
#' @importFrom ggplot2 ggplot aes geom_col ggtitle theme_bw last_plot ggsave
#' @importFrom magrittr %>%
plotModuleEigenproteinsBySample <- function(MEPs, module, filename){
  selectedEigenproteinModule <- MEPs%>%
    filter(`Module` == module)
  ggplot(data = selectedEigenproteinModule,
         mapping = aes(x = `ModuleEigenprotein`,
                       y = `Experiment`)) +
    geom_col() +
    ggtitle(paste(module)) +
    theme_bw(base_family = "Helvetica")
  ggsave(filename = filename, plot = last_plot(), device = "pdf")
  message("Module eigenproteins exported to Results folder")
}

#' plots the eigenprotein clustering after module merging
#' @name plotEigenproteinClusteringPostMerge
#' @param preMergeDendroData dendrogram of module eigenproteins before module merging
#' @param postMergeDendroData dendrogram of module eigenprotiens after module merging
#' @param cutHeight the user input module merge cut height
#' @param filenName file name for plot saving
#' @param widthInches Default to 10 inches
#' @return Outputs pdf file of the denrograms of the eigenproteins both before
#' and after module merging.
#' @export
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par abline
plotEigenproteinClusteringPostMerge <- function(preMergeDendroData,
                                                postMergeDendroData,
                                                cutHeight,
                                                fileName = "Results/ModuleEigenproteinMergeDendrogram.pdf",
                                                widthInches = 10){

  pdf(file = fileName, width = widthInches)
  par(mfrow = c(2,1))
  par(cex = 0.6)
  plot(preMergeDendroData,
       main = "Clustering of module eigenproteins, pre-merging",
       xlab = "",
       sub = "")
  abline(h = cutHeight, col = "red")
  plot(postMergeDendroData,
       main = "Clustering of module eigenproteins, post-merging",
       xlab = "", sub = "")
  dev.off()
  message("Eigenproteins dendrograms pre- and post-merge cretaed")
}

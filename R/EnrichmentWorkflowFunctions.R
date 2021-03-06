## Allez Enrichment Functions

#' Converts uniprot accessions via biomaRt
#' @param wgcnaResults The complete results from the WGCNA portion of the shiny app
#' @param convertTo What the accessions will be converted to. Defaults to gene symbol
#' @param filterBy What the accession currently are. Default is UniProt accessions
#' @param mart Object of type biomart
#' @return A list of length(wgcnaResults) where each list contains the gene symbols
#'  converted from UniProt accessions
#'  @name convertAccessions
#'  @export
convertAccessions <- function(wgcnaResults,
                              convertTo = "external_gene_name",
                              filterBy = "uniprot_gn_id",
                              mart){
  sheetNumber <- length(wgcnaResults)
  if(sheetNumber == 1){
    UniprotAcessions <- wgcnaResults[,1]
  }else{
    UniprotAcessions <- list()
    for(i in seq_len(sheetNumber)){
      UniprotAcessions[[i]] <- wgcnaResults[[i]][,1]
    }
  }
  ModuleNames <- names(wgcnaResults)
  names(UniprotAcessions) <- ModuleNames

  ConvertedGeneSymbols <- list()
  for(i in seq_len(sheetNumber)){
    ConvertedGeneSymbols[[i]] <- biomaRt::getBM(attributes = convertTo,
                                       filters = filterBy, mart = mart,
                                       values = UniprotAcessions[[i]])$external_gene_name
  }
  return(ConvertedGeneSymbols)
}

## Adapted from the enrich_shiny application found at https://github.com/lengning/Enrich_shiny
#' Helper function to add p-values to the final output of the Allez workflow
#' @param outputAllez The object created at the end of the EnrichAllez function
#' @param Lowersetsize terms with less than this many instances are dropped from the list
#' @param Uppersetsize Terms with greater than this many instances are dropped from the list
#' @param side If True, indicates that the p-values distribution is two-tailed
#' @return Combines the original outputAllez object with the p-values
#' @importFrom magrittr %>%
#' @name addPvaluesToAllezOutput
#' @export

addPvaluesToAllezOutput <- function(outputAllez,
                                    Lowersetsize = 5,
                                    Uppersetsize = 500,
                                    side = "T"){
  if(side=="F")outputAllez$p.value <- stats::pnorm(-abs(outputAllez$z.score))# two tailed
  if(side=="T"){
    prb <- stats::pnorm(outputAllez$z.score)# one tailed
    outputAllez$p.value <- ifelse(1-prb>prb, prb, 1-prb)*2
  }
  outputAllez$p.adj <- stats::p.adjust(outputAllez$p.value, method="BH")
  outputAllez <- outputAllez[which(outputAllez$set.size>Lowersetsize),]
  outputAllez <- outputAllez[which(outputAllez$set.size<Uppersetsize),]
  outputAllezOut <- dplyr::arrange(outputAllez$p.adj)

  message("sets with size < ", Lowersetsize, " or > ", Uppersetsize, " are not considered" )

  message("Successfully added p-values to Allez Output")
  return(outputAllezOut)
}

#' Writes the Allez Enrichment excel worksheet to the results folder
#' @param geneUniverse A vector containing all gene symbols used as the background in enrichment
#' @param AllezPvalues A list containing the calculated p-values from the addPvaluesToAllezOutput function
#' @return Writes the results of the Allez enrichment workbook in .xlsx format
#' @name createAllezEnrichmentXLSX
#' @export
createAllezEnrichmentXLSX <- function(geneUniverse, AllezPvalues){
  wb = openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, paste("GeneUniverse"))
  openxlsx::writeData(wb, sheet = 1, geneUniverse)

  # write modules only tabs
  for(ii in 1:length(AllezPvalues)) {
    openxlsx::addWorksheet(wb, names(AllezPvalues)[[ii]])
    openxlsx::writeData(wb, names(AllezPvalues)[[ii]], AllezPvalues[[ii]][[1]])
  }
  openxlsx::saveWorkbook(wb, file = "Results/AllezEnrichment.xlsx", overwrite=TRUE)
  message("Allez enrichment workbook saved successfully")
}

## Adapted from the enrich_shiny application found at https://github.com/lengning/Enrich_shiny
#' Performs Allez enrichment
#' @param GeneSymbols Gene symbols
#' @param GeneUniverse The total genes identified in this experiment
#' @param SpeciesLibrary The species corresponding to the samples.
#' @param idtype Defaults to "SYMBOL," meaning the standard gene symbol.
#' @param alter Defaults to true. Filters high and low set size values.
#' @param Lowersetsize Default 5. Terms below this count are filtered out .
#' @param Uppsersetsize Default 500. Terms above this count are filtered out.
#' @param outprefix The file written.
#' @param ... Passes arguments
#' @return Returns Allez gene ontology enrichment results.
#' @importFrom BiocGenerics %in%
#' @name enrichAllez
#' @export
enrichAllez <- function(GeneSymbols,
                        GeneUniverse,
                        SpeciesLibrary,
                        idtype = "SYMBOL",
                        alter = TRUE,
                        Lowersetsize = 5,
                        Uppsersetsize = 500,
                        outprefix,
                        ...){
  Scores <- rep(0, length(GeneUniverse))
  names(Scores) <- as.vector(GeneUniverse)
  Scores[base::intersect(names(Scores), GeneSymbols)] <- 1
  allezOutput <- allez::allez(scores = Scores,
                       lib = SpeciesLibrary,
                       idtype = "SYMBOL", library.loc = GeneSymbols)
  allezOutput$setscores <- tibble::rownames_to_column(data.frame(allezOutput$setscores),
                                              var = "GO_Term")
  allezOutput$setscores <- addPvaluesToAllezOutput(outputAllez = allezOutput$setscores)

  if(alter == TRUE){

    CategoryMatrix <- allezOutput$setscores$GO_Term
    AuxMatrix <- allezOutput$aux$set.data
    gInData <- names(which(Scores == 1))
    AuxMatrix <- AuxMatrix[which(AuxMatrix$symbol %in% gInData),]

    MaxInCategory <- pmin(length(CategoryMatrix))

    GenesInCategories <- lapply(1:MaxInCategory, function(x) {
      useg <- AuxMatrix[which(AuxMatrix$go_id == CategoryMatrix[x]), "symbol"]
      numOL <- length(useg)
      gcats <- paste0(useg, collapse = ", ")
      return(list(gcats, numOL))
    })
    NumInCats <- unlist(sapply(GenesInCategories, function(x) x[2]))
    if(length(CategoryMatrix) > 1000){
      NumberInCats <- c(NumInCats, rep(0, length(CategoryMatrix) - MaxInCategory))
    }
    names(NumInCats) <- CategoryMatrix
    GenesInCategories <- unlist(sapply(GenesInCategories, function(x) x[1]))
    if(length(CategoryMatrix > 1000)) {c(GenesInCategories,
                                         rep("",length(CategoryMatrix) - MaxInCategory))}
    names(GenesInCategories) <- CategoryMatrix

    allezOutput$setscores$num.overlap <- NumInCats
    allezOutput$setscores$Genes <- GenesInCategories


  }

  message("Allez enrichment completed")
  if(length(allezOutput$setscores) == 0) message("Warning: output is blank")
  return(allezOutput)
}
#' Gets the correct mart from biomaRt
#' @param dataSetIndex An index input by the shiny reactive values.
#' @return An object of type "mart".
#' @name fetchMart
#' @export
fetchMart <- function(dataSetIndex){
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = dataSetIndex, version = dataSetIndex)

  if(is.null(mart) == FALSE){
    message("Successfully retrieved database")}else{
      stop("Error in database retrieval")
    }
  return(mart)
}

# @description Obsolete function.
#' getGOEnrichment <- function(df, fullGOTermsList, plot = TRUE, file.name){
#'
#'   #get the module color for plotting
#'   modcolor <- names(df)
#'   #
#'   fisher.vec <- vector()
#'   for(i in seq_along(df$`GO:ID`)){
#'     C_a <- df$`GO:ID`[[i]]
#'     index.val <- which(fullGOTermsList$`GO:ID` %in% C_a)
#'     A <- df$n[i]
#'     B <- sum(df$n)
#'     C <- fullGOTermsList$n[index.val] - A
#'     D <- sum(fullGOTermsList$n) - B
#'     fisher.tib <- tibble(rbind(A,C), rbind(B,D))
#'     inter1 <- fisher.test(fisher.tib)
#'     fisher.vec[i] <- inter1$p.value
#'   }
#'
#'
#'   GOterms.pvalue <- df %>%
#'     add_column(fisher.vec)%>%
#'     arrange(fisher.vec)
#'   top.20.pvalues <- GOterms.pvalue[1:20,]
#'
#'   if(plot == TRUE){
#'     p <- ggplot(data = top.20.pvalues, mapping = aes(x = -log10(fisher.vec), y = reorder(GOTermName, -fisher.vec)))+
#'       geom_col()+
#'       geom_vline(xintercept = -log10(0.05/length(df$`GO:ID`)), col = "red")+
#'       scale_y_discrete()+
#'       ggtitle(file.name)+
#'       ylab("Top 20 GO Terms")+
#'       xlab("-log10(p-value)")+
#'       theme_bw()
#'     ggsave(file.name, plot = last_plot())
#'   }
#' }
#'
#  @description Obsolete function.
#' getGOterms <- function(df, GOTerms, plot = TRUE, file.name){
#'   ##get the module color:
#'   mod.color <- df$moduleColors[1]
#'   ## Get the number of unique proteins
#'   protein.count <- length(unique(df$accession))
#'   #get the GO terms from GOTerms, using the df.
#'   cluster.uniprot <- df$accession
#'   index <- which(GOTerms$`Input Accession` %in% cluster.uniprot, useNames = TRUE, arr.ind = TRUE)
#'
#'   GO.listModule <- GOTerms$`Input Accession`[index]
#'   GO.listID <- GOTerms$`GO:ID`[index]
#'   GO.listTerms <- GOTerms$`GO Term Name`[index]
#'   GO.listAspect <- GOTerms$Aspect[index]
#'   GO.df <- tibble(GO.listModule, GO.listID, GO.listTerms, GO.listAspect)
#'   ##Filtering the GO terms for P-type GO terms only
#'   GO.count <- GO.df%>%
#'     group_by(`GO.listID`)%>%
#'     filter(GO.listAspect == "F")%>% ## for both process and function: GO.listAspect == "F"|GO.listAspect == "P
#'     count(`GO.listTerms`)%>%
#'     arrange(-n)
#'   ##Plotting function
#'   if(plot == TRUE){
#'     p <- ggplot(data = GO.count[1:20,], mapping = aes(x = n, y = reorder(GO.listTerms, -n)))+
#'       geom_col()+
#'       ggtitle(paste("GO Terms for", mod.color, "module. ", protein.count, "proteins", sep = " "))+
#'       theme_classic()
#'     ggsave(filename = file.name, last_plot())
#'   }else{return(GO.count)}
#' }

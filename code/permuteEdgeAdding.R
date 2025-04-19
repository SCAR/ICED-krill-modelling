
## Edge adding permutations

## Libraries
library(QPress)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(colorspace)



## Set up --------------------------------------------------------

modVersion <- "v16.2"

allInteractionsExplanTable <- read.csv(
    file.path(getwd(),
              "models",
              paste0("combinedModel_",modVersion,"_pairNo_and_interactionExplanations.csv")))

# rename "Depth of benthos" to just "Depth"
allInteractionsExplanTable$From <- replace(
  allInteractionsExplanTable$From, allInteractionsExplanTable$From== "Depth of benthos", "Depth")
allInteractionsExplanTable$To <- replace(
  allInteractionsExplanTable$To, allInteractionsExplanTable$To== "Depth of benthos", "Depth")

## Factoring information

# Levels of classifications of interactions
classificationLevels <- c("Certain","Uncertain")
classificationLabels <- seq(0,length(classificationLevels)-1)


# Rename climate scenarios and create levels for order of appearance
climateScenarioLevels <- c("A","B","B.2","C","D")
climateScenarioLabels <- c(
  "+Temp & +PP",
  "+Temp & -PP",
  "no sea ice & -PP",
  "-Temp & +PP",
  "-Temp & -PP"
)

# ordering nodes
allNodes <- unique(c(allInteractionsExplanTable$From, allInteractionsExplanTable$To))
responseNodes <- c("Growth","Survival")
regionDriverNodes <- c("Depth","Latitude")
physicalDriverNodes <- c("Autumn primary production", "Temperature","Wind stress")
otherNodes <- allNodes[
  which(!(allNodes %in% c(responseNodes, regionDriverNodes, physicalDriverNodes)))]
otherNodes <- otherNodes[order(otherNodes)]
nodeNames <- c(responseNodes, otherNodes, regionDriverNodes, physicalDriverNodes)


# ordering regions
regionNames <- c(RT1="Open ocean",RT2="Continental shelf",RT3="North Antarctic Peninsula")

# ordering model structure
modelStructureNames <- c("Certain","+ Uncertain")


## Lookup table for press perturbations to make regional type
regionTypeLookupDriver <- data.frame(
  region_id <- names(regionNames),
  region = regionNames,
  Latitude = factor(c("Negative","Positive","Negative")),
  Depth = c("Positive","Negative","Negative"),
  glacialMelt = c(NA, "+ Glacial melt", "+ Glacial melt"),
  hist_model = c("exclGlacialMelt", "all","all"),
  scenarioB2_model = c("no sea ice",NA,"no sea ice")
)
regionTypeLookupDriver$regionLabel <- factor(regionTypeLookupDriver$region,
                                             levels = regionNames,
                                             labels = c("Open ocean", "North\nAntarctic\nPeninsula","Continental\nshelf"))


## Lookup table of which regions and scenarios to not include
regionClimScenarioRm <- data.frame(
  Region = regionNames[c(2,3,3)],
  `Climate scenario` = climateScenarioLevels[3:5]
)
colnames(regionClimScenarioRm)[2] <- c("Climate scenario")


# remove explanations and format for reading into QPress
allEdges <- allInteractionsExplanTable %>%
  mutate(Group = factor(Class.,
                        levels = classificationLevels,
                        labels = classificationLabels)) %>%
  dplyr::select(From, To, Group, Type, Pair) %>%
  mutate(Group = as.numeric(Group)-1,
         From = factor(From, levels = nodeNames),
         To = factor(To, levels=nodeNames),
         Type = factor(Type, levels = c("N","P","U","Z")))


# then get edge labels and order factoring for those
allEdgesLabels <- edge.labels(allEdges) # b4 self limiting
certainEdges <- allEdgesLabels[which(allInteractionsExplanTable$Class. == "Certain")]
certainEdges <- certainEdges[order(certainEdges)]
uncertainEdges <- allEdgesLabels[which(allInteractionsExplanTable$Class. == "Uncertain")]
uncertainEdges <- uncertainEdges[order(uncertainEdges)]


## SIMULATION PARAMETERS --------------------------------

nSims <- 5000

## Regional press perturbations
R1.pp <- c("Latitude"= -1, "Depth"=1)
R2.pp <- c("Latitude" =1, "Depth"=-1)
R3.pp <- c("Latitude"=-1, "Depth"=-1)
regionalPressPerturbations <- list(R1.pp, R2.pp, R3.pp)
names(regionalPressPerturbations) <- regionNames

## Climate change scenario
CCA.pp <- c("Temperature"=1, "Wind stress"=1, "Autumn primary production"=1)
CCB.pp <- c("Temperature"=1, "Wind stress"=1, "Autumn primary production"=-1)
CCC.pp <- c("Temperature"=-1, "Wind stress"=1, "Autumn primary production"=1)
CCD.pp <- c("Temperature"=-1, "Wind stress"=1, "Autumn primary production"=-1)
climatePressPerturbations <- list(CCA.pp, CCB.pp,CCB.pp, CCC.pp, CCD.pp)
names(climatePressPerturbations) <- climateScenarioLevels

source(here::here("R/extend.vector"))
source(here::here("R/impact.barplot.myMod"))

## Certain models

# All certain interactions
certainMod <- retain.groups(allEdges, groups=0) # all certain interactions
certainMod <- enforce.limitation(certainMod) # enforce self-limiting effect
certainModSim <- system.simulate(nSims, certainMod) # simulate

# excluding glacial melt (for open-ocean)
tempEdges <- allEdges
tempEdges$Group[which(tempEdges$From == "Glacial melt")] <- 1
tempEdges$Group[which(tempEdges$To == "Glacial melt")] <- 1
certainMod.noGM <- retain.groups(tempEdges, groups=0)
certainMod.noGM <- enforce.limitation(certainMod.noGM) # enforce self-limiting effect
certainModSim.noGM <- system.simulate(nSims, certainMod.noGM) # simulate

# excluding sea ice (for CC scenario B2 north antarctic peninsula)
tempEdges <- allEdges
tempEdges$Group[grepl("ice", tempEdges$From, ignore.case = TRUE)] <- 1
tempEdges$Group[grepl("ice", tempEdges$To, ignore.case = TRUE)] <- 1
certainMod.noIce <- retain.groups(tempEdges, groups=0)
certainMod.noIce <- enforce.limitation(certainMod.noIce) # enforce self-limiting effect
certainModSim.noIce <- system.simulate(nSims, certainMod.noIce) # simulate

# excluding glacial melt AND sea ice (for CC scenario B2 open ocean)
tempEdges <- allEdges
tempEdges$Group[c(grep("ice", tempEdges$From, ignore.case = TRUE),
                  which(tempEdges$From == "Glacial melt"))] <- 1
tempEdges$Group[c(grep("ice", tempEdges$To, ignore.case = TRUE),
                  which(tempEdges$To == "Glacial melt"))] <- 1
certainMod.noGMIce <- retain.groups(tempEdges, groups=0)
certainMod.noGMIce <- enforce.limitation(certainMod.noGMIce) # enforce self-limiting effect
certainModSim.noGMIce <- system.simulate(nSims, certainMod.noGMIce) # simulate


## add Uncertain interactions models

# All interactions
allMod <- allEdges # all all interactions
allMod <- enforce.limitation(allMod) # enforce self-limiting effect
allModSim <- system.simulate(nSims, allMod) # simulate

# excluding glacial melt (for open-ocean)
tempEdges <- allEdges
tempEdges$Group <- 0
tempEdges$Group[which(tempEdges$From == "Glacial melt")] <- 1
tempEdges$Group[which(tempEdges$To == "Glacial melt")] <- 1
allMod.noGM <- retain.groups(tempEdges, groups=0)
allMod.noGM <- enforce.limitation(allMod.noGM) # enforce self-limiting effect
allModSim.noGM <- system.simulate(nSims, allMod.noGM) # simulate

# excluding sea ice (for CC scenario B2 north antarctic peninsula)
tempEdges <- allEdges
tempEdges$Group <- 0
tempEdges$Group[grepl("ice", tempEdges$From, ignore.case = TRUE)] <- 1
tempEdges$Group[grepl("ice", tempEdges$To, ignore.case = TRUE)] <- 1
allMod.noIce <- retain.groups(tempEdges, groups=0)
allMod.noIce <- enforce.limitation(allMod.noIce) # enforce self-limiting effect
allModSim.noIce <- system.simulate(nSims, allMod.noIce) # simulate

# excluding glacial melt AND sea ice (for CC scenario B2 open ocean)
tempEdges <- allEdges
tempEdges$Group <- 0
tempEdges$Group[c(grep("ice", tempEdges$From, ignore.case = TRUE),
                  which(tempEdges$From == "Glacial melt"))] <- 1
tempEdges$Group[c(grep("ice", tempEdges$To, ignore.case = TRUE),
                  which(tempEdges$To == "Glacial melt"))] <- 1
allMod.noGMIce <- retain.groups(tempEdges, groups=0)
allMod.noGMIce <- enforce.limitation(allMod.noGMIce) # enforce self-limiting effect
allModSim.noGMIce <- system.simulate(nSims, allMod.noGMIce) # simulate



## Join altogether in a list
modelList <- list(
  Certain = list(
    all = certainModSim,
    exclGlacialMelt = certainModSim.noGM,
    exclSeaIce = certainModSim.noIce,
    exclGlacialAndIce = certainModSim.noGMIce
  ),
  All = list(
    all = allModSim,
    exclGlacialMelt = allModSim.noGM,
    exclSeaIce = allModSim.noIce,
    exclGlacialAndIce = allModSim.noGMIce
  )
)


## EDGE ADDING PERMUTATIONS ---------------------------------------------


keepSnow <- TRUE ## keep snow?
if (keepSnow) climatePressPerturbations <- lapply(climatePressPerturbations, function(x) c(x, Snow=1))

set.seed(123)
require(combinat)


## Set up objects for loop

# identify the permuation orders
perms <- permn(uncertainEdges)
nPerm <- length(perms)

# identify regions and scenarios to loop through,
# removing combinations of region x climate scenario that aren't investigated
permCombinations <- expand.grid(
  regionNames, names(climatePressPerturbations)
)
colnames(permCombinations) <- c("Region","Climate scenario")
permCombinations <- anti_join(permCombinations, regionClimScenarioRm)

# identify node outcomes to record
keepNode <- responseNodes



## Loop through all regions and climate change scenarios
for(i in 1:nrow(permCombinations)){

  print(paste(i,"/", nrow(permCombinations)))

  # identify model simulation to use for this combination
  regionMod <- regionTypeLookupDriver$hist_model[
    regionTypeLookupDriver$region == permCombinations$Region[i]]
  noIce <- permCombinations$`Climate scenario`[i] == "B.2"
  noGlacialMelt <- grepl("GlacialMelt", regionMod)
  if(noIce & noGlacialMelt){
    modName <- "exclGlacialAndIce"
  }else if(noIce){
    modName <- "exclSeaIce"
  }else if(noGlacialMelt){
    modName <- "exclGlacialMelt"
  }else{
    modName <- "all"
  }
  modObj = modelList[["All"]][[modName]]

  # identify press perturbations
  # get the climate scenario press perturbation
  climPress <- climatePressPerturbations[[permCombinations$`Climate scenario`[i]]]
  # get the regional driver perturbations
  regionPress <- regionalPressPerturbations[[permCombinations[i,"Region"]]]


  ## Loop through all the different orders of permuting
  for(p in 1:nPerm){

    # Using all the edges in the model, assign different group names to
    # uncertain interactions so they can be added one at a time
    # Make a new edge list with groups that I can subset based on edge classification
    tempEdges <- modObj$edges
    ind <- which(edge.labels(modObj$edges) %in% uncertainEdges)
    tempEdges[ind,"Group"] <- edge.labels(modObj$edges)[ind]

    ## Loop through each step of the permutation order.
    for(s in 0:length(perms[[p]])){ # Start with 0 because 0 == certain model

      # if s==0, it's just the certain model, if not, add the appropriate edge
      if(s==0){
        temp <- retain.groups(tempEdges, 0)
      }else{
        temp <- retain.groups(tempEdges, c(0,perms[[p]][1:s]))
      } # end if step != 0

      # simulate model
      temp$Group <- rep(0, nrow(temp))
      temp <- enforce.limitation(temp)
      tempSim <- system.simulate(nSims, temp)
      temp <- impact.barplot.myMod(
        tempSim, c(climPress, regionPress),
        plot=FALSE, percentage = TRUE)

      # subset temp Results to just the response nodes and format data frame
      temp <- temp[which(rownames(temp) %in% keepNode),]
      temp <- as.data.frame(temp)
      temp$Node <- rownames(temp); rownames(temp) <- NULL
      temp$Region <- permCombinations[i,"Region"]
      temp$Experiment <- paste("Future climate")
      temp$`Climate scenario` <- permCombinations[i,"Climate scenario"]
      temp$permID <- rep(p, nrow(temp)) # ID equivalent to the row of perms
      temp$stepID <- rep(s, nrow(temp)) # ID equivalent to col index of perms
      if(s==0){
        temp$edgeAdded <- "Certain"
      }else{
        temp$edgeAdded = rep(perms[[p]][s], nrow(temp))
      }

      ## rbind dataframe
      if(p==1 & s==0){
        edgePermutations <- temp
      }else{
        edgePermutations <- rbind(edgePermutations, temp)
      }

    } # end of looping through s (edge adding)

  } # end of looping through p (permutations)

  ## Save
  filename <- paste0("edgeAdd_",
                     permCombinations[i,"Region"],"-scenario",
                     permCombinations[i,"Climate scenario"],
                     ".csv")
  filename<- gsub(" ","-", filename)

  dirname <- if (keepSnow) "outputs/future_uncertainEdgeAddingPermutations" else "outputs/future_uncertainEdgeAddingPermutations_noSnowPress"

  write.csv(edgePermutations, here::here(dirName, filename))

} # end of looping through climate and region scenarios


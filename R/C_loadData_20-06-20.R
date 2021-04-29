# read libraries and functions --------------------------------------------
library(ape)
library(phytools)
library(SYNCSA)
library(picante)
library(geiger)
library(ade4)
library(phylobase)
library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(letsR)
source(here::here("R", "functions", "assembly_test_07-8-2019.R"))
source(here::here("R", "functions", "anova.1way.R"))
source(here::here("R", "functions", "GridFilter.R"))


# read data ---------------------------------------------------------------

tree <- ape::read.tree(here::here("data", "tree")) #read tree
W <- as.matrix(read.table(here::here("data", "matrixW.txt"), header= TRUE)) #read community matrix
ancestral.area<- read.table(here::here("data", "matrixEcoNodes.txt"), header= TRUE) #ancestral area data - from 
biogeo<- read.table(here::here("data", "matrixEco.txt"), header= TRUE) #ecoregions of each point in the map
Eco<- read.table(here::here("data", "matrixEco.txt"),header=TRUE) 
rownames(ancestral.area)<- ancestral.area$rows_name #renaming ancestral area data
ancestral.area<- data.frame(matrix(ancestral.area[,2],nrow=length(ancestral.area[,2]),ncol=1, dimnames= list(ancestral.area[,1],"EcoRegion")))
Eco<- data.frame(matrix(Eco[,2], nrow= length(Eco[,2]), ncol= 1, dimnames= list(Eco[,1],"EcoRegion")))
biogeo<- Eco #ecoregions of each point in the map
coords <- read.table("data/coords.txt", sep = ";") #coordinates for all points
temp_trop<- c(rep("temperate", length(1:2248)), rep("tropical", length(2249:nrow(W)))) #categorizing the coordinates


# read results NRI and Ages -----------------------------------------------

ages<- readRDS(here::here("R", "agesResult.rds"))
nri<- readRDS(here::here("R", "nriRes.rds"))

#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
rm(list=ls())
#
# define working directory
wd00 <-"/home/hal9000/Documents/shrfldubuntu18/ggplot_barchart_points"
# I was unable to install 'RcppParallel' as I had 
# ' /etc/profile.d/rstudio.sh' adding  an extra unneeded path to my $PATH
# the $PATH that bash is using.
# On stackoverflow I found this suggestion, to find out which file was
# adding this unneeded extra path
# $ grep --color -H 'PATH=' ~/.bashrc ~/.profile ~/.bash_profile ~/bash.login                      ~/.bash_aliases /etc/bash.bashrc /etc/profile                      /etc/profile.d/* /etc/environment 2> /dev/null
#This showed me the problem was in '/etc/profile.d/rstudio.sh'
# so I edited this file using
# $ sudo nano /etc/profile.d/rstudio.sh
# and just deleted all paths here, so they are not added to my $PATH
# and this allowed Rstudio to find Rscript where it is placed, instead #
# of looking for it in some odd temporary directory that Rstudio once installed
#
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    #
library(tidyr)  
if(!require("phyloseq")){
  install.packages("devtools")
  install.packages("https://github.com/joey711/phyloseq", repos = NULL, type="source")
  BiocManager::install("phyloseq", force = T)
}
# package littler is required for the Rscript
# Rscript is required for the 'RcppCore/RcppParallel'
if(!require("littler")){
  install.packages("littler")
}
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("DESeq2")}

library("DESeq2"); packageVersion("DESeq2")
# get the library for all packages
if (!require("RcppParallel", quietly = TRUE)){
  devtools::install_github("RcppCore/RcppParallel")}
library(RcppParallel)
# check which libraries R draws upon
#.libPaths()
# load packages
library("RcppParallel")
library("ape")
library("maps")
library("phyloseq")
library("polyclip")
library("stringr")
# load packages
library("ggforce")
library("ggimage")
library("ggplot2")
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
library("sf")
#load the library, and check the version
library(BiocManager); packageVersion("BiocManager")
# https://github.com/escamero/mirlyn/
# get the biogeo package 
if(!require(mirlyn)){
  install.packages("BiocManager")
  BiocManager::install("escamero/mirlyn")
}
library(mirlyn)

#require('S4Vectors')

library(ggplot2)
library(phyloseq)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    #
library(tidyr)  

# get libraries for making plots
library(ggplot2)
library(scales)
library(ggrepel)
#________________________________________________________________________
# section 04 - start - import phyloseq data frames
#________________________________________________________________________

psexl06.sam <- read.table(paste0(wd00,"/psexl06.sam.csv"))
psexl06.otu <- read.table(paste0(wd00,"/psexl06.otu.csv"))
psexl06.tax <- read.table(paste0(wd00,"/psexl06.tax.csv"))
# make the tables into phyloseq  data frames
psexl06.tax01 <- phyloseq::tax_table(as.matrix(psexl06.tax))
psexl06.otu <- phyloseq::otu_table(psexl06.otu,taxa_are_rows = T)
psexl06.sam01 <- phyloseq::sample_data(psexl06.sam)
# merge all phyloseq objects. Note, that these needs to be converted into
# phyloseq objects before they can be merged
ps06 <- phyloseq::merge_phyloseq(psexl06.otu, psexl06.tax01,psexl06.sam01)
psexl06 <- ps06
#________________________________________________________________________
# section 04 - end - import phyloseq data frames
#________________________________________________________________________

# get the phyloseq object data frames and get the sampling data frame
# and the taxonomy data frame

psexl06.tax <- as.data.frame(psexl06@tax_table)
psexl06.otu <- as.data.frame(psexl06@otu_table)
psexl06.sam <- as.data.frame(as.matrix(psexl06@sam_data))

# get the row names as this holds the samplNoID - i.e. 'samplNoID'
# which is used in the following steps for combining data frames
# using 'left_join'
psexl06.sam$samplNoID <- row.names(psexl06.sam)
# make the phyloseq-otu-table a data frame, and pivot it longer
# by using the row names (which holds the seq-read-ID - i.e. 'seqNoID')
# as reference - notice that this is also the common column
# for using 'left_join' in the next step
psexl06.otu <- as.data.frame(psexl06@otu_table) %>% 
  tibble::rownames_to_column(var = "seqNoID") %>%
  tidyr::pivot_longer(-seqNoID, names_to = "smplNm",
                      values_to = "seqrd.cnt") %>%
  dplyr::group_by(smplNm)
# replace in the column that has the sample no ID
colnames(psexl06.otu)[grepl("smpl",colnames(psexl06.otu))] <- "samplNoID"
# use dplyr::left_join to match between data frames 
psexl06.tax.otu <- dplyr::left_join(psexl06.otu,
                                    psexl06.tax %>% dplyr::select(
                                      seqNoID,
                                      species,
                                      genus,
                                      family,
                                      order,
                                      class,
                                      phylum,
                                      kingdom
                                    ),
                                    by = "seqNoID")
# define columns to keep
ckeep <- c( "sample_location_no",
            "Sample_number",
            "Location",
            "Latitude",
            "Longitude",
            "shsmpTp")
# ensure the phyloseq matrix is a data frame
psexl06.sam01 <- as.data.frame(as.matrix(psexl06.sam))
# copy a column with the sample type
psexl06.sam01$shsmpTp <-  psexl06.sam01$Sampletype
# use 'global substitute' to remove the part of the string
# that is not needed
psexl06.sam01$shsmpTp <- gsub("Sterivex filter, ","",psexl06.sam01$shsmpTp)
# subset the data frame to only have the columns required for the 
# sample stations
psexl06.sam01<- psexl06.sam01[ckeep]
# copy columns  to have columns that match the sations data frame
psexl06.sam01$id <- psexl06.sam01$sample_location_no
psexl06.sam01$lon <- psexl06.sam01$Longitude
psexl06.sam01$lat <- psexl06.sam01$Latitude
psexl06.sam01$samplNoID <- psexl06.sam01$Sample_number
# use dplyr::left_join to match between data frames 
psexl06.tax.otu01 <- 
  dplyr::left_join(psexl06.tax.otu,
                   psexl06.sam01 %>% dplyr::select(
                     shsmpTp,
                     samplNoID,
                     id,
                     Location
                   ),
                   by = "samplNoID")

#View(psexl06.tax.otu01)
psexl06.tax.otu01$Category <- psexl06.tax.otu01$shsmpTp
psexl06.tax.otu01$Species <- psexl06.tax.otu01$family
psexl06.tax.otu01$n <- psexl06.tax.otu01$seqrd.cnt
#colnames(df_bar_data)
ckeep <- c("id","Species","Category","n" )
df_bar_data2 <- psexl06.tax.otu01[ckeep]
# use the 'substr' function to get only the first letter of the text string
df_bar_data2$Category <- substr(df_bar_data2$Category,1,1)
# define columns to keep
ckeep <- c("id", "lon", "lat")
df_stns2 <- psexl06.sam01[ckeep]
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# - 01 # make 'bar_chart' function - start
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# function to make individual bar charts
bar_chart <- function(selected_id=NA,
                      df,
                      xvar="Category",
                      yvar="n", 
                      group="Species",
                      imgdir="img/",
                      size_x=100, 
                      size_y=200){
  
  df <- df %>%
    filter(id==selected_id)
  
  p <- ggplot(df, aes(x=!!as.name(quo_name(xvar)), 
                      y=!!as.name(quo_name(yvar)), 
                      fill=!!as.name(quo_name(group)))) +
    geom_bar(stat = "identity", position="stack") +
    theme_void(base_size=7) +
    scale_fill_discrete(guide=NULL) +
    theme(axis.text.x = element_text(vjust=3))
  
  return(p)
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# - 01 # make 'bar_chart' function - end
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# remove rows where the sample location is duplicated- duplicates arise
# because each sampling location holds tow ID samples
df_stns2 <- df_stns2[!duplicated(df_stns2),]
# make specific colunms numeric, to ensure they can be plotted later on
library(dplyr)
df_stns2 <- df_stns2 %>% 
  dplyr::mutate_at(vars(lon,lat), as.numeric)
#copy the data frames
df_bar_data <- df_bar_data2
df_stns <- df_stns2
# do the bar charts for each station
id_list <- unique(df_bar_data$id)
plot_list <-  purrr::map(id_list, bar_chart, df=df_bar_data)
# add the image information to the stations dataframe
df_stns$bar_img <- plot_list
# define positions for the images so that they are slightly offset 
# from the station positions
# the value of the offset will depend on the scale of the map , i.e.
# a more zoomed in map will require  a lower offset
df_stns <- df_stns %>%
  dplyr::mutate(x_img = lon+0, y_img=lat+0.16)
# # get a shape file for the Danish EEZ
# eez <- terra::vect("shp/EEZ_polygon.shp")
# terra::crs(eez) <-  terra::crs("EPSG:4326") # add the CRS
# # 
library("ggplot2")
theme_set(theme_bw())
library("sf")
library(ggforce)
library(ggimage)
#install.packages("sf", repos = "http://packages.ropensci.org", type = "source")
#devtools::install_github("r-lib/devtools")
if(!require(sf)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  devtools::install_github("r-spatial/sf")
}
library(sf)
if(!require(rnaturalearthhires)){
  install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source",force=T)
  install.packages("rnaturalearth", repos = "http://packages.ropensci.org", type = "source", force=T)
  install.packages("rnaturalearthdata", repos = "http://packages.ropensci.org", type = "source",force=T)
  devtools::install_github("ropensci/rnaturalearthhires")
}
library("rnaturalearthhires")
library("rnaturalearth")
library("rnaturalearthdata")
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
world <- ne_countries(scale = 10, returnclass = "sf")
library(ggplot2)

# ---- simple plot with map that is limited and points  ----
#make plot
p05 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.4) +
  
  geom_point(data=df_stns,aes(x=lon,y=lat), colour="#FF0000") +
  geom_text(data=df_stns,aes(x=lon,y=lat,label=id),
            hjust=1, vjust=1) +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(10.6, 13.2),
                    ylim = c(54.4, 56.4), 
                    expand = FALSE) 


# ---- make a dataframe to use as legend for the map  ----
df_legend <- df_bar_data %>%
  distinct(Species) %>%
  mutate(n=1)

df_legend <- df_stns %>% 
  merge(df_legend, all=T)


#  ---- plot with map and bars  ----

# the geom_rect plots rectangles with fill according to species
# because we use the same x-values for xmin and xmax and the samew
# y-values for ymin and ymax, they should not be visible
# but they *will* trigger a legend with species colours

p06 <- ggplot(data = world) +
  geom_sf(color = "black", fill = "azure3", lwd=0.4) +
  #define limits of the plot 
  ggplot2::coord_sf(xlim = c(10.6, 13.2),
                    ylim = c(54.4, 56.4), 
                    expand = FALSE) +
  geom_rect(data=df_legend,
            aes(xmin=lon, xmax=lon, ymin=lat, ymax=lat, fill=Species)) +
  geom_point(data=df_stns,aes(x=lon,y=lat)) +
  geom_text(data=df_stns,aes(x=lon,y=lat,label=id),
            hjust=1, vjust=1) +
  ggimage::geom_subview(data=df_stns,aes(x=x_img,y=y_img,subview=bar_img),
                        width=0.1, height=0.36) 
# remove the legend
p06 <- p06 + theme(legend.position="none")

#p06
#set variable to define if figures are to be saved
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p06, 
         filename = paste0(wd00,"/Fig28_v01_stacked_bars_on_map.png"),
         width=210,height=297,
         #width=2*297,height=2*210,
         units="mm",dpi=300)
}
#__________________________________________________________________________
# section 02 - end - make stacked bar plots on map 
#__________________________________________________________________________

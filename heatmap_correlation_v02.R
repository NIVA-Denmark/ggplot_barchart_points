


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

#:
ps06.otu <- as.data.frame(ps06@otu_table)
ps06.tax <- as.data.frame(ps06@tax_table)
# get family per row names, nad order and class
ps06.otu_family <- ps06.tax$family[match(rownames(ps06.otu),ps06.tax$seqNoID)]
ps06.otu_order <- ps06.tax$order[match(rownames(ps06.otu),ps06.tax$seqNoID)]
ps06.otu_class <- ps06.tax$class[match(rownames(ps06.otu),ps06.tax$seqNoID)]
# paste together class order  and family to get a taxonomic group
ps06.otu_txgrp <- paste(ps06.otu_class,
                        ps06.otu_order,
                        ps06.otu_family,
                        sep="_")
# identify not assignable catecories
ps06.otu_txgrp[grepl("NA_NA_NA",ps06.otu_txgrp)] <- "not_assignable"
# add as a column
ps06.otu$txgrp  <- ps06.otu_txgrp
# sum multiple columns per group
# https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
ps06.otu <- ps06.otu %>% 
  group_by(txgrp) %>% 
  summarise(across(everything(), sum))
# make the tibble a data frame
ps06.otu <- as.data.frame(ps06.otu)
# replace the row names
rownames(ps06.otu) <- ps06.otu$txgrp
# remove the row names
ps06.otu$txgrp <- NULL

# get the sample data - this will serve as the 'meta_table'
ps06.sam <- as.data.frame(ps06@sam_data)
# get the sampling time and the freezing time
Tm.st <- ps06.sam$Time
Tm.en <- ps06.sam$Freezing_time
# calculate the difference in seconds between the sampling  and the
# freezing of the sample
#https://stackoverflow.com/questions/68825149/calculate-time-difference-between-2-timestamps-in-hours-using-r
library(lubridate)
differ_hrs.in.sec<- as.numeric(difftime(Tm.en, Tm.st, 
                                        units = "secs"))
out <- seconds_to_period(differ_hrs.in.sec)
differ_HHMM <- sprintf('%02d:%02d', out@hour, out$minute)
# add back difference in seconds onto the sample data frame
ps06.sam$dif.sec <- differ_hrs.in.sec
# and calculate the difference in hours
ps06.sam$dif.hrs<- ps06.sam$dif.sec/(60*60)
# 
ps06.sam$smplTp <- gsub("Sterivex filter, ","",ps06.sam$Sampletype)
ps06.sam$smplLcNo_Type <- paste0(ps06.sam$sample_location_no,"_",ps06.sam$smplTp)
# replace the column names with location numbers plus 
# the sampling type
colnames(ps06.otu) <- ps06.sam$smplLcNo_Type[match(colnames(ps06.otu),
                                                   ps06.sam$Sample_number)]



#________________________________________________________________________
# section 01 -  start - Try correlation map
#________________________________________________________________________
# Trying to use the part with:
# Correlation.R
# ============================================================
# Tutorial on drawing a correlation map using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
# Presented Here

#  https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html


# see this question on advice for getting a file from the web
# https://stackoverflow.com/questions/6299220/access-a-url-and-read-data-with-r
# I went to the example web page : https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
# and right-clicked the files, to be allowed to 'copy the link'
# I could then paste this path for each file in an object
# put the path for the file on the web 
SPE_pitlatrine <- "https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological/SPE_pitlatrine.csv"
# and then read in the file
abund_table<-read.csv(SPE_pitlatrine,row.names=1,check.names=FALSE)
abund_table2 <- ps06.otu
# also get the path for the 'ENV_pitlatrine' file
ENV_pitlatrine <- "https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological/ENV_pitlatrine.csv"
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
abund_table2<-t(abund_table2)


meta_table<-read.csv(ENV_pitlatrine,row.names=1,check.names=FALSE)
meta_table2 <- ps06.sam
# ensure all columns have numeric values

rownames(meta_table2) <- ps06.sam$smplLcNo_Type[match(rownames(meta_table2),
                                                      ps06.sam$Sample_number)]
#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table<-meta_table[rownames(abund_table),]
meta_table2<-meta_table2[rownames(abund_table2),]

#Filter out samples with fewer counts
#abund_table<-abund_table[rowSums(abund_table)>200,]

#Extract the corresponding meta_table for the samples in abund_table
meta_table2<-meta_table2[rownames(abund_table2),]

#You can use sel_env to specify the variables you want to use and sel_env_label to specify the labes for the pannel
sel_env <- c( #"Date", 
  #"Time",
  #"Sampletype",
  "Water_volume_filtered_.mL.",
  "pH",
  "Temperature_.ºC.", 
  #"Freezing_time",
  "Qubit_tube_concentration_with_2ul_sample_.ng.ml.", 
  "Dilution",
  "Sample_no",
  "long",
  "lat",
  "DVFI", 
  #"DV",
  "dif.hrs")
#sel_env<-c("pH","Temp","TS","VS","VFA","CODt","CODs","perCODsbyt","NH4","Prot","Carbo")
sel_env_label <- list(
  'pH'="PH",
  'Temp'="Temperature",
  'TS'="TS",
  'VS'="VS",
  'VFA'="VFA",
  'CODt'="CODt",
  'CODs'="CODs",
  'perCODsbyt'="%CODs/t",
  'NH4'="NH4",
  'Prot'="Protein",
  'Carbo'="Carbon"
)
#
sel_env_label <- list(
  #'Date'="Date", 
  #'Time'="Time",
  #'smplTp'="Sampletype",
  'FvolmL'="Water_volume_filtered_.mL.",
  'pH'="pH",
  'Temp'="Temperature_.ºC.",
  #'FreezT'="Freezing_time",
  'QuConcngmL'= "Qubit_tube_concentration_with_2ul_sample_.ng.ml.",
  'Dilution'="Dilution",
  'Sampleno'="Sample_no",
  'long'="long",
  'lat'="lat",
  'DVFI'= "DVFI",
  #'DV'="DV",
  'difhrs'='dif.hrs')
#
sel_env_label<-t(as.data.frame(sel_env_label))
sel_env_label<-as.data.frame(sel_env_label)
colnames(sel_env_label)<-c("Trans")
sel_env_label$Trans<-as.character(sel_env_label$Trans)

#Now get a filtered table based on sel_env
meta_table_filtered<-meta_table2[,sel_env]
abund_table_filtered<-abund_table2[rownames(meta_table_filtered),]

# replace 'too low' with zero
meta_table_filtered$Qubit_tube_concentration_with_2ul_sample_.ng.ml.[grep("low",meta_table_filtered$Qubit_tube_concentration_with_2ul_sample_.ng.ml.)] <- 0
# ensure all columns have numeric values
meta_table_filtered[] <- lapply(meta_table_filtered, as.numeric)
# replace NAs with zero for pH - not really apropriate, but ...
meta_table_filtered$pH[is.na(meta_table_filtered$pH)] <- 0
# ensure all columns have numeric values
meta_table_filtered[] <- lapply(meta_table_filtered, as.numeric)

#Apply normalisation (either use relative or log-relative transformation)
#x<-abund_table_filtered/rowSums(abund_table_filtered)
x<-log((abund_table_filtered+1)/(rowSums(abund_table_filtered)+
                                   dim(abund_table_filtered)[2]))

x<-x[,order(colSums(x),decreasing=TRUE)]
#Extract list of top N Taxa
N<-51
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)
x<-data.frame(x[,colnames(x) %in% taxa_list])
y<-meta_table_filtered

# ensure all columns have numeric values
#abund_table2[] <- lapply(abund_table2, as.numeric)

abund_table <- abund_table2
#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),
                          t(as.data.frame(strsplit(rownames(abund_table),
                                                   "_"))))
# > head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# T_2_2   T  2  2
# T_2_3   T  2  3
# T_2_6   T  2  6

#Let us group on countries
groups<-grouping_info[,1]

#You can use kendall, spearman, or pearson below:
method<-"kendall"


#Now calculate the correlation between individual Taxa and the environmental data
df<-NULL
for(i in colnames(x)){
  for(j in colnames(y)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],
                     use="everything",method=method),
             cor.test(a[complete.cases(b),],b[complete.cases(b),],
                      method=method)$p.value,k)
      if(is.null(df)){
        df<-tmp  
      }
      else{
        df<-rbind(df,tmp)
      }    
    }
  }
}

df<-data.frame(row.names=NULL,df)
colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
df$Pvalue<-as.numeric(as.character(df$Pvalue))
df$AdjPvalue<-rep(0,dim(df)[1])
df$Correlation<-as.numeric(as.character(df$Correlation))

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-5

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Env)){
    for(j in unique(df$Type)){
      sel<-df$Env==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Env)){
    sel<-df$Env==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for signifant values
df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

#We want to reorganize the Env data based on they appear
# df$Env<-factor(df$Env,as.character(df$Env))

#We use the function to change the labels for facet_grid in ggplot2
Env_labeller <- function(variable,value){
  return(sel_env_label[as.character(value),"Trans"])
}

p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x",labeller=Env_labeller)
#pdf(paste("Correlation_",adjustment_label[adjustment],".pdf",sep=""),height=8,width=22)
print(p)


#________________________________________________________________________
# section 01 -  end - Try correlation map
#________________________________________________________________________

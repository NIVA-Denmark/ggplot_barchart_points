


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
# read in otu and sam data frames from csv files
ps06.otu.gh <- "https://github.com/NIVA-Denmark/ggplot_barchart_points/raw/main/psexl06.otu.csv"
ps06.sam.gh <- "https://github.com/NIVA-Denmark/ggplot_barchart_points/raw/main/psexl06.sam.csv"
ps06.tax.gh <- "https://github.com/NIVA-Denmark/ggplot_barchart_points/raw/main/psexl06.tax.csv"
# read in csv files from github
abund_table2<-read.csv(ps06.otu.gh,row.names=1,check.names=FALSE, sep=" ")
tax_table2<-read.csv(ps06.tax.gh,row.names=1,check.names=FALSE, sep=" ")
# make a column that has phylum, class, order, family, genus and species
tax_table2$pcofgs  <- paste(tax_table2$phylum,
                            tax_table2$class,
                            tax_table2$order,
                            tax_table2$family,
                            tax_table2$genus,
                            tax_table2$species,
                            sep="_")
# get a new column with name comprising 'phylum, class, order, family, genus and species'
abund_table2$pcofgs <- tax_table2$pcofgs[match(row.names(abund_table2),tax_table2$qseqid)]
# sum across groups - see this question: https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
abund_table2 <- abund_table2 %>% 
  group_by(pcofgs) %>% 
  summarise(across(everything(), sum))
# make the tibble a data frame
abund_table2 <- as.data.frame(abund_table2)
# then replace the row names with name comprising 'phylum, class, order, family, genus and species'
row.names(abund_table2) <- abund_table2$pcofgs
# and get rid of the column that holds the name comprising 'phylum, class, order, family, genus and species'
abund_table2$pcofgs <- NULL
# also get the path for the 'ENV_pitlatrine' file
ENV_pitlatrine <- "https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological/ENV_pitlatrine.csv"
# read in csv file from github page
meta_table<-read.csv(ENV_pitlatrine,row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
abund_table2<-t(abund_table2)
# read in the sample data frame as the meta data
meta_table2<-read.csv(ps06.sam.gh,row.names=1,check.names=FALSE, sep=" ")
# make a column that has both the sample location number and
# the sampling type
meta_table2$smpl.tp.lcno <- paste0(gsub("Sterivex filter, ",
                                        "",meta_table2$Sampletype),
                                   "_",meta_table2$sample_location_no)
# replace row names
row.names(meta_table2) <- meta_table2$smpl.tp.lcno
# arrange row names to have the sample type first then followed
# by sampling location number
row.names(abund_table2) <- meta_table2$smpl.tp.lcno[match(row.names(abund_table2),meta_table2$Sample_number)]
#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table==meta_table[rownames(abund_table),]
meta_table2==meta_table2[rownames(abund_table2),]

#Filter out samples with fewer counts
#abund_table<-abund_table[rowSums(abund_table)>200,]

#Extract the corresponding meta_table for the samples in abund_table
meta_table2<-meta_table2[rownames(abund_table2),]
# replace in column names as these perhaps might cause problems later on
colnames(meta_table2) <- gsub("\\.","_",colnames(meta_table2))
colnames(meta_table2) <- gsub("\\º","",colnames(meta_table2))
# change some of the long odd column names
colnames(meta_table2)[grepl("volume_filt",colnames(meta_table2))] <- 'FvolmL'
colnames(meta_table2)[grepl("Temperature",colnames(meta_table2))] <- 'TmpC'
colnames(meta_table2)[grepl("Qubit_tube",colnames(meta_table2))] <- 'QuCnngmL'

#You can use sel_env to specify the variables you want to use 
# and sel_env_label to specify the labes for the pannel
sel_env2 <- c(
  "FvolmL",
  "pH",
  "TmpC",
  "QuCnngmL",
  "Dilution",
  "Sample_no",
  "long",
  "lat",
  "DVFI"
)
# 
sel_env<-c("pH","Temp","TS","VS","VFA","CODt","CODs","perCODsbyt","NH4","Prot","Carbo")
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
sel_env_label2 <- list(
  "FvolmL"='FvolmL', 
  "pH"='pH',
  "TmpC"='TmpC',
  "QuCnngmL"='QuCnngmL', 
  "Dilution"='Dilution',
  "Sample_no"='Sampleno', 
  "Longitude"='long', 
  "Latitude"='lat',
  "DVFI"='DVFI')
#

#sel_env_label<-t(as.data.frame(sel_env_label))
#sel_env_label2<-t(as.data.frame(sel_env_label2))
sel_env_label<-as.data.frame(sel_env_label)
sel_env_label2<-as.data.frame(sel_env_label2)
colnames(sel_env_label)<-c("Trans")
colnames(sel_env_label2)<-c("Trans")
sel_env_label$Trans<-as.character(sel_env_label$Trans)
sel_env_label2$Trans<-as.character(sel_env_label2$Trans)
# head(sel_env_label,4)
# head(sel_env_label2,4)
#Now get a filtered table based on sel_env
meta_table_filtered<-meta_table[,sel_env]
meta_table_filtered2<-meta_table2[,sel_env2]

abund_table_filtered<-abund_table[rownames(meta_table_filtered),]
abund_table_filtered2<-abund_table2[rownames(meta_table_filtered2),]
# ensure all columns have numeric values
meta_table_filtered[] <- lapply(meta_table_filtered, as.numeric)
meta_table_filtered2[] <- lapply(meta_table_filtered2, as.numeric)
# replace NAs with zero for pH - not really apropriate, but ...
meta_table_filtered$pH[is.na(meta_table_filtered$pH)] <- 0
meta_table_filtered2$pH[is.na(meta_table_filtered2$pH)] <- 0
meta_table_filtered2$QuCnngmL[is.na(meta_table_filtered2$QuCnngmL)] <- 0

#Apply normalisation (either use relative or log-relative transformation)
#x<-abund_table_filtered/rowSums(abund_table_filtered)
x<-log((abund_table_filtered+1)/(rowSums(abund_table_filtered)+
                                   dim(abund_table_filtered)[2]))

x2<-log((abund_table_filtered2+1)/(rowSums(abund_table_filtered2)+
                                     dim(abund_table_filtered2)[2]))

x<-x[,order(colSums(x),decreasing=TRUE)]
x2<-x2[,order(colSums(x2),decreasing=TRUE)]
#Extract list of top N Taxa
N<-51
N2<- ncol(abund_table2)

taxa_list<-colnames(x)[1:N]
taxa_list2<-colnames(x2)[1:N2]

#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
taxa_list2<-taxa_list2[!grepl("Unknown",taxa_list2)]
N<-length(taxa_list)
N2<-length(taxa_list2)
x<-data.frame(x[,colnames(x) %in% taxa_list])
x2<-data.frame(x2[,colnames(x2) %in% taxa_list2])
y<-meta_table_filtered
y2<-meta_table_filtered2

# ensure all columns have numeric values
#abund_table2[] <- lapply(abund_table2, as.numeric)
#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),
                          t(as.data.frame(strsplit(rownames(abund_table),
                                                   "_"))))

#Get grouping information
grouping_info2<-data.frame(row.names=rownames(abund_table2),
                           t(as.data.frame(strsplit(rownames(abund_table2),
                                                    "_"))))
# 
# head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# head(grouping_info2)

#Let us group on countries
groups<-grouping_info[,1]
groups2<-grouping_info2[,1]

x <- x2
y <- y2
groups <- groups2

#You can use kendall, spearman, or pearson below:
method<-"kendall" 
method<-"pearson" 
#method<-"spearman"

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
df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                     label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

#We want to reorganize the Env data based on they appear
# df$Env<-factor(df$Env,as.character(df$Env))

#We use the function to change the labels for facet_grid in ggplot2
Env_labeller <- function(variable,value){
  return(sel_env_label2[as.character(value),"Trans"])
}

p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
#p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x",labeller=Env_labeller)
p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")
#p<-p+facet_grid(Env ~ ., labeller = labeller(day_count=labels_param))

#pdf(paste("Correlation_",adjustment_label[adjustment],".pdf",sep=""),height=8,width=22)
print(p)

#set variable to define if figures are to be saved
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = p, 
         filename = paste0(wd00,"/Fig29_v01_heatmap_correlat.png"),
         width=210*2,height=297*2.4,
         #width=297,height=210,
         units="mm",dpi=300)
}
# end 'bSaveFigures' test


#________________________________________________________________________
# section 01 -  end - Try correlation map
#________________________________________________________________________

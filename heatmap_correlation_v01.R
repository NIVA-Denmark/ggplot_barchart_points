
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________

# Correlation.R
if(!require("ggplot2")){
  install.packages("ggplot2")
}
library(ggplot2)
# ============================================================
# Tutorial on drawing a correlation map using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
# see this question on advice for getting a file from the web
# https://stackoverflow.com/questions/6299220/access-a-url-and-read-data-with-r
# I went to the example web page : https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
# and right-clicked the files, to be allowed to 'copy the link'
# I could then paste this path for each file in an object
# put the path for the file on the web 
SPE_pitlatrine <- "https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological/SPE_pitlatrine.csv"
# and then read in the file
abund_table<-read.csv(SPE_pitlatrine,row.names=1,check.names=FALSE)
# also get the path for the 'ENV_pitlatrine' file
ENV_pitlatrine <- "https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological/ENV_pitlatrine.csv"
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
meta_table<-read.csv(ENV_pitlatrine,row.names=1,check.names=FALSE)


abund_table<-read.csv(SPE_pitlatrine,row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
meta_table<-read.csv(ENV_pitlatrine,row.names=1,check.names=FALSE)

#Filter out samples with fewer counts
abund_table<-abund_table[rowSums(abund_table)>200,]

#Extract the corresponding meta_table for the samples in abund_table
meta_table<-meta_table[rownames(abund_table),]

#You can use sel_env to specify the variables you want to use and sel_env_label to specify the labes for the pannel
sel_env<-c("pH","Temp","TS","VS","VFA","CODt","CODs","perCODsbyt",
           "NH4","Prot","Carbo")
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

sel_env_label<-t(as.data.frame(sel_env_label))
sel_env_label<-as.data.frame(sel_env_label)
colnames(sel_env_label)<-c("Trans")
sel_env_label$Trans<-as.character(sel_env_label$Trans)

#Now get a filtered table based on sel_env
meta_table_filtered<-meta_table[,sel_env]
abund_table_filtered<-abund_table[rownames(meta_table_filtered),]

#Apply normalisation (either use relative or log-relative transformation)
#x<-abund_table_filtered/rowSums(abund_table_filtered)
x<-log((abund_table_filtered+1)/(rowSums(abund_table_filtered)+dim(abund_table_filtered)[2]))

x<-x[,order(colSums(x),decreasing=TRUE)]
#Extract list of top N Taxa
N<-51
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)
x<-data.frame(x[,colnames(x) %in% taxa_list])
y<-meta_table_filtered

#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
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
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
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
#df$Env<-factor(df$Env,as.character(df$Env))
df <- df[order(df$Env),]
#We use the function to change the labels for facet_grid in ggplot2
Env_labeller <- function(variable,value){
  return(sel_env_label[as.character(value),"Trans"])
}

p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
p <- p + geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
p<-p+theme(axis.text.x = element_text(angle = 90, 
                                      hjust = 1, vjust=0.5))
p<-p+geom_text(aes(label=Significance), 
               color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x",
                labeller=Env_labeller)

# see the plot
p


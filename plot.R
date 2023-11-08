library(renv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggimage)
library(terra)
library(tidyterra)

# -------- test data -----------------------

# create test data frame for stations
id <- c(1,2,3,4) # station id
lon <- c( 6.1,  8.7,  9.4, 11.9)
lat <- c(55.9, 55.2, 57.1, 54.9)
df_stns <- data.frame(id=id, lon=lon, lat=lat)

# create test data.frame for bar charts

set.seed(99)

Species <- c("a","b","c","d","e")
Category <- c("x","y","z")

df_bar_data <- merge(data.frame(Species=factor(Species,levels=Species)),
                 data.frame(Category=factor(Category, levels=Category)),
                 all=T)

# we include station id to match the positions
df_bar_data <- merge(data.frame(id=id),
                     df_bar_data,
                 all=T)

# add random values for our y-variable "n"
df_bar_data <- df_bar_data %>%
  mutate(n=rnorm(nrow(df_bar_data), mean=10, sd=4)) %>%
  mutate(n=ifelse(n<0,0,n))


# calculate as fractions of the total within each category/stn
df_bar_data <- df_bar_data %>%
  group_by(id, Category) %>%
  mutate(n=n/sum(n,na.rm=T)) %>%
  ungroup()


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

# do the bar charts for each station
id_list <- unique(df_bar_data$id)
plot_list <-  purrr::map(id_list, bar_chart, df=df_bar_data)

# add the image information to the stations dataframe
df_stns$bar_img <- plot_list

# define positions for the images so that they are slightly offset 
# from the station positions

df_stns <- df_stns %>%
  mutate(x_img = lon+0.3, y_img=lat+0.3)


# get a shape file for the Danish EEZ
eez <- terra::vect("shp/EEZ_polygon.shp")
terra::crs(eez) <-  terra::crs("EPSG:4326") # add the CRS

# ---- simple plot with EEZ and points  ----
ggplot() +
  geom_spatvector(data=eez) +
  geom_point(data=df_stns,aes(x=lon,y=lat), colour="#FF0000") +
  geom_text(data=df_stns,aes(x=lon,y=lat,label=id),
            hjust=1, vjust=1)


# ---- make a dataframe to use as legend for the map  ----
df_legend <- df_bar_data %>%
  distinct(Species) %>%
  mutate(n=1)

df_legend <- df_stns %>% 
  merge(df_legend, all=T)


#  ---- plot with EEZ and bars  ----

# the geom_rect plots rectangles with fill according to species
# because we use the same x-values for xmin and xmax and the samew
# y-values for ymin and ymax, they should not be visible
# but they *will* trigger a legend with species colours

ggplot() +
  geom_rect(data=df_legend,
            aes(xmin=lon, xmax=lon, ymin=lat, ymax=lat, fill=Species)) +
  geom_spatvector(data=eez) +
  geom_point(data=df_stns,aes(x=lon,y=lat)) +
  geom_text(data=df_stns,aes(x=lon,y=lat,label=id),
            hjust=1, vjust=1) +
  geom_subview(data=df_stns,aes(x=x_img,y=y_img,subview=bar_img),
               width=0.5, height=0.5) 


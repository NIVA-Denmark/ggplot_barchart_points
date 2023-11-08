ggplot - plot bars as points
================

A plain R script with the code is here: [plot.R](./plot.R)

# Test data

Create test dataframe for stations

``` r
id <- c(1,2,3,4) # station id
lon <- c( 6.1,  8.7,  9.4, 11.9)
lat <- c(55.9, 55.2, 57.1, 54.9)
df_stns <- data.frame(id=id, lon=lon, lat=lat)
```

Create test dataframe for bar chart data

``` r
set.seed(99)

Species <- c("a","b","c","d","e")
Category <- c("x","y","z")

df_bar_data <- merge(data.frame(Species=factor(Species,levels=Species)),
                 data.frame(Category=factor(Category, levels=Category)),
                 all=T)
```

we include station id with the test data to give a postion to each bar
chart

``` r
df_bar_data <- merge(data.frame(id=id),
                     df_bar_data,
                 all=T)
```

add random values for our y-variable “n”

``` r
df_bar_data <- df_bar_data %>%
  mutate(n=rnorm(nrow(df_bar_data), mean=10, sd=4)) %>%
  mutate(n=ifelse(n<0,0,n))
```

normalise values to as fractions of the total within each
category/station combination

``` r
df_bar_data <- df_bar_data %>%
  group_by(id, Category) %>%
  mutate(n=n/sum(n,na.rm=T)) %>%
  ungroup()
```

# Function to make individual bar charts

``` r
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
```

# Make the bar charts

Make the bar chart for each station and save them locally

``` r
id_list <- unique(df_bar_data$id)
plot_list <-  purrr::map(id_list, bar_chart, df=df_bar_data) 
```

Add the list of image files to the station data frame

``` r
df_stns$bar_img <- plot_list
```

Define positions for the images so that they are slightly offset from
the station positions.

I haven’t done it here but We could also specify separate sizes for each
bar chart.

``` r
df_stns <- df_stns %>%
  mutate(x_img = lon+0.3, y_img=lat+0.3)
```

# Make the plot

Get a shape file for the Danish EEZ

``` r
eez <- terra::vect("shp/EEZ_polygon.shp")
terra::crs(eez) <-  terra::crs("EPSG:4326") # add the CRS
```

### Simple map with EEZ and points

``` r
ggplot() +
  geom_spatvector(data=eez) +
  geom_point(data=df_stns,aes(x=lon,y=lat), colour="#FF0000") +
  geom_text(data=df_stns,aes(x=lon,y=lat,label=id),
            hjust=1, vjust=1)
```

![](README_files/figure-gfm/plot%20stations-1.png)<!-- -->

### Map with bar charts

The bars are introduced to the ggplot as images

We need to find another way to add the legend for the bars to the map.

We can do this by plotting a layer which uses the same fill variable as
the bar charts (Species). We will use a geom_rect layer - this plots
rectangles with fill according to species. Because we use the same
x-values for xmin and xmax and the same y-values for ymin and ymax, then
the rectangles will have zero height and width. Therefore they should
not be visible but they *will* trigger the creation of a legend with
species colours.

``` r
df_legend <- df_bar_data %>%
  distinct(Species) %>%
  mutate(n=1)

df_legend <- df_stns %>% 
  merge(df_legend, all=T)
```

And here is the plot:

``` r
ggplot() +
  geom_rect(data=df_legend,
            aes(xmin=lon, xmax=lon, ymin=lat, ymax=lat, fill=Species)) +
  geom_spatvector(data=eez) +
  geom_point(data=df_stns,aes(x=lon,y=lat)) +
  geom_text(data=df_stns,aes(x=lon,y=lat,label=id),
            hjust=1, vjust=1) +
  geom_subview(data=df_stns,aes(x=x_img,y=y_img,subview=bar_img),
               width=0.5, height=0.5)
```

![](README_files/figure-gfm/final%20plot-1.png)<!-- -->

You will need to mess around with the themes and formatting in the bar
chart plot function and with the size and positioning of the bars.

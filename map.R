## input: directed edge list + node data matrix with virus strain and geographical location
## output: infection map


## input directories
virus_dir = "virus.csv" # virus strain information; 
# columns: 
# long (longitude), 
# lat (latitude), 
# strain (strain name), 
# sig (size of virus impact for plot size), 
# time (tiem of infection?), etc.
vphy_dir = "phylogeny.csv" # virus phylogeny information;
# columns:
# start_long (longitude), 
# start_lat (latitude), 
# end_long (longitude), 
# end_lat (latitude), 


## libraries
install.packages(c("plotly","dplyr"))
library("plotly")
library("dplyr")


# viral strains
virus = read.csv(virus_dir)

# viral phylogeny
vphy = read.csv(vphy_dir)
vphy$id = seq_len(nrow(vphy))

# map attributes
geo = list(
  # scope = 'north america',
  projection = list(type = 'azimuthal equal area'),
  showland = TRUE,
  landcolor = toRGB("gray95"),
  countrycolor = toRGB("gray80")
)

# map plot
p = plot_geo(locationmode = 'USA-states', color = I("red")) %>%
  add_markers(
    data = virus, x = ~long, y = ~lat, text = ~strain,
    size = ~sig, hoverinfo = "text", alpha = 0.5
  ) %>%
  add_segments(
    data = group_by(vphy, id),
    x = ~start_lon, xend = ~end_lon,
    y = ~start_lat, yend = ~end_lat,
    alpha = 0.3, size = I(1), hoverinfo = "none"
  ) %>%
  layout(
    title = 'Influenza Virus Infection Visualization',
    geo = geo, showlegend = FALSE, height=800
  )

# shareable link to map
chart_link = api_create(p, filename="virus")









# Shiny app
library(shiny); library(ggiraph); library(ukcovidtools); library(ggplot2); library(tidyverse)
library(ggspatial); library(maptools); library(sp); library(rgdal); library(rgeos); library(patchwork); library(ggrepel)

load('R0timeseries.RData')
data("UKCovidMaps")
r0shapes = UKCovidMaps$unitaryAuthority %>% 
  left_join(R0timeseries, by=c("ctyua19cd"="GSS_CD")) %>% 
  mutate(ago=difftime(date,lubridate::now(),units="days")) %>% 
  filter(!is.na(date))
r0shapes = r0shapes %>% mutate(`Median(R)` = ifelse(`Median(R)`>10, 9.999,`Median(R)`))

keyDates = tibble(
  date = as.Date(c("2020-03-13","2020-03-16","2020-03-19","2020-03-23")), #max(r0shapes$date-1, na.rm=TRUE)),
  impactDate = as.Date(c("2020-03-14","2020-03-21","2020-03-24","2020-03-28")), #max(r0shapes$date-1, na.rm=TRUE)),
  event = c("Inpatient only testing","Social isolation of vulnerable","Travel ban / school closure","Stay at home") #,"Latest")
) %>% mutate(label = paste0(date,": \n",event))
r0shapes_key = r0shapes %>% inner_join(keyDates, by="date")

# Function for creating the map
createMap <- function(r0shapes, Location, Date){
  
  # Subset data by date
  r0shapes <- subset(r0shapes, date == Date)  
  
  # Map
  mm <- ggplot(r0shapes)+
    geom_sf(aes(fill=`Median(R)`), data=r0shapes)+
    scale_fill_gradient2(
      low="green",
      mid="white",
      high="red",
      midpoint=0,
      trans="log",
      na.value = "grey80", 
      limits=c(0.1,10), 
      breaks=c(0.1,0.4,1,2.5,10), 
      labels=c("<0.1","0.4","1","2.5",">10"))
  
  if (Location == "East of England"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-1, 2), ylim = c(51.25, 53.25), expand = FALSE)
  }
  
  if (Location == "London"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-0.7, 0.5), ylim = c(51.25, 51.75), expand = FALSE)
  }
  
  if (Location == "Midlands"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-3.5, 0.5), ylim = c(51.5, 53.75), expand = FALSE)
  }
  
  if (Location == "North East and Yorkshire"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-2.75, 0.25), ylim = c(53.5, 56), expand = FALSE)
  }
  
  if (Location == "North West"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-3.75, -1.75), ylim = c(52.75, 55.25), expand = FALSE)
  }
  
  if (Location == "South East"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-1.5, 1.5), ylim = c(50.5, 51.75), expand = FALSE)
  }
  
  if (Location == "South West"){
    mm = mm + coord_sf(crs = 4326, xlim = c(-6, -1), ylim = c(49.75, 51.75), expand = FALSE)
  }

  return(mm)
}





# Define the UI
ui = fluidPage(
  
  # App title
  titlePanel("Rt by unitary authority over time"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs 
    sidebarPanel(
      selectInput(inputId = "Location",
                  label = "Select a region:",
                  choices = list("England" = "England", 
                                 "East of England" = "East of England", 
                                 "London" = "London", 
                                 "Midlands" = "Midlands",
                                 "North East and Yorkshire" = "North East and Yorkshire",
                                 "North West" = "North West",
                                 "South East" = "South East",
                                 "South West" = "South West")),
      
      sliderInput(inputId = "Date",
                  label = "Date:",
                  min = min(r0shapes$date),
                  max = max(r0shapes$date),
                  value = max(r0shapes$date),
                  step = 1),
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      plotOutput("map")
    )
  )
)

# Define the server
server = function(input, output) {
  
  # Create the map
  output$map <- renderPlot({createMap(r0shapes, input$Location, input$Date) })
  
  #output$map <- renderGirafe({
  #  ggiraph(code = print(createMap(r0shapes, input$Location, input$Date)), width_svg = 8, height_svg = 8)
  #})

}

# Finally, we can run our app by either clicking "Run App" in the top of our RStudio IDE, or by running
shinyApp(ui = ui, server = server)

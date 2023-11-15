################### Load Libraries ###################

#data wrangling packages
library(data.table)
library(tidyverse)
library(dplyr)
library(lubridate)

#plotting packages
library(ggplot2)
library(plotly)
library(gganimate)

library(png)
library(gifski)

#spatial packages
library(maps)
library(sf)


################### Read Data Sets ###################

#set working directory to current directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#sif data 
sifdatUSA<-fread(file="../data/SIF_USA.csv")
sifdatCA<-fread(file="../data/SIF_CA_updated.csv")

#Mendocino forest boundary shape file
Mboundary<-st_read("../data/Perimeters/Mendocino/FS_National_Forests_Dataset_(US_Forest_Service_Proclaimed_Forests).shp")
#Fire boundary shape file
Bboundary<-st_read("../data/Perimeters/Ranch-fire-perimeter/California_Wildland_Fire_Perimeters_(All).shp")


################### Add time variables ###################

# For full USA 
#sif data frame of USA for manipulation
df<-sifdatUSA
#add month variable
df$month<-month(df$date) 
#add unique pixel identification by lon & lat
df$lonlat<-paste0("x",round(df$lon,2),"y",round(df$lat,2)) #add unique coordinate for pixel location
#add new "normalized" day variable --- days since first observation
time<-as.numeric(df$date) #numeric time 
mintime<-min(time) #April 1 2018
df$dayssinceapril1<-time-mintime #Add new day variable as days since 
dfUSA<-df

#For CA
df<-sifdatCA
#add month variable
df$month<-month(df$date) 
#add unique pixel identification by lon & lat
df$lonlat<-paste0("x",round(df$lon,2),"y",round(df$lat,2)) #add unique coordinate for pixel location
#add new "normalized" day variable --- days since first observation
time<-as.numeric(df$date) #numeric time 
mintime<-min(time) #April 1 2018
df$dayssinceapril1<-time-mintime #Add new day variable as days since 
dfCA<-df


################### Obtain Mean Monthly SIF over full US  ###################

mean_sif <- dfUSA %>% 
  group_by(year, month, lonlat) %>% 
  summarise(sif = mean(sif, na.rm=TRUE), lon=unique(lon), lat=unique(lat), time=unique(floor_date(date, "month")))
mean_sif<-data.frame(mean_sif)

################### Visualize Mean Monthly SIF over the US  ###################

#Obtain coordinates for US map
sf::sf_use_s2(FALSE)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))

#plot base layer US map
state_plot<-ggplot(data = states) +
  geom_sf() +
  geom_sf(data = states, fill = "black", colour="black") + 
  labs(x="longitude", y="latitude")

#obtain unique sif locations from data
loc<-mean_sif %>% group_by(lonlat) %>% summarise(lon=unique(lon), lat=unique(lat))

#convert sif locations into sf object 
sif_loc<-st_as_sf(loc, coords = c("lon", "lat"), crs = 4326)

#check the coordinates reference system for the us map and sif locations
st_crs(sif_loc)==st_crs(states) #different coordinate systems, need to convert to same before joining data sets

#change states crs to sif_loc crs to combine the data sets
states<-st_transform(states,st_crs(sif_loc))

#confirm coordinate systems are the same
st_crs(states)==st_crs(sif_loc) #TRUE

#spatial join 2 data sets. This gives unique pixel locations along with state ID
state_info<-st_join(sif_loc, states)

#covert to data frame 
state_info<-data.frame(state_info)

#data frame with pixel ID and state ID
group_state<-state_info[c(1,3)] 

#new data set with state IDs 
df_new <- merge(mean_sif,group_state,by="lonlat") 

# remove the sif data with "NA" states IDs (water bodies)
df_new2<-df_new[-which(is.na(df_new$ID)),]

#subset to desired year, if any 
df_new3<-subset(df_new2, year==2019)

#create new month year variable for animation titles 
months_names<-c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
time_my<-paste0(months_names[df_new3$month]," ", df_new3$year)
df_new3$time_my<-factor(time_my, levels=c("January 2019", "February 2019", "March 2019", "April 2019", "May 2019", "June 2019", "July 2019", "August 2019", "September 2019", "October 2019", "November 2019", "December 2019"))

# faceted plot 
facet_map<-state_plot+geom_tile(data=df_new3, aes(x = lon, y=lat, fill=sif))+scale_fill_gradient(low = "black", high = "green", na.value="white")+
  facet_wrap(~factor(time_my))

# animated plot
sif_map1 <- state_plot+
  geom_tile(data=df_new3, aes(x = lon, y=lat, fill=sif))+
  scale_fill_gradient(low = "black", high = "green")+
  labs(title = 'Mean SIF: {closest_state}') +
  gganimate::transition_states(time_my, wrap=TRUE) 

#render animation
animate(sif_map1, renderer = gifski_renderer()) #view animation

anim_save("figures/sif_anim.gif", sif_map1,height = 400, width =700) #view animation 


################### Visualize Map of Forest and Burned Regions  ###################

# calif bounding box
rlon <- c(-125,-113) 
rlat <- c(31.5,42.5)

Forest_plot<-state_plot+
  geom_sf(data = states, fill = "white")+
  geom_sf(data = Mboundary, color = NA, fill = "gold", alpha=0.7)+
  coord_sf(xlim = c(rlon[1], rlon[2]-6), ylim = c(rlat[1]+5, rlat[2]), expand = FALSE)+
  theme(panel.background = element_rect(fill = "slategray3",colour = "slategray3"))

Burned_Forest_Plot<-state_plot+
  geom_sf(data = states, fill = "white")+
  geom_sf(data = Mboundary, color = NA, fill = "gold", alpha=0.7)+
  geom_sf(data = Bboundary, color = "black", fill="black", alpha=0.8)+
  coord_sf(xlim = c(rlon[1], rlon[2]-6), ylim = c(rlat[1]+5, rlat[2]), expand = FALSE)+
  theme(panel.background = element_rect(fill = "slategray3",colour = "slategray3"))

################### Determine whether points within the forest are burned or unburned  ###################

#add forest membership to data
m<-st_transform(Mboundary, crs = 4326)
#join forest info and sif location
m_info<-st_join(sif_loc, m)
m_info<-data.frame(m_info)

#sif location and forest membership
group_m<-m_info[c(1,4)]
group_m<-data.frame(group_m)
names(group_m)[2]<-"forest"

#merge california sif info with forest info
df_m <- merge(dfCA,group_m,by="lonlat")

#change forest name to "Mendocino"
df_m$forest[which(!is.na(df_m$forest))]<-"Mendocino"
head(df_m)

#add burn membership to data
b<-st_transform(Bboundary, crs = 4326)
b_info<-st_join(sif_loc, b)
b_info<-data.frame(b_info)

#sif locations and burn membership
group_b<-b_info[c(1,3)]
group_b<-data.frame(group_b)
names(group_b)[2]<-"group"
df_b <- merge(df_m,group_b,by="lonlat")
head(df_b)
df_b$group[which(!is.na(df_b$group))]<-"burned" #replace object id with label "burned" to indicated fire 
head(df_b)

#obtain burn membership locations
df_b_loc = df_b %>% group_by(lonlat) %>% summarise(lon=unique(lon), lat=unique(lat), group=unique(group))

################### Visualize Map of Forest and SIf locations  ###################
Forest_burned_sif_map<-state_plot+
  geom_sf(data = states, fill = "white")+
  geom_sf(data = Mboundary, color = "black", fill = "lightgreen", alpha=0.8)+
  geom_sf(data = Bboundary, color = "black", fill = "firebrick", alpha=0.8)+
  geom_point(data=df_b_loc, aes(x=lon, y=lat, col=group), alpha=0.8)+
  coord_sf(xlim = c(rlon[1], rlon[2]-6), ylim = c(rlat[1]+5, rlat[2]), expand = FALSE)+
  theme(panel.background = element_rect(fill = "azure3",colour = "azure3"))

ggplotly(Forest_burned_sif_map)

#use above interactive plot to set burned/unburned memberships. Cross check with SIF time series 
df_b$group[df_b$lonlat=="x-122.9y39.1"]<-"burned" 

df_b$group[df_b$lonlat=="x-122.9y40.1"]<-"unburned"
df_b$group[df_b$lonlat=="x-122.7y40.1"]<-"unburned"
df_b$group[df_b$lonlat=="x-123.1y39.7"]<-"unburned"


## isolate just the mendocino forest 
df_m2<-subset(df_b, forest=="Mendocino" | group =="burned" | group == "unburned") 
tail(df_m2)
second_fire<-as.Date("2020-08-16") #date before the next fire in the region -- cut off here 
df_m2<-subset(df_m2, date < second_fire)
#set remaining locations in the forest and outside the burn area to "unburned"
df_m2$group[which(is.na(df_m2$group))]<-"unburned"

#obtain final burned/unburned sif locations
df_b_loc_final = df_m2 %>% group_by(lonlat) %>% summarise(lon=unique(lon), lat=unique(lat), group=unique(group))
table(df_b_loc_final$group)

#visualize the map with the sif locations corresponding to burn/unburn group
Forest_burned_sif_map<-state_plot+
  geom_sf(data = states, fill = "white")+
  geom_sf(data = Mboundary, color = NA, fill = "grey", alpha=1)+
  geom_sf(data = Bboundary, color = "black", fill = NA, alpha=0.8)+
  geom_point(data=df_b_loc_final, aes(x=lon, y=lat, colour=group,))+
  coord_sf(xlim = c(rlon[1], rlon[2]-6), ylim = c(rlat[1]+5, rlat[2]), expand = FALSE)+
  scale_color_manual(values = c("burned" = "black", "unburned" = "lightgreen"))+
  theme(panel.background = element_rect(fill = "slategray3",colour = "slategray3"))+mytheme

ggplotly(Forest_burned_sif_map)

################### Visualize the SIF time series for burned/unburned locations  ###################
#set the date of the fire
fireday<-paste("2018-07-27")
#set theme for plotting 
mytheme=theme_light()+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

sif_time_plot<-ggplot(data=df_m2, aes(x=date, y=sif, group=lonlat, color = group, label=lon, label2=lat), alpha=0.5)+
  geom_line()+geom_vline(xintercept = as.Date(fireday), linetype = "dashed")+labs(title="", y="SIF", x="Time", color="Group")+scale_color_manual(labels = c("Burned", "Unburned"), values=c('black','lightgreen'))+
  ylim(-1,1.5)+mytheme+
  scale_x_date(breaks = "4 months")
  

sif_time_plot

ggplotly(sif_time_plot)

#combined plots 
sif_interactive<-subplot(Forest_burned_sif_map, sif_time_plot, nrows = 1, margin = 0.05, heights=0.7, widths=c(0.5,0.5))


############# Visualize Mean SIF over time between the two groups #################
df_mean_sif_time<-df_m2 %>% group_by(group, date) %>% summarise(sif=mean(sif, na.rm=TRUE))
df_mean_sif_time<-data.frame(mean_sif_time)

mean_sif_time_plot<-ggplot(data=df_mean_sif_time, aes(x=date, y=sif, color = group), alpha=1)+
  geom_line()+labs(title="Mean SIF", y="SIF", x="Time", color="Group")+scale_color_manual(labels = c("Burned", "Unburned"), values=c('black','green'))+
  ylim(-1,1.5)+mytheme+scale_x_date(breaks = "4 months")+
  geom_vline(xintercept = as.Date(fireday), linetype = "dashed")


#### save subsetted and labeled data for statistical analysis ####
savedat<-df_m2[,-c(11)] #remove forest column since we already subsetted just the forest region
#fwrite(savedat, file="../data/mendocino_sif.csv", sep=",", row.names=FALSE, col.names = TRUE)


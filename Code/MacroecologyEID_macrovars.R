# ------------------------------------------------------------------------------
# Advancing the use of macroecological approaches for the prediction of zoonotic disease risk'
# Corresponding authors: moreno.dimarco@uniroma1.it; andrea.tonelli@uniroma1.it

library(terra); library(sf); library(raster); library(dplyr); library(rnaturalearth); library(funspace); library(ggplot2)

setwd("")

#Load multipolygon shapefile of Brazil's municipalities
mun_shp_path <- "BR_Municipios_2021.shp" # The shapefile of Brazilian municipalities (2021 edition) is publicly available from IBGE at the following link: https://www.ibge.gov.br/en/geosciences/territorial-organization/territorial-meshes/18890-municipal-mesh.html?edicao=33161&t=downloads
mun_shp <- read_sf(mun_shp_path) 
mun_shp$ID <- 1:nrow(mun_shp)
mun_shp
crs(mun_shp)

#Take Brazil extent 
brazil <- ne_countries(country = "Brazil")
ext(brazil)
crs(brazil)

# Paths to AOH data for selected mammal orders
# The Area of Habitat (AOH) maps used here come from: Lumbierres et al. (2022) – https://doi.org/10.1038/s41597-022-01838-w
# Download and extract the species-level AOH rasters for the relevant mammalian orders into the following subfolders under: 

pp <- c(
  "Mammals_carnivora",
  "Mammals_primates",
  "Mammals_cingulata",
  "Mammals_didelphimorphia",
  "Mammals_pilosa",
  "Mammals_cetartiodactyla",
  "Mammals_chiroptera",
  "Mammals_lagomorpha",
  "Mammals_perissodactyla"
)


# ------------------------------------------------------------------------------
# Host richness
# ------------------------------------------------------------------------------

for (i in 1:length(pp)){
  n=0
  lr <- list()
  #analyzed order name
  analyzed_order <- pp[i]
  analyzed_order <- strsplit(analyzed_order, split="Mammals_")
  analyzed_order <- analyzed_order[[1]][2]
  print(paste("analyzing:",analyzed_order))
  print(paste("num order",i,"over",length(pp)))
  
  #initialize output raster for Host Richness
  order_sum<- rast()
  terra::values(order_sum) <- NA
  terra::ext(order_sum) <- ext(brazil)
  terra::res(order_sum) <- c(0.009920722,0.009920722)
  
  #list of species 
  filelist_analyzed_order <- list.files(path = pp[i], pattern="\\.tif",
                                        recursive = T, full.names = TRUE)
  
  
  for (j in 1:length(filelist_analyzed_order)){
    
    print(j)
    suppressWarnings(rm(int_rshp))
    
    #load raster
    r <- rast(filelist_analyzed_order[[j]]) 
    
    #aggregate of fact 10: from 100 m to 1 km resolution
    r_resample <- terra::aggregate(r, fact=10, fun=max, na.rm = TRUE)
    crs(r_resample) <- crs(brazil)
    
    #intersection
    int_rshp <- terra::intersect(ext(r_resample), ext(brazil))
    
    #If there are NA values only, skip 
    if (all(is.na(terra::values(r_resample))) == TRUE){
      next
    }
    #If the extents do not overlap, skip
    if (is.null(int_rshp)==T) {
      next
    }
    
    #crop and aggregate rasters
    r_crop <- terra::crop(r_resample, ext(brazil))
    r_crop <- terra::resample(r_crop, order_sum, "near")
    stack <- c(order_sum, r_crop)
    order_sum <- sum(stack, na.rm = TRUE)
    
    #fix NaN and Infinite values
    order_sum[values(order_sum)<1] <- 0
    order_sum[!is.finite(values(order_sum))] <- 0
    
    lr <- c(lr,r_crop)
    n=n+1
    
  }
  
  #export Order Richness
  order_sum <- crop(order_sum, ext(brazil), extend=F)
  rast_out_name <- paste0(getwd(), "/Mammals/", sep="")
  rast_out_name <- paste0(rast_out_name,analyzed_order,sep="")
  rast_out_name <- paste(rast_out_name,".tif",sep="")
  terra::writeRaster(order_sum, rast_out_name, overwrite=T)
  print(n)
  print(paste("finish", analyzed_order))
  
}


# list and load all rasters
hostlist<- list.files(path = paste(getwd(), "/Mammals",sep=""), 
                      pattern="\\.tif",
                      recursive = F, full.names = TRUE)

# initialize df for host richness
df_HR <- data.frame(mun_shp$CD_MUN)
df_HR$ID <- mun_shp$ID
head(df_HR)
colnames(df_HR)[1] <- "CD_MUN"


for (i in 1:length(hostlist)){
  r <- rast(hostlist[i])
  crs(r) <- crs(mun_shp)
  
  #extract values
  r_in_mun <- terra::extract(r, mun_shp, method='simple', df=T)
  
  #averaging by ID mun 
  r_in_mun_averaged <- r_in_mun %>% group_by(ID) %>% 
    summarise(across(everything(),  ~mean(.x, na.rm=TRUE)))
  
  df_HR <- left_join(df_HR, r_in_mun_averaged, by="ID")
  a <- strsplit((hostlist[i]), split="/")
  a <- strsplit(a[[1]][7], split=".tif")
  colnames(df_HR)[i+2] <- a[[1]][1]
}


#host richness dataframe
df_HR$Host_rich <- apply(df_HR[,-c(1,2)], 1, FUN=sum, na.rm=T)
df_HR$Host_rich <- round(df_HR$Host_rich,2)



# ------------------------------------------------------------------------------
# Functional richness
# ------------------------------------------------------------------------------

#create output dataframe for FR
df_species <- st_drop_geometry(mun_shp[, c("CD_MUN", "ID")])


for (i in 1:length(pp)){
  n=0
  lr <- list()
  # analyzed order name
  analyzed_order <- pp[i]
  analyzed_order <- strsplit(analyzed_order, split="Mammals_")
  analyzed_order <- analyzed_order[[1]][2]
  print(paste("analyzing:",analyzed_order))
  print(paste("num order",i,"over",length(pp)))
  
  # list of species 
  filelist_analyzed_order <- list.files(path = pp[i], pattern="\\.tif",
                                        recursive = T, full.names = TRUE)
  
  
  for (j in 1:length(filelist_analyzed_order)){
    
    print(j)
    suppressWarnings(rm(int_rshp))
    
    # load raster
    r <- rast(filelist_analyzed_order[[j]])
    
    # aggregate of fact 10: from 100 m to 1 km resolution
    r_resample <- terra::aggregate(r, fact=10, fun=max, na.rm = T)
    
    # intersection
    int_rshp <- terra::intersect(ext(r_resample), ext(brazil))
    
    # if there are NA values only, skip 
    if (all(is.na(terra::values(r_resample))) == TRUE){
      next
    }
    # if the extents do not overlap, skip
    if (is.null(int_rshp)==T) {
      next
    }
    
    
    # extract values within municipalities
    r_in_mun <- terra::extract(r_resample, mun_shp, method='simple', df=TRUE)#, ID=T)
    
    # get whether species is present in each municipality
    r_in_mun_presence <- r_in_mun %>%    distinct() %>%    na.omit()
    
    # retrieve MUN_ID
    r_in_mun_presence_munid <- merge(mun_shp[c(1, 6)] %>% st_drop_geometry(), r_in_mun_presence, by = "ID") %>% 
      mutate()
    
    # combine to the final df
    df_species <- left_join(df_species, r_in_mun_presence_munid %>% dplyr::select(-ID), by = "CD_MUN")
    n=n+1
  }
  
  
}


### FR COMPUTATION ###

#Fix species name
names(df_species) <- gsub("ex_", "", names(df_species))
names(df_species) <- gsub("_", " ", names(df_species))
names(df_species)
traits_all <- read.csv("trait_data.csv") #from COMBINE database (Soria et al., 2021 -  https://doi.org/10.1002/ecy.3344)

#select the desired traits
traits_toconsider <- c("iucn2020_binomial",
                       "adult_mass_g",
                       "adult_body_length_mm",
                       "max_longevity_d",
                       "age_first_reproduction_d",
                       "gestation_length_d",
                       "litter_size_n",
                       "litters_per_year_n",
                       "interbirth_interval_d",
                       "weaning_age_d",
                       "generation_length_d",
                       "dispersal_km",
                       "det_diet_breadth_n",
                       "trophic_level",
                       "activity_cycle",
                       "habitat_breadth_n")

#subsetting for desired traits
traits <- traits_all[,colnames(traits_all) %in% traits_toconsider]


# subsetting for Brazilian host species 
# this bit is to remove species with NAs in all municipalities (i.e., not present in Brazil) that slipped through the extent-match filter: colSums(is.na(df_species)) < nrow(df_species)])
traits <- traits[traits$iucn2020_binomial %in% colnames(df_species[, colSums(is.na(df_species)) < nrow(df_species)]),]


# checking for duplicates
length(unique(traits$iucn2020_binomial))
# take unique value based on species name
traits <- traits %>%
  group_by(iucn2020_binomial) %>%
  summarise_if(is.numeric,mean, na.rm = TRUE)#, order=(order))
head(traits)
df_traitsBrazil <- traits

# scaling
df_traitsBrazil_scaled <- scale(df_traitsBrazil[,-1])

# imputation of missing traits
imputed_traits <- funspace::impute(df_traitsBrazil_scaled, messages = F)
imputed_traits <- as.data.frame(imputed_traits$imputed)
rownames(imputed_traits) <- df_traitsBrazil$iucn2020_binomial
head(imputed_traits)
dim(imputed_traits)

# PCA of Brazil
BRAZIL_PCA <- princomp(imputed_traits)
summary(BRAZIL_PCA)


# explained variance
ExpVar <- BRAZIL_PCA$sdev^2/sum(BRAZIL_PCA$sdev^2)*100

sum(ExpVar[1:2])
sum(ExpVar[1:3])
sum(ExpVar[1:4]) ##need to consider 4 PC

# create output dataframe
df_FR <- df_species[,1:2]
df_FR$FR12 <- NA
df_FR$FR13 <- NA
df_FR$FR23 <- NA
df_FR$FR14 <- NA
df_FR$FR24 <- NA
df_FR$FR34 <- NA


# MAIN CYCLE (loop over municipalities)
for (i in 1:nrow(df_FR)) { 
  #subset for species that are actually present
  c_df <- as.data.frame(df_species[i,])
  c_df <- c_df[,!is.na(c_df)]
  
  list.of.species.present <- colnames(c_df)[-c(1,2)]
  imputed_traits$inMunicipality <- ifelse((rownames(imputed_traits) %in%  list.of.species.present),
                                          "Yes", "No")
  #FD using PC1 and PC2
  FR_municipality12 <- funspace(BRAZIL_PCA, PCs=c(1,2), group.vec = imputed_traits$inMunicipality)
  df_FR$FR12[i] <- round(FR_municipality12$FD$groups$FRich[2],2)
  
  #FD using PC1 and PC3
  FR_municipality13 <- funspace(BRAZIL_PCA, PCs=c(1,3), group.vec = imputed_traits$inMunicipality)
  df_FR$FR13[i] <- round(FR_municipality13$FD$groups$FRich[2],2)
  
  #FD using PC2 and PC3
  FR_municipality23 <- funspace(BRAZIL_PCA, PCs=c(2,3), group.vec = imputed_traits$inMunicipality)
  df_FR$FR23[i] <-  round(FR_municipality23$FD$groups$FRich[2],2)
  
  #FD using PC1 and PC4
  FR_municipality14 <- funspace(BRAZIL_PCA, PCs=c(1,4), group.vec = imputed_traits$inMunicipality)
  df_FR$FR14[i] <-  round(FR_municipality14$FD$groups$FRich[2],2)
  
  #FD using PC2 and PC4
  FR_municipality24 <- funspace(BRAZIL_PCA, PCs=c(2,4), group.vec = imputed_traits$inMunicipality)
  df_FR$FR24[i] <-  round(FR_municipality24$FD$groups$FRich[2],2)
  
  #FD using PC3 and PC4
  FR_municipality34 <- funspace(BRAZIL_PCA, PCs=c(3,4), group.vec = imputed_traits$inMunicipality)
  df_FR$FR34[i] <-  round(FR_municipality34$FD$groups$FRich[2],2)
  
  print(i)
}


#sum functional richness within each municipality

df_FR <- df_FR %>%
  rowwise() %>%
  mutate(Fun_rich = sum(c_across(starts_with("FR")), na.rm = TRUE)) %>%
  ungroup()

df_FR <- df_FR[c("CD_MUN", "ID", "Fun_rich")]





df = readRDS("~/Desktop/nSPC.site.2008_2023.rds")

df <- df %>%
  mutate(YEAR = coalesce(as.integer(ANALYSIS_YEAR), OBS_YEAR))

df = df %>% 
  filter(REGION == "MHI") %>% 
  select(REGION, ISLAND, YEAR, LONGITUDE, LATITUDE, MEAN_DEPTH, SPECIES, abund.site)

df_all = NULL

for (s in 1:length(unique(df$SPECIES))) {
  
  df_i = df
  
  df_i$response = ifelse(df_i$SPECIES == unique(df$SPECIES)[s], df_i$abund.site, 0)
  
  df_i = df_i %>% 
    group_by(ISLAND, YEAR, LONGITUDE, LATITUDE, MEAN_DEPTH) %>% 
    summarise(response = sum(response, na.rm = T))
  
  colnames(df_i) = c("island", "year", "lon", "lat", "depth", unique(df$SPECIES)[s])
  
  
  if (s == 1) {
    
    df_all = df_i
    
  }else{
    
    df_all = left_join(df_all, df_i)
    
  }
  
}

df = df_all

plot(df$lon, df$lat)

save(df, file = "~/al_ensemble/ncrmp_all_sp_clean.Rdata")


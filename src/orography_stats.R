head(topologyInfo)
statsTop<- bind_rows(read.csv("../../Structural change/dmel_10kb/wd_mountains.bed",
                              sep = "\t", header = FALSE) %>% 
                              mutate(species = "Dmel", resolution = "10kb"),
                     read.csv("../../Structural change/dmel_5kb/wd_mountains.bed",
                              sep = "\t", header = FALSE)%>% 
                              mutate(species = "Dmel", resolution = "5kb"),
                     read.csv("../../Structural change/dsim_10kb/wd_mountains.bed",
                              sep = "\t", header = FALSE) %>% 
                              mutate(species = "Dsim", resolution = "10kb"),
                     read.csv("../../Structural change/dsim_5kb/wd_mountains.bed",
                              sep = "\t", header = FALSE)%>% 
                              mutate(species = "Dsim", resolution = "5kb"),
                     read.csv("../../Structural change/dmel_10kb/wd_valleys.bed",
                              sep = "\t", header = FALSE) %>% 
                       mutate(species = "Dmel", resolution = "10kb"),
                     read.csv("../../Structural change/dmel_5kb/wd_valleys.bed",
                              sep = "\t", header = FALSE)%>% 
                       mutate(species = "Dmel", resolution = "5kb"),
                     read.csv("../../Structural change/dsim_10kb/wd_valleys.bed",
                              sep = "\t", header = FALSE) %>% 
                       mutate(species = "Dsim", resolution = "10kb"),
                     read.csv("../../Structural change/dsim_5kb/wd_valleys.bed",
                              sep = "\t", header = FALSE)%>% 
                       mutate(species = "Dsim", resolution = "5kb")) %>%
                    select(V1,V2,V3,V4,species,resolution) %>%
                    rename(chr = V1, start = V2, end = V3, region_type = V4) %>%
                    mutate(size = end - start)

totalSize <- statsTop %>% group_by(species, resolution, region_type) %>% summarise(sum(size)/1e6)
valleySizes <- statsTop %>% 
                filter(region_type == "Valley") %>%
                group_by(species, resolution) %>%
                summarise(quantile(size, 0.50)/1e3)



library(plotly)


ggplot(statsTop %>% filter(region_type == "Valley"),
                aes(y = size/1e3, fill = resolution)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cadetblue4","darkorchid4"))+
  theme_pubclean() + 
  labs(y = "size (kb)") +
  facet_wrap(~species) +
  theme(legend.position = "bottom")

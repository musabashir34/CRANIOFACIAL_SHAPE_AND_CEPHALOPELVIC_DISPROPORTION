# error_bubble_plot.R

# 1. Load the per-landmark error metrics from your directory CSV file
ant_error_metrics <- read.csv("anterior_per_landmark_error_geomorph.csv")

# 2. Extract spatial 2D baseline tracking anchors from your adjusted control mean configuration
# This provides the physical X and Y coordinate planes matching the face profile layout
ant_spatial_coordinates <- as.data.frame(ant_C_adj)
colnames(ant_spatial_coordinates) <- c("X_Coord", "Y_Coord")
ant_spatial_coordinates$landmark <- 1:nrow(ant_spatial_coordinates)

# 3. Create the definitive localized anterior anatomical label dictionary map vector
# Matches your 15 confirmed labels precisely row-for-row
ant_landmark_anatomical_labels <- c(
  "1: Glabella", "2: Pronasale", "3: Labiale sup.", "4: Labiale inf.", 
  "5: Gnathion", "6: L.exocanthion", "7: R.exocanthion", "8: L.endocanthion", 
  "9: R.endocanthion", "10: L.alare", "11: R.alare", "12: L.cheilion", 
  "13: R.cheilion", "14: L.zygion", "15: R.zygion"
)

# 4. Bind the tables together into a single unified plotting dataframe
plot_ant_error_df <- ant_error_metrics %>%
  left_join(ant_spatial_coordinates, by = "landmark") %>%
  mutate(Anatomical_Label = ant_landmark_anatomical_labels[landmark])

# 5. Build the Publication-Grade Bubble Diagram using ggplot2
gg_ant_bubble_err <- ggplot(plot_ant_error_df, aes(x = X_Coord, y = Y_Coord)) +
  # Render the data point bubbles with custom border weights
  geom_point(aes(size = mean_sq_replicate_disp, fill = mean_sq_replicate_disp), 
             shape = 21, color = "black", alpha = 0.85) +
  
  # Scale colors dynamically from low noise (blue) to peak error (red) matching the sample image
  scale_fill_gradientn(colors = c("dodgerblue3", "purple", "hotpink3", "tomato1")) +
  
  # Adjust bubble diameters to preserve layout spacing across your data ranges
  scale_size_continuous(range = c(4, 14), name = "Mean sq.\nreplicate\ndisplacement") +
  
  # Inject text flags using ggrepel to prevent overlapping labels near the center face clusters
  geom_text_repel(aes(label = Anatomical_Label), 
                  size = 3.2, fontface = "bold", box.padding = 0.4, 
                  point.padding = 0.3, segment.color = "grey60") +
  
  labs(
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    panel.grid.major = element_line(color = "grey93"),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(), # Suppresses raw coordinate metrics on margins
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(fill = "none") # Unifies size and color channels into one legend layout



# 6. Save the high-resolution vector layout file 
if(!dir.exists("figures_publication")) dir.create("figures_publication")
ggsave("figures_publication/Figure_04_Anterior_Per_Landmark_Error_Map.png", 
       plot = gg_ant_bubble_err, width = 10, height = 8, dpi = 300)

#=================LATERAL ERROR PLOT====================================

# 1. Load the per-landmark error metrics from your directory CSV file
lat_error_metrics <- read.csv("lateral_per_landmark_error_geomorph.csv")

# 2. Extract spatial 2D baseline tracking anchors from your adjusted control mean configuration
# This provides the physical X and Y coordinate planes matching the skull profile layout
lat_spatial_coordinates <- as.data.frame(lat_C_adj)
colnames(lat_spatial_coordinates) <- c("X_Coord", "Y_Coord")
lat_spatial_coordinates$landmark <- 1:nrow(lat_spatial_coordinates)

# 3. Create the definitive localized anatomical label dictionary map vector
landmark_anatomical_labels <- c(
  "1: Glabella", "2: Nasion", "3: Pronasale", "4: Subnasale", 
  "5: Labiale sup.", "6: Labiale inf.", "7: Pogonion", "8: Gnathion", 
  "9: Tragus", "10: Sup. ear attach.", "11: Opisthocranion", 
  "12: Lambda", "13: Vertex", "14: Bregma"
)

# 4. Bind the tables together into a single unified plotting dataframe
plot_error_df <- lat_error_metrics %>%
  left_join(lat_spatial_coordinates, by = "landmark") %>%
  mutate(Anatomical_Label = landmark_anatomical_labels[landmark])

# 5. Build the Publication-Grade Bubble Diagram using ggplot2
gg_bubble_err <- ggplot(plot_error_df, aes(x = X_Coord, y = Y_Coord)) +
  # Gridlines are generated via panel grids; background fields are kept clean
  geom_point(aes(size = mean_sq_replicate_disp, fill = mean_sq_replicate_disp), 
             shape = 21, color = "black", alpha = 0.85) +
  
  # Scale colors dynamically from low noise (blue/purple) to peak error (red/orange)
  scale_fill_gradientn(colors = c("dodgerblue3", "purple", "hotpink3", "tomato1")) +
  
  # Adjust bubble diameters to prevent overlay clutter while keeping sizes distinguishable
  scale_size_continuous(range = c(4, 14), name = "Mean sq.\nreplicate\ndisplacement") +
  
  # Inject text flags using ggrepel to handle close clusters like 7: Pogonion and 8: Gnathion safely
  geom_text_repel(aes(label = Anatomical_Label), 
                  size = 3.2, fontface = "bold", box.padding = 0.4, 
                  point.padding = 0.3, segment.color = "grey60") +
  
  labs(
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    panel.grid.major = element_line(color = "grey93"),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(), # Cleans up raw spatial pixel unit tags from the margins
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(fill = "none") # Combines size and color indicators into a single clean legend array

# Display plot in your knitted report


# 6. Save the high-resolution vector layout file 
if(!dir.exists("figures_publication")) dir.create("figures_publication")
ggsave("figures_publication/Figure_03_Lateral_Per_Landmark_Error_Map.png", 
       plot = gg_bubble_err, width = 10, height = 8, dpi = 300)

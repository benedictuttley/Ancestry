# Load Dependencies
library(tidyr)
packages <- c("tidyverse",  "RColorBrewer")
install.packages(packages)
lapply(packages, library, character.only = TRUE)

# Core file paths
eigenvalues_fp <- "../results/combined.eigenval"
eigenvectors_fp <- "../results/combined.eigenvec"
population_fp <- paste("../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chr",
"all-mac5-v2.population", sep="")
super_population_fp <- paste("../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phas",
"e3-chrall-mac5-v2.super-population", sep = "")

# Produce Scree plot
eigenvalues <- read.table(eigenvalues_fp, 
                          header = F, 
                          col.names = c("Eigenvalue")) %>% 
  as_tibble()

  scree_plot <- ggplot(data=eigenvalues, 
                       aes(x = seq(1, nrow(eigenvalues), 1), y = Eigenvalue)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(1, nrow(eigenvalues), 1)) +
    labs(x = "Principal Component", y = "Eigenvalue") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  ggsave(file = "scree.png", plot = scree_plot, width=12, height=8)
  
  
########## MODULE 5 START ##########

# STEP 1: Calculate the mean and SD for each PC in individuals tagged as GBR
eigenvectors <- read.table(eigenvectors_fp, header = T) %>% 
  as_tibble()

population_codes <- read.table(population_fp, header = F, 
                               col.names = c("FID", "IID", "POP")) %>% 
  as_tibble()

super_population_codes <- read.table(super_population_fp, header = F, 
                                     col.names = c("FID", "IID", "SUPPOP")) %>% 
  as_tibble()


GBR_tagged_individuals <- inner_join(eigenvectors,
                                     population_codes,
                                     by="FID") %>% 
  filter(POP=="GBR")

# Generate inclusion criteria for a principal component (mean +/- 4SD)
CalculateBounds <- function(pc_column) {
  mean <- mean(pc_column)
  sd <- sd(pc_column)
  lower_bound <- mean - (4*sd)
  upper_bound <- mean + (4*sd)
  tibble(UpperBound=upper_bound, LowerBound=lower_bound)
}

pc_subtable <- GBR_tagged_individuals %>% 
  select(PC1,PC2,PC3)

pc_bounds <- pc_subtable %>% 
  purrr::map_df(CalculateBounds) %>% 
  add_column(PC=colnames(pc_subtable), .before = "UpperBound")

# Exclude any test individual who's PC's are outside the GBR-like limits
test_datset_eigenvectors <- eigenvectors %>% 
  filter(IID==1)
for (pc in colnames(pc_subtable)){
  target <- filter(pc_bounds, PC==pc)
  test_datset_eigenvectors <- filter(test_datset_eigenvectors, 
                                     !!as.symbol(pc) <= target$UpperBound,
                                     !!as.symbol(pc) >= target$LowerBound)
  }

# STEP 2: Create the .keep file
# NOTE: You may need to manually change permissions on file to view it
# TODO: Must be tab-delimited
test_datset_eigenvectors %>% 
  select(FID, IID) %>% 
  write_delim("GBR-like.keep")

sprintf("%d individuals with ancestry similar to the GBR reference genomes", 
        nrow(test_datset_eigenvectors))

merged <- merge(eigenvectors, population_codes, by=c("FID", "IID"), all = T) %>% 
  merge(super_population_codes, by=c("FID", "IID"), all = T) %>%
  replace(., is.na(.), "test") %>% 
  arrange(SUPPOP, POP)


# STEP 3: Plot the PCs
core_colors <- c("Oranges", "Blues", "Greens", "Purples", "Reds")
color_vector <- vector()

# All ancestry
all_ancestry <- merged %>% filter(POP != "test")
# All test
all_test <- merged %>% filter(POP == "test")
# GBR-like test
GBR_like_test <- test_datset_eigenvectors
nrow(GBR_like_test)
all_sup_pops <- all_ancestry %>% 
  select(SUPPOP) %>% 
  distinct() %>% 
  as_vector()

all_pops <- all_ancestry %>% 
  select(POP) %>% 
  distinct() %>% 
  as_vector()

for(sup_pop in all_sup_pops){
  
  num_pops <- all_ancestry %>% 
    filter(SUPPOP==sup_pop) %>% 
    select(POP) %>% 
    distinct() %>% 
    count() %>%
    as.integer()

  color_vector <- c(color_vector, 
                    rev(brewer.pal(n=9, name = core_colors[1]))[1:num_pops])
  core_colors <- core_colors[-1]
}

groups_to_plot <- c(all_pops, "test", "GBR-like")
colors_vector <- c(color_vector, "#FFFFFF", "#FFFF00")
shapes_vector <- c(rep(c(21), length(all_pops)), 22, 22)

GBR_like_test$POP <- rep(c("GBR-like"), times=nrow(GBR_like_test))

GenerateAncestryPlot <- function(x, y){
 plt <- ggplot()+ 
    geom_point(data=all_ancestry,
               color="black",
               size=3,
               aes_string(x=x, y=y, fill="POP", shape="POP"))+
   
    geom_point(data=all_test,
               color="black",
               size=3,
               aes_string(x=x, y=y, fill="POP", shape="POP"))+

   geom_point(data=GBR_like_test,
              color="black",
              size=3,
              aes_string(x=x, y=y, fill="POP", shape="POP"))+

   scale_fill_manual(values=colors_vector, breaks=groups_to_plot) +
   scale_shape_manual(values=shapes_vector, breaks=groups_to_plot) +
   
   theme(
     plot.background = element_rect(fill = "white"),
     panel.background = element_rect(fill = "white"),
     legend.title = element_blank())
 return(plt)
 }

# Generate plots
pc1_vs_pc2 <- GenerateAncestryPlot("PC1", "PC2")
pc2_vs_pc3 <- GenerateAncestryPlot("PC2", "PC3")
pc1_vs_pc3 <- GenerateAncestryPlot("PC1", "PC3")

# Save plots
ggsave(file = "pc1_vs_pc2.png", plot = pc1_vs_pc2, width=12, height=8)
ggsave(file = "pc2_vs_pc3.png", plot = pc2_vs_pc3, width=12, height=8)
ggsave(file = "pc1_vs_pc3.png", plot = pc1_vs_pc3, width=12, height=8)
########## MODULE 5 END ##########

########## MODULE 6 START ##########

test_females <- read.table("../results/met583-test-female.keep",
                           header = TRUE) %>%
  as_tibble()

GBR_like_females <- merge(test_females, GBR_like_test,
                          by=c("FID", "IID"), 
                          all = FALSE)

sprintf("%d females with ancestry similar to the GBR reference genomes", 
        nrow(GBR_like_females))

GBR_like_females %>% 
  select(FID, IID) %>% 
  write_delim("GBR-like-females.keep", delim = "\t")

########## MODULE 6 END ##########



########## ADDITIONAL ANALYSIS FOR REPORT START ##########

# EUR ONLY PLOT
eur_ancestry <- all_ancestry %>% filter(SUPPOP=="EUR")
pops <- eur_ancestry %>% select(POP) %>% distinct() %>% as_vector()
num_pops <- length(pops)
color_vector <- rev(brewer.pal(n=8, name = "Dark2"))[1:num_pops]
shapes_vector <- rep(c(21), num_pops)
GenerateEurOnlyAncestryPlot <- function(x, y){
  plt <- ggplot()+ 
    geom_point(data=eur_ancestry,
               color="black",
               size=3,
               aes_string(x=x, y=y, fill="POP", shape="POP"))+
    scale_fill_manual(values=color_vector, breaks = pops) +
    scale_shape_manual(values= shapes_vector) + 
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.title = element_blank())
  return(plt)
  }


pc1_vs_pc2_eur_only <- GenerateEurOnlyAncestryPlot("PC2", "PC1")
pc2_vs_pc3_eur_only <- GenerateEurOnlyAncestryPlot("PC2", "PC3")
pc1_vs_pc3_eur_only <- GenerateEurOnlyAncestryPlot("PC1", "PC3")

# Save plots
ggsave(file = "pc1_vs_pc2_eur_only.png", plot = pc1_vs_pc2_eur_only, width=12,
       height=8)
ggsave(file = "pc2_vs_pc3_eur_only.png", plot = pc2_vs_pc3_eur_only, width=12,
       height=8)
ggsave(file = "pc1_vs_pc3_eur_only.png", plot = pc1_vs_pc3_eur_only, width=12,
       height=8)

super_population_codes %>%group_by(SUPPOP) %>% summarize(count=n())
########## ADDITIONAL ANALYSIS FOR REPORT END ##########

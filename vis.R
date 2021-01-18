# Load Dependencies
library(tidyr)
# TODO: We have redundancy here
# TODO: Set up GitHub
packages <- c("tidyverse","dplyr", "ggplot2", "purrr", "RColorBrewer")
install.packages(packages)
lapply(packages, library, character.only = TRUE)
getwd()
v <- read.table("as.raw", header = TRUE,)
prcomp(na.omit(v[7:10]))
sum(diag(cov(na.omit(v[7:10]))))
eigen(cov(na.omit(v[7:10])))$values
# Core file paths
eigenvalues_fp <- "../results/combined_alt.eigenval"
eigenvectors_fp <- "../results/combined_alt.eigenvec"
population_fp <- paste("../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chr",
"all-mac5-v2.population", sep="")
super_population_fp <- paste("../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phas",
"e3-chrall-mac5-v2.super-population", sep = "")

# <<< MODULE 4 >>>, Step 2, Produce Scree plot
eigenvalues <- read.table(eigenvalues_fp, 
                          header = F, 
                          col.names = c("Eigenvalue")) %>% 
  as_tibble()

print(eigenvalues)
print(nrow(eigenvalues))
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
  
  scree_plot
  ggsave(file = "scree.png", plot = scree_plot, width=12, height=8)



# <<< MODULE 5 >>>, Step 1, PC mean and SD

# Calculate the mean and SD for each PC in individuals tagged as GBR

eigenvectors <- read.table(eigenvectors_fp, header = T) %>% 
  as_tibble()

population_codes <- read.table(population_fp, header = F, 
                               col.names = c("FID", "IID", "POP")) %>% 
  as_tibble()

super_population_codes <- read.table(super_population_fp, header = F, 
                                     col.names = c("FID", "IID", "SUPPOP")) %>% 
  as_tibble()



# 91 GBR tagged individuals identified
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
nrow(test_datset_eigenvectors)
# <<< MODULE 5 >>>, Step 2, Create the .keep file
# NOTE: You may need to manually change permissions on file to view it
# TODO: Must be tab-delimited
test_datset_eigenvectors %>% 
  select(FID, IID) %>% 
  write_delim("GBR-like.keep")

merged <- merge(eigenvectors, population_codes, by=c("FID", "IID"), all = T) %>% 
  merge(super_population_codes, by=c("FID", "IID"), all = T) %>%
  replace(., is.na(.), "test") %>% 
  arrange(SUPPOP, POP)


# <<< MODULE 5 >>>, Step 3, Plot the PCs
core_colors <- c("Oranges", "Blues", "Greens", "Purples", "Reds")
color_vector <- vector()

# TODO: Unify these better
# All ancestry
all_ancestry <- merged %>% filter(POP != "test")
# All test
all_test <- merged %>% filter(POP == "test")
# GBR-like test
GBR_like_test <- test_datset_eigenvectors

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
pc1_vs_pc2 <- GenerateAncestryPlot("PC1", "PC2") # PC1 vs PC2
pc2_vs_pc3 <- GenerateAncestryPlot("PC2", "PC3") # PC2 vs PC3
pc1_vs_pc3 <- GenerateAncestryPlot("PC1", "PC3") # PC1 vs PC3

# Save plots
ggsave(file = "pc1_vs_pc2.png", plot = pc1_vs_pc2, width=12, height=8)
ggsave(file = "pc2_vs_pc3.png", plot = pc2_vs_pc3, width=12, height=8)
ggsave(file = "pc1_vs_pc3.png", plot = pc1_vs_pc3, width=12, height=8)

# Current: 190 lines

# TODO
# <<<< MODULE 6 >>>
# Keep intersect of GBR-like with identified as female

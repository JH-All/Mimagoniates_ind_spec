# Packages --------------
library(readxl)
library(tidyverse)
library(vegan)
library(RInSp)
library(ggthemes)
library(igraph)
library(pairwiseAdonis)
library(R2jags)
library(loo)
library(ggthemes)
library(see)
library(betareg)
library(car)
library(adespatial)
library(purrr)
library(FSA)
library(pairwiseAdonis)
library(RVAideMemoire)

# Getting data ready ------
data = read_excel("data.xlsx")
str(data)
data$S = specnumber(data[,13:34])

data_filtered = 
  data %>% 
  filter(S > 0 )

data_filtered %>%
  group_by(Env) %>%
  summarise(n = n())

data_filtered$Env <- factor(data_filtered$Env, levels = c("Stream", "Pool", "Ditch"))

data_filtered <- data_filtered %>%
  dplyr::filter(Env == "Stream") %>%
  bind_rows(data %>% dplyr::filter(Env == "Pool")) %>%
  bind_rows(data %>% dplyr::filter(Env == "Ditch"))

data_filtered = 
  data_filtered %>% 
  filter(S > 0 )

data_filtered %>%
  group_by(Env) %>%
  summarise(
    min_SL = min(SL, na.rm = TRUE),
    max_SL = max(SL, na.rm = TRUE),
    mean_SL = mean(SL, na.rm = TRUE),
    sd_SL = sd(SL, na.rm = TRUE),
    n = n()
  )

stream = data_filtered %>% 
  filter(Env == "Stream")

stream_diet = stream[,13:34]

pool = data_filtered %>% 
  filter(Env == "Pool")

pool_diet = pool[,13:34]

ditch  = data_filtered %>% 
  filter(Env == "Ditch")

ditch_diet = ditch[,13:34]

stream_InSp = import.RInSp( stream_diet,  col.header = TRUE ,  row.names = 0 ,  
                           dec= ".") 
pool_InSp = import.RInSp( pool_diet ,  col.header = TRUE ,  row.names = 0 ,  
                         dec= ".") 
ditch_InSp = import.RInSp(ditch_diet ,  col.header = TRUE ,  row.names = 0 ,  
                         dec= ".") 

# Frequency of occurrence ---------------------
data_filtered$Env <- factor(data_filtered$Env)
diet_df <- cbind(Env = data_filtered$Env, diet_total)
diet_long <- diet_df %>%
  pivot_longer(-Env, names_to = "Item", values_to = "Presence")
diet_fo <- diet_long %>%
  group_by(Env, Item) %>%
  summarise(FO = 100 * mean(Presence), .groups = "drop") 

diet_fo$Item = as.factor(diet_fo$Item)
levels(diet_fo$Item)

diet_fo %>%
  filter(FO > 0) %>%   
  group_by(Env) %>%
  summarise(n_itens = n())


## Figure 2 -------------
freq = ggplot(diet_fo, aes(x = reorder(Item, -FO), y = FO, fill = Env)) +
  geom_bar(stat = "identity", position = "dodge", color = "black",
           show.legend = T) +
  labs(x = "Diet item", y = "Frequency of occurrence (%)",
       fill = NULL) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00")) +
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits = c(0,80))+
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.8, "lines")      
  )

ggsave("freq_occurrence.png", plot = freq,
       width = 10, height = 8, dpi = 300, units = "in")

# PERMANOVA & Niche breadth --------------------
diet_total = decostand(data_filtered[,13:34], method = "pa")
diss_mat <- vegdist(diet_total, method = "jaccard", binary = TRUE)
perm = adonis2(diss_mat ~ Env, data = data_filtered, permutations = 999)
perm
pairwise.perm.manova(diss_mat, data_filtered$Env, nperm = 999)

betadisp <- betadisper(diss_mat, group = data_filtered$Env)
anova(betadisp) 
permutest(betadisp, permutations = 999) 

levins_B <- function(x) {
  p <- x / sum(x)
  B <- 1 / sum(p^2)
  Bsta <- (B - 1) / (length(p) - 1)
  return(Bsta)
}

levins_values <- apply(diet_total, 1, levins_B)

df_levins <- data.frame(
  Env = data_filtered$Env,            
  Levins = levins_values  
)

summary_levins <- df_levins %>%
  group_by(Env) %>%
  summarise(
    mean_levins = mean(Levins, na.rm = TRUE),
    se_levins = sd(Levins, na.rm = TRUE) / sqrt(n())
  )

## Figure 3 -----------------
levins = ggplot(summary_levins, aes(x = Env, y = mean_levins, fill = Env)) +
  geom_col(width = 0.6, color = "black",
           show.legend = F, size = 1.5) +
  geom_errorbar(aes(ymin = mean_levins - se_levins, ymax = mean_levins + se_levins),
                width = 0.2, size = 1.1) +
  labs(x = NULL, y = "Levins' standardized niche breadth") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.11))

ggsave("Levins.png", plot = levins,
       width = 6, height = 4, dpi = 300, units = "in")

kruskal.test(Levins ~ Env, data = df_levins)
dunnTest(Levins ~ Env, data = df_levins, method = "bonferroni")


# IS ------------------------
IS_stream <- PSicalc(stream_InSp, pop.diet = "average", exclude = FALSE,
                      replicates = 1000)
# Stream: IS = 0,37, p < 0,001 

IS_pool <- PSicalc(pool_InSp, pop.diet = "average", exclude = FALSE,
                     replicates = 1000)
# Pool: IS = 0,44, p < 0,001

IS_ditch <- PSicalc(ditch_InSp, pop.diet = "average", exclude = FALSE,
                   replicates = 1000)
# Ditch: IS = 0,39, p < 0,001

env_stream <- rep("Stream", length(IS_stream$PSi))
env_pool   <- rep("Pool", length(IS_pool$PSi))
env_ditch  <- rep("Ditch", length(IS_ditch$PSi))

df_PSi <- data.frame(
  Habitat = c(env_stream, env_pool, env_ditch),
  PSi = c(IS_stream$PSi, IS_pool$PSi, IS_ditch$PSi)
)

df_PSi$Habitat <- factor(df_PSi$Habitat, 
                          levels = c("Stream", "Pool", "Ditch"))

data_filtered$PSi = df_PSi$PSi


anova_sex <- aov(PSi ~ Sex, data = data_filtered)
summary(anova_sex)

data_filtered$SL = as.numeric(data_filtered$SL)
mod_SL <- lm(PSi ~ SL, data = data_filtered)
summary(mod_SL)

anova_model <- aov(PSi ~ Habitat, data = df_PSi)
summary(anova_model)
TukeyHSD(anova_model)

## Figure 4 --------------------
hist = ggplot(df_PSi, aes(x = PSi, fill = Habitat)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7,
                 show.legend = F, size = 1.4 ) +
  facet_wrap(~ Habitat, scales = "free_y", nrow = 3) +
  labs(x = expression(PS[i]), y = "Frequency") +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00"))

ggsave("PSi_histogram_by_habitat.png", plot = hist,
       width = 8, height = 8, dpi = 300, units = "in")

# Morphometrics x PSi  ------------------
data_filtered$SL = as.numeric(data_filtered$SL)
data_filtered$Env <- fct_rev(data_filtered$Env)

data_filtered %>% 
  ggplot(aes(x = SL, y = PSi, fill = Env))+
  geom_point(shape = 21, size = 3, alpha = 0.8,
             show.legend = F)+
  geom_smooth(method = "glm", se = F,
              color = "black", linetype = "dashed",
              show.legend = F)+
  facet_wrap(~Env, nrow = 3)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00"))

data_filtered %>% 
  ggplot(aes(x = index_vental_flat, y = PSi, fill = Env))+
  geom_point(shape = 21, size = 3, alpha = 0.8,
             show.legend = F)+
  geom_smooth(method = "glm", se = F,
              color = "black", linetype = "dashed",
              show.legend = F)+
  facet_wrap(~Env, nrow = 3)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00"))

## Figure 5 ------------
rel_eye = data_filtered %>% 
  ggplot(aes(x = rel_eye_posit, y = PSi, fill = Env))+
  geom_point(shape = 21, size = 4.5, alpha = 0.8,
             show.legend = F, stroke = 1.2)+
  geom_smooth(method = "glm", se = F,
              color = "black", linetype = "dashed",
              show.legend = F, linewidth = 1.5)+
  facet_wrap(~Env, nrow = 3, scales = "free_x")+
  theme_classic(base_size = 18)+
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00"))+
  labs(x = "Relative eye position", y = expression(PS[i]))

ggsave("PSi_rel_eye_position.png", plot = rel_eye,
       width = 5, height = 8, dpi = 300, units = "in")

data_filtered$SL
data_filtered$index_vental_flat
data_filtered$rel_eye_posit

stream_data <- filter(data_filtered, Env == "Stream")

mod_stream <- betareg(PSi ~ SL + index_vental_flat + rel_eye_posit,
                      data = stream_data)
summary(mod_stream)

vif(mod_stream)

pool_data <- filter(data_filtered, Env == "Pool")

mod_pool <- betareg(PSi ~ SL + index_vental_flat + rel_eye_posit,
                    data = pool_data)
summary(mod_pool)

vif(mod_pool)

ditch_data <- filter(data_filtered, Env == "Ditch")

mod_ditch <- betareg(PSi ~ SL + index_vental_flat + rel_eye_posit,
                     data = ditch_data)
summary(mod_ditch)

vif(mod_ditch)
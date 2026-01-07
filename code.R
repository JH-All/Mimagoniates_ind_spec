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
library(R2jags)
library(loo)
library(ggdist)
library(purrr)
library(loo)
library(diptest)
library(mclust)

# Getting data ready ------
data = read_excel("data_final.xlsx")
str(data)

data <- data %>%
  mutate(across(
    .cols = -c(Env, Sex),    
    .fns  = ~ as.numeric(.x)  
  ))

data %>%
  group_by(Env) %>%
  summarise(n = n())

data$Env <- factor(data$Env, 
                            levels = c("Stream", "Pool", "Ditch"))

data %>%
  group_by(Env) %>%
  summarise(
    min_SL = min(SL, na.rm = TRUE),
    max_SL = max(SL, na.rm = TRUE),
    mean_SL = mean(SL, na.rm = TRUE),
    sd_SL = sd(SL, na.rm = TRUE),
    n = n()
  )

stream = data %>% 
  filter(Env == "Stream")

stream_diet = stream[,13:25]

pool = data %>% 
  filter(Env == "Pool")

pool_diet = pool[,13:25]

ditch  = data %>% 
  filter(Env == "Ditch")

ditch_diet = ditch[,13:25]

# Frequency of occurrence ---------------------
data$Env <- factor(data$Env)
diet_total = decostand(data[,13:25], method = "pa")
diet_df <- cbind(Env = data$Env, diet_total)
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

diet_fo$Env <- fct_rev(diet_fo$Env)

## Figure 2 -------------
fig2 = ggplot(diet_fo, aes(x = reorder(Item, -FO), y = FO, fill = Env)) +
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

ggsave("Figure_2.png", plot = fig2,
       width = 10, height = 8, dpi = 300, units = "in")

# PERMANOVA & Niche breadth --------------------
diet_total = decostand(data[,13:25], method = "pa")
diss_mat <- vegdist(diet_total, method = "jaccard", binary = TRUE)
perm = adonis2(diss_mat ~ Env, data = data, permutations = 999)
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
  Env = data$Env,            
  Levins = levins_values  
)

summary_levins <- df_levins %>%
  group_by(Env) %>%
  summarise(
    mean_levins = mean(Levins, na.rm = TRUE),
    se_levins = sd(Levins, na.rm = TRUE) / sqrt(n())
  )

summary_levins$Env <- fct_rev(summary_levins$Env)

## Figure 3 -----------------
fig3 = ggplot(summary_levins, aes(x = Env, y = mean_levins, fill = Env)) +
  geom_col(width = 0.6, color = "black",
           show.legend = F, size = 1.5) +
  geom_errorbar(aes(ymin = mean_levins - se_levins, ymax = mean_levins + se_levins),
                width = 0.2, size = 1.1) +
  labs(x = NULL, y = "Levins' standardized niche breadth") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c("Stream" = "#009E73", 
                               "Pool" = "#0072B2", 
                               "Ditch" = "#E69F00"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.17))

ggsave("Figure_3.jpg", plot = fig3,
       width = 6, height = 4, dpi = 300, units = "in")

kruskal.test(Levins ~ Env, data = df_levins)
dunnTest(Levins ~ Env, data = df_levins, method = "bonferroni")

# Bayesian PSi Models --------------------
diet_mat_perc <- as.matrix(data[, 13:25])
diet_mat_counts <- round(diet_mat_perc)
row_sums <- rowSums(diet_mat_counts)

for (i in seq_len(nrow(diet_mat_counts))) {
  diff <- 100 - row_sums[i]
  if (diff != 0) {
    j <- which.max(diet_mat_counts[i, ])
    diet_mat_counts[i, j] <- diet_mat_counts[i, j] + diff
  }
}

rowSums(diet_mat_counts)
diet_mat <- diet_mat_counts

N  <- nrow(diet_mat)
J  <- ncol(diet_mat)

data$Env <- factor(data$Env,
                            levels = c("Stream", "Pool", "Ditch"))
hab <- as.numeric(data$Env)
H <- length(levels(data$Env))  

data$Sex <- factor(data$Sex,
                            levels = c("M", "F", "I"))
sex <- as.numeric(data$Sex)
S <- length(levels(data$Sex)) 

ni <- apply(diet_mat, 1, sum)

SL <- data$SL

hSL <- hist(SL, breaks = "Sturges", plot = FALSE)
breaks_SL <- hSL$breaks

size_class <- cut(SL,
                  breaks = breaks_SL,
                  include.lowest = TRUE,
                  labels = FALSE)

size_class[is.na(size_class)] <- max(size_class, na.rm = TRUE)

K <- max(size_class)   

jags_data_base <- list(
  yi   = diet_mat,
  N    = N,
  J    = J,
  ni   = ni,
  hab  = hab,
  H    = H,
  sex  = sex,
  S    = S,
  size = size_class,
  K    = K
)

## Habitat Model ------------------------
sink("multinom_habitat.txt")
cat("
model {
  for (i in 1:N) {
    yi[i,1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i,1:J] ~ ddirch(alphahab[hab[i], 1:J])
  }

  for (h in 1:H) {
    for (j in 1:J) {
      alphahab[h, j] <- q[h, j] * w[h] + 0.05
    }
    q[h, 1:J] ~ ddirch(alpha[])
    w[h] ~ dunif(0.1, 30)
  }

  for (j in 1:J) {
    alpha[j] <- 1
  }

  for (i in 1:N) {
    for (j in 1:J) {
      diff.prop[i,j] <- abs(pi[i,j] - q[hab[i], j])
    }
    PS[i] <- 1 - 0.5 * sum(diff.prop[i,1:J])
    log_lik[i] <- logdensity.multi(yi[i,], pi[i,], ni[i])
  }

  mean.PS <- mean(PS[])
}
", fill = TRUE)
sink()

inits_fun <- function() {
  q_init <- matrix(NA, nrow = jags_data_base$H, ncol = jags_data_base$J)
  for (h in 1:jags_data_base$H) {
    idx <- which(jags_data_base$hab == h)
    sub_y <- jags_data_base$yi[idx, , drop = FALSE]
    tot_ind <- apply(sub_y, 1, sum)
    q_obs <- apply(sub_y, 2, sum) / sum(tot_ind)
    q_obs[q_obs == 0] <- 1e-6
    q_obs <- q_obs / sum(q_obs)
    q_init[h, ] <- q_obs
  }
  list(q = q_init, w = rep(1, jags_data_base$H))
}

params  <- c("pi", "q", "w", "PS", "log_lik", "mean.PS")

mod_hab <- jags(
  data = jags_data_base,
  inits = inits_fun,
  parameters.to.save = params,
  model.file = "multinom_habitat.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 100,
  DIC = TRUE
)

## Habitat + Sex --------------------------
sink("multinom_habitat_sex.txt")
cat("
model {
  for (i in 1:N) {
    yi[i,1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i,1:J] ~ ddirch(alphasex[hab[i], sex[i], 1:J])
  }

  for (h in 1:H) {
    for (s in 1:S) {
      for (j in 1:J) {
        alphasex[h,s,j] <- q[h,s,j] * w[h,s] + 0.05
      }
      q[h,s,1:J] ~ ddirch(alpha[])
      w[h,s] ~ dunif(0.1, 30)
    }
  }

  for (j in 1:J) {
    alpha[j] <- 1
  }

  for (i in 1:N) {
    for (j in 1:J) {
      diff.prop[i,j] <- abs(pi[i,j] - q[hab[i], sex[i], j])
    }
    PS[i] <- 1 - 0.5 * sum(diff.prop[i,1:J])
    log_lik[i] <- logdensity.multi(yi[i,], pi[i,], ni[i])
  }

  mean.PS <- mean(PS[])
}
", fill = TRUE)
sink()

inits_fun2 <- function() {
  H <- jags_data_base$H
  S <- jags_data_base$S
  J <- jags_data_base$J
  
  q_init <- array(NA, dim = c(H, S, J))
  w_init <- array(1, dim = c(H, S))
  
  for (h in 1:H) {
    for (s in 1:S) {
      idx <- which(jags_data_base$hab == h & jags_data_base$sex == s)
      if (length(idx) > 0) {
        sub_y <- jags_data_base$yi[idx, , drop = FALSE]
        tot_ind <- apply(sub_y, 1, sum)
        q_obs <- apply(sub_y, 2, sum) / sum(tot_ind)
        q_obs[q_obs == 0] <- 1e-6
        q_obs <- q_obs / sum(q_obs)
      } else {
        q_obs <- rep(1/J, J)
      }
      q_init[h, s, ] <- q_obs
    }
  }
  
  list(q = q_init, w = w_init)
}

params2 <- c("pi", "q", "w", "PS", "log_lik", "mean.PS")

mod_hab_sex <- jags(
  data = jags_data_base,
  inits = inits_fun2,
  parameters.to.save = params2,
  model.file = "multinom_habitat_sex.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 100,
  DIC = TRUE
)

## Habitat + Sex + Standard length -------------------
sink("multinom_habitat_sex_size.txt")
cat("
model {
  for (i in 1:N) {
    yi[i,1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i,1:J] ~ ddirch(alphasize[hab[i], sex[i], size[i], 1:J])
  }


  for (h in 1:H) {
    for (s in 1:S) {
      for (k in 1:K) {
        for (j in 1:J) {
          alphasize[h,s,k,j] <- q[h,s,k,j] * w[h,s,k] + 0.05
        }
      
        q[h,s,k,1:J] ~ ddirch(alpha[])
      
        w[h,s,k] ~ dunif(0.1, 30)
      }
    }
  }

  # prior Dirichlet base
  for (j in 1:J) {
    alpha[j] <- 1
  }

  for (i in 1:N) {
    for (j in 1:J) {
      diff.prop[i,j] <- abs(pi[i,j] - q[hab[i], sex[i], size[i], j])
    }
    PS[i] <- 1 - 0.5 * sum(diff.prop[i,1:J])
    log_lik[i] <- logdensity.multi(yi[i,], pi[i,], ni[i])
  }

  mean.PS <- mean(PS[])
}
", fill = TRUE)
sink()

inits_fun3 <- function() {
  H <- jags_data_base$H
  S <- jags_data_base$S
  K <- jags_data_base$K
  J <- jags_data_base$J
  
  q_init <- array(NA, dim = c(H, S, K, J))
  w_init <- array(1, dim = c(H, S, K))
  
  for (h in 1:H) {
    for (s in 1:S) {
      for (k in 1:K) {
        idx <- which(jags_data_base$hab == h &
                       jags_data_base$sex == s &
                       jags_data_base$size == k)
        if (length(idx) > 0) {
          sub_y <- jags_data_base$yi[idx, , drop = FALSE]
          tot_ind <- apply(sub_y, 1, sum)
          q_obs <- apply(sub_y, 2, sum) / sum(tot_ind)
          q_obs[q_obs == 0] <- 1e-6
          q_obs <- q_obs / sum(q_obs)
        } else {
          
          q_obs <- rep(1/J, J)
        }
        q_init[h, s, k, ] <- q_obs
      }
    }
  }
  
  list(q = q_init, w = w_init)
}

params3 <- c("pi", "q", "w", "PS", "log_lik", "mean.PS")

mod_hab_sex_size <- jags(
  data = jags_data_base,
  inits = inits_fun3,
  parameters.to.save = params3,
  model.file = "multinom_habitat_sex_size.txt",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 100,
  DIC = TRUE
)

## Comparing models -----------------------------
loglik1 <- mod_hab$BUGSoutput$sims.list$log_lik
loglik2 <- mod_hab_sex$BUGSoutput$sims.list$log_lik
loglik3 <- mod_hab_sex_size$BUGSoutput$sims.list$log_lik

waic1 <- waic(loglik1)
waic2 <- waic(loglik2)
waic3 <- waic(loglik3)

waic1
waic2
waic3

loo_compare(waic1, waic2, waic3)

## PSi and IS values --------------------
sumj <- mod_hab$BUGSoutput$summary
idx  <- grep("^PS\\[", rownames(sumj))
PSibayes <- sumj[idx, "mean"]  
PSibayes_df <- data.frame(
  individuo = 1:length(PSibayes),
  PS = PSibayes,
  Env = data$Env
)

IS_by_hab <- PSibayes_df %>%
  group_by(Env) %>%
  summarise(IS = mean(PS, na.rm = TRUE))

IS_by_hab 

## Figure 4 --------------------------
df_PSi <- data.frame(
  PS = PSibayes,
  Env = data$Env   
)

df_PSi <- df_PSi %>%
  rename(Habitat = Env) %>%
  mutate(PS = as.numeric(PS))

df_PSi_summary <- df_PSi %>%
  group_by(Habitat) %>%
  summarise(
    mean_PS = mean(PS, na.rm = TRUE),
    sd_PS   = sd(PS, na.rm = TRUE),
    var_PS  = var(PS, na.rm = TRUE),
    n       = n(),
    .groups = "drop"
  )

df_PSi_summary

data$PS <- df_PSi$PS

data$Env <- factor(data$Env,
                            levels = rev(c("Stream", "Pool", "Ditch")))

fig4 <- data %>%
  ggplot(aes(x = PS, y = Env, fill = Env)) +
  ggdist::stat_halfeye(alpha = 0.8, adjust = 1.5, width = 0.7,
                       show.legend = FALSE) +
  scale_fill_manual(values = c("Stream" = "#009E73",
                               "Pool"   = "#0072B2",
                               "Ditch"  = "#E69F00")) +
  labs(x = expression(PS[i]), y = NULL) +
  theme_classic(base_size = 18)

fig4

ggsave("Figure_4.png", plot = fig4,
       width = 6, height = 5, dpi = 300, units = "in")


PS_stream <- data$PS[data$Env == "Stream"]
PS_pool   <- data$PS[data$Env == "Pool"]
PS_ditch  <- data$PS[data$Env == "Ditch"]

dip_stream <- dip.test(PS_stream)
dip_pool   <- dip.test(PS_pool)
dip_ditch  <- dip.test(PS_ditch)

dip_stream
dip_pool
dip_ditch

M_stream <- Mclust(PS_stream)
M_pool   <- Mclust(PS_pool)
M_ditch  <- Mclust(PS_ditch)

summary(M_stream)
summary(M_pool)
summary(M_ditch)

# Morphometrics X PSi -----------------------------
data$SL = as.numeric(data$SL)
data$Env <- fct_rev(data$Env)
data$SL
data$index_vental_flat
data$rel_eye_posit

stream_data <- filter(data, Env == "Stream")
mod_stream <- betareg(PS ~ SL + index_vental_flat + rel_eye_posit,
                      data = stream_data)
summary(mod_stream)
vif(mod_stream)

pool_data <- filter(data, Env == "Pool")
mod_pool <- betareg(PS ~ SL + index_vental_flat + rel_eye_posit,
                    data = pool_data)
summary(mod_pool)
vif(mod_pool)

ditch_data <- filter(data, Env == "Ditch")
mod_ditch <- betareg(PS ~ SL + index_vental_flat + rel_eye_posit,
                     data = ditch_data)
summary(mod_ditch)
vif(mod_ditch)

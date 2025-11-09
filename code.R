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

data_filtered$Env <- factor(data_filtered$Env, 
                            levels = c("Stream", "Pool", "Ditch"))

data_filtered <- data_filtered %>%
  dplyr::filter(Env == "Stream") %>%
  bind_rows(data %>% dplyr::filter(Env == "Pool")) %>%
  bind_rows(data %>% dplyr::filter(Env == "Ditch"))

data_filtered = 
  data_filtered %>% 
  filter(S > 0 )

data_filtered <- data_filtered |>
  mutate(
    SL = na_if(SL, "*"), 
    SL = as.numeric(SL)     
  )

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

# Frequency of occurrence ---------------------
data_filtered$Env <- factor(data_filtered$Env)
diet_total = decostand(data_filtered[,13:34], method = "pa")
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

ggsave("freq_occurrence.png", plot = fig2,
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.11))

ggsave("Levins.png", plot = fig3,
       width = 6, height = 4, dpi = 300, units = "in")

kruskal.test(Levins ~ Env, data = df_levins)
dunnTest(Levins ~ Env, data = df_levins, method = "bonferroni")

# Bayesian PSi values ---------------------------------------
## Stream IS Bayes -----------------------

sink("multinom.dir.hier.txt")
cat("
    data{
    for (i in 1:N) { 
    ni[i]  <- sum(yi[i,1:J]) # calcula o numero de presas de cada individuo
    }} model {
    for (i in 1:N) {
    yi[i, 1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i, 1:J] ~ ddirch(alphapop[1:J])
    }
    for (j in 1:J) {
    alphapop[j]  <- q[j] * w + 0.05 
    # os components de alphapop sao determinados pela dieta media da pop (q)
    # e o parametro de concentracao (w)
    }
    q[1:J] ~ ddirch(alpha[]) # # define uma distribuicao 
                             # uniforme a priori para q
    for (j in 1:J){
    alpha[j]  <- 1 }
    w ~ dunif(0.1, 30) # define o prior de w
    ### calculate PSi 
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pi[i, j] - q[j]) }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J]) }
    mean.PS  <- mean(PS[])
    for (i in 1:N){
    log_lik[i] <- logdensity.multi(yi[i, ], pi[i, ], ni[i]) }
    }
    ",fill = TRUE)
sink()

stream.Jags = list(
  yi = as.matrix(stream_diet), 
  N  = nrow(stream_diet),     
  J  = 22
)

Nprey.ind  = apply(stream_diet, 1, sum)            
q = as.numeric(apply(stream_diet, 2, 
                     function(x){sum(x)/sum(Nprey.ind)})) 

q_obs <- as.numeric(apply(stream_diet, 2, function(x) sum(x) / sum(Nprey.ind)))
q_obs[q_obs == 0] <- 1e-6                
q_obs <- q_obs / sum(q_obs)            

inits = function(){ list(q = q_obs, w = 1) }

params = c("pi", "q", "w", "PS") 

bayes.mod.stream = jags(stream.Jags, inits, params, model.file="multinom.dir.hier.txt",
                        n.chains=3, n.iter=2000, n.burnin=1000,
                        n.thin=100, DIC=TRUE, progress.bar = "text", digits=5)

N <- nrow(stream_diet)
sumj <- bayes.mod.stream$BUGSoutput$summary
idx  <- grep("^PS\\[", rownames(sumj))
PSibayes_stream <- sumj[idx, "mean"]

PSibayes.df_stream <- data.frame(
  individuo = 1:length(PSibayes_stream),
  PS = PSibayes_stream
)

View(PSibayes.df_stream)

ISbayes_stream = mean(PSibayes_stream)
ISbayes_stream # 0.6015

## Pool PSi Bayes -----------------------

sink("multinom.dir.hier.txt")
cat("
    data{
    for (i in 1:N) { 
    ni[i]  <- sum(yi[i,1:J]) # calcula o numero de presas de cada individuo
    }} model {
    for (i in 1:N) {
    yi[i, 1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i, 1:J] ~ ddirch(alphapop[1:J])
    }
    for (j in 1:J) {
    alphapop[j]  <- q[j] * w + 0.05 
    # os components de alphapop sao determinados pela dieta media da pop (q)
    # e o parametro de concentracao (w)
    }
    q[1:J] ~ ddirch(alpha[]) # # define uma distribuicao 
                             # uniforme a priori para q
    for (j in 1:J){
    alpha[j]  <- 1 }
    w ~ dunif(0.1, 30) # define o prior de w
    ### calculate PSi 
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pi[i, j] - q[j]) }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J]) }
    mean.PS  <- mean(PS[])
    for (i in 1:N){
    log_lik[i] <- logdensity.multi(yi[i, ], pi[i, ], ni[i]) }
    }
    ",fill = TRUE)
sink()

pool.Jags = list(
  yi = as.matrix(pool_diet), 
  N  = nrow(pool_diet),     
  J  = 22
)

Nprey.ind  = apply(pool_diet, 1, sum)            
q = as.numeric(apply(pool_diet, 2, 
                     function(x){sum(x)/sum(Nprey.ind)})) 

q_obs <- as.numeric(apply(pool_diet, 2, function(x) sum(x) / sum(Nprey.ind)))
q_obs[q_obs == 0] <- 1e-6                
q_obs <- q_obs / sum(q_obs)            

inits = function(){ list(q = q_obs, w = 1) }

params = c("pi", "q", "w", "PS") 

bayes.mod.pool = jags(pool.Jags, inits, params, model.file="multinom.dir.hier.txt",
                      n.chains=3, n.iter=2000, n.burnin=1000,
                      n.thin=100, DIC=TRUE, progress.bar = "text", digits=5)

N <- nrow(pool_diet)
sumj <- bayes.mod.pool$BUGSoutput$summary
idx  <- grep("^PS\\[", rownames(sumj))
PSibayes_pool <- sumj[idx, "mean"]

PSibayes.df_pool <- data.frame(
  individuo = 1:length(PSibayes_pool),
  PS = PSibayes_pool
)

View(PSibayes.df_pool)

ISbayes_pool = mean(PSibayes_pool)
ISbayes_pool # IS = 0.6078

## Ditch PSi Bayes --------------------------

sink("multinom.dir.hier.txt")
cat("
    data{
    for (i in 1:N) { 
    ni[i]  <- sum(yi[i,1:J]) # calcula o numero de presas de cada individuo
    }} model {
    for (i in 1:N) {
    yi[i, 1:J] ~ dmulti(pi[i,1:J], ni[i])
    pi[i, 1:J] ~ ddirch(alphapop[1:J])
    }
    for (j in 1:J) {
    alphapop[j]  <- q[j] * w + 0.05 
    # os components de alphapop sao determinados pela dieta media da pop (q)
    # e o parametro de concentracao (w)
    }
    q[1:J] ~ ddirch(alpha[]) # # define uma distribuicao 
                             # uniforme a priori para q
    for (j in 1:J){
    alpha[j]  <- 1 }
    w ~ dunif(0.1, 30) # define o prior de w
    ### calculate PSi 
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pi[i, j] - q[j]) }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J]) }
    mean.PS  <- mean(PS[])
    for (i in 1:N){
    log_lik[i] <- logdensity.multi(yi[i, ], pi[i, ], ni[i]) }
    }
    ",fill = TRUE)
sink()

ditch.Jags = list(
  yi = as.matrix(ditch_diet), 
  N  = nrow(ditch_diet),     
  J  = 22
)

Nprey.ind  = apply(ditch_diet, 1, sum)            
q = as.numeric(apply(ditch_diet, 2, 
                     function(x){sum(x)/sum(Nprey.ind)})) 

q_obs <- as.numeric(apply(ditch_diet, 2, function(x) sum(x) / sum(Nprey.ind)))
q_obs[q_obs == 0] <- 1e-6                
q_obs <- q_obs / sum(q_obs)            

inits = function(){ list(q = q_obs, w = 1) }

params = c("pi", "q", "w", "PS") 

bayes.mod.ditch = jags(ditch.Jags, inits, params, model.file="multinom.dir.hier.txt",
                       n.chains=3, n.iter=2000, n.burnin=1000,
                       n.thin=100, DIC=TRUE, progress.bar = "text", digits=5)

N <- nrow(ditch_diet)
sumj <- bayes.mod.ditch$BUGSoutput$summary
idx  <- grep("^PS\\[", rownames(sumj))
PSibayes_ditch<- sumj[idx, "mean"]

PSibayes.df_ditch <- data.frame(
  individuo = 1:length(PSibayes_ditch),
  PS = PSibayes_ditch
)

View(PSibayes.df_ditch)

ISbayes_ditch = mean(PSibayes_ditch)
ISbayes_ditch # IS = 0.5850

# Figure 4 -------------------------
PS_stream <- PSibayes.df_stream %>%
  mutate(Habitat = "Stream")
PS_pool <- PSibayes.df_pool %>%
  mutate(Habitat = "Pool")
PS_ditch <- PSibayes.df_ditch %>%
  mutate(Habitat = "Ditch")

df_PSi <- bind_rows(PS_stream, PS_pool, PS_ditch)

df_PSi <- df_PSi %>%
  mutate(PS = as.numeric(PS))

df_PSi_summary <- df_PSi %>%
  group_by(Habitat) %>%
  summarise(
    mean_PS = mean(PS, na.rm = TRUE),
    sd_PS   = sd(PS, na.rm = TRUE),
    var_PS  = var(PS, na.rm = TRUE),
    n       = n()
  )

df_PSi_summary

data_filtered$PS = df_PSi$PS

fig4 = data_filtered %>%
  ggplot(aes(x = PS, y = Env, fill = Env)) +
  ggdist::stat_halfeye(alpha = 0.8, adjust = 1.5, width = 0.7,
                       show.legend = F) +
  scale_fill_manual(values = c("Stream" = "#009E73", "Pool" = "#0072B2",
                               "Ditch" = "#E69F00")) +
  labs(x = expression(PS[i]), y = NULL) +
  theme_classic(base_size = 18)

fig4

ggsave("Figure_4.png", plot = fig4,
       width = 6, height = 5, dpi = 300, units = "in")

# Morphometrics X PSi -----------------------------
data_filtered$SL = as.numeric(data_filtered$SL)
data_filtered$Env <- fct_rev(data_filtered$Env)
data_filtered$SL
data_filtered$index_vental_flat
data_filtered$rel_eye_posit

stream_data <- filter(data_filtered, Env == "Stream")
mod_stream <- betareg(PS ~ SL + index_vental_flat + rel_eye_posit,
                      data = stream_data)
summary(mod_stream)
vif(mod_stream)

pool_data <- filter(data_filtered, Env == "Pool")
mod_pool <- betareg(PS ~ SL + index_vental_flat + rel_eye_posit,
                    data = pool_data)
summary(mod_pool)
vif(mod_pool)

ditch_data <- filter(data_filtered, Env == "Ditch")
mod_ditch <- betareg(PS ~ SL + index_vental_flat + rel_eye_posit,
                     data = ditch_data)
summary(mod_ditch)
vif(mod_ditch)

## Figure 5 -------------------------
mods <- list(
  Stream = mod_stream,
  Pool   = mod_pool,
  Ditch  = mod_ditch
)

coef_df <- imap_dfr(mods, ~{
  sm <- summary(.x)
  tab <- sm$coefficients$mean 
  as.data.frame(tab) |>
    mutate(term = rownames(tab),
           Habitat = .y)
})

coef_df <- coef_df |>
  filter(term != "(Intercept)") |>
  mutate(
    lower = Estimate - 1.96 * `Std. Error`,
    upper = Estimate + 1.96 * `Std. Error`,
    sig   = `Pr(>|z|)` < 0.05
  )

coef_df$term <- factor(coef_df$term,
                       levels = c("SL", "index_vental_flat", "rel_eye_posit"),
                       labels = c("SL", "Index ventral", "Rel. eye posit"))

xmax <- max(abs(c(coef_df$lower, coef_df$upper)), na.rm = TRUE)

fig5 = ggplot(coef_df,
       aes(x = Estimate, y = Habitat)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0, color = "black") +
  geom_point(aes(fill = sig),
             size = 3.8, shape = 21, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40")) +
  facet_wrap(
    ~ term,
    nrow = 1,
    scales = "fixed",
    labeller = as_labeller(c(
      "SL" = "SL",
      "Index ventral" = "IVF",
      "Rel. eye posit" = "REP"
    ))
  ) +
  coord_cartesian(xlim = c(-xmax, xmax)) +
  labs(
    x = "Coefficient estimate (logit link)",
    y = NULL
  )+
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )

fig5

ggsave("Figure_5.png", plot = fig5,
       width = 8, height = 4, dpi = 300, units = "in")

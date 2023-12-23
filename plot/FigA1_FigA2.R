rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source("s.R")
c4 <- rev(brewer.pal(4, "Set1"))
if(!dir.exists(paste(arg[1], '/supp/', sep = ''))) {dir.create(paste(arg[1], '/supp/', sep = ''))}

################################################################################ LOAD RESULTS
S <- 5
L <- 20
C <- c('ML', 'FL', 'BN')
res <- data.frame()
for(u in C) {
  df <- read.table(paste('../model/res/', u, '.txt', sep = ''), sep = '\t', header = T)
  df$uni_var <- NA
  df$uni_var[which(df$nfa == S & df$nma == S & df$nfb == S & df$nmb == S)] <- 'uni_s'
  df$uni_var[which(df$nfa == S & df$nma == S & df$nfb == L & df$nmb == L)] <- 'var_sl'
  df$uni_var[which(df$nfa == L & df$nma == L & df$nfb == L & df$nmb == L)] <- 'uni_l'
  df$case <- u
  res <- rbind(res, df)
}

res$uni_var[which(res$uni_var == 'var_sl' & res$patch_type == 'ALPHA')] <- 'var_s'
res$uni_var[which(res$uni_var == 'var_sl' & res$patch_type == 'BETA')] <- 'var_l'

res$r <- NA
for(v in 1:nrow(res)) {
  if(res$sexes[v] == 'F' & res$patch_type[v] == 'ALPHA') {
    res$r[v] <- ((res$nfa[v] - 1) * res$r_ff[v] + res$nma[v] * res$r_mf[v]) / (res$nfa[v] + res$nma[v] - 1)
  } else if(res$sexes[v] == 'F' & res$patch_type[v] == 'BETA') {
    res$r[v] <- ((res$nfb[v] - 1) * res$r_ff[v] + res$nmb[v] * res$r_mf[v]) / (res$nfb[v] + res$nmb[v] - 1)
  } else if(res$sexes[v] == 'M' & res$patch_type[v] == 'ALPHA') {
    res$r[v] <- ((res$nma[v] - 1) * res$r_mm[v] + res$nfa[v] * res$r_fm[v]) / (res$nfa[v] + res$nma[v] - 1)
  } else {
    res$r[v] <- ((res$nmb[v] - 1) * res$r_mm[v] + res$nfb[v] * res$r_fm[v]) / (res$nfb[v] + res$nmb[v] - 1)
  }
}

################################################################################ KINSHIP DYNAMICS
df_0 <- res[which(res$uni_var == 'uni_s' & res$fec_or_mor == 'FECUNDITY'), ]
df_0 <- rbind(df_0, res[which(res$uni_var == 'uni_l' & res$fec_or_mor == 'FECUNDITY'), ])
df_1 <- res[which(res$uni_var == 'var_s' & res$fec_or_mor == 'FECUNDITY'), ]
df_1 <- rbind(df_1, res[which(res$uni_var == 'var_l' & res$fec_or_mor == 'FECUNDITY'), ])
d <- rbind(df_0, df_1)

# complement data for the counterpart homogeneous population
temp <- res[which(res$uni_var == 'uni_s' & res$fec_or_mor == 'FECUNDITY'), ]
temp$u <- .1
d <- rbind(d, temp)
temp <- res[which(res$uni_var == 'uni_s' & res$fec_or_mor == 'FECUNDITY'), ]
temp$u <- .9
d <- rbind(d, temp)
#
temp <- res[which(res$uni_var == 'uni_l' & res$fec_or_mor == 'FECUNDITY'), ]
temp$u <- .1
d <- rbind(d, temp)
temp <- res[which(res$uni_var == 'uni_l' & res$fec_or_mor == 'FECUNDITY'), ]
temp$u <- .9
d <- rbind(d, temp)

d$sexes[which(d$sexes == 'F')] <- 'females'
d$sexes[which(d$sexes == 'M')] <- 'males'
d$case[which(d$case == 'ML')] <- 'typical mammals'
d$case[which(d$case == 'FL')] <- 'apes'
d$case[which(d$case == 'BN')] <- 'whales'
d$case <- factor(d$case, levels = c('whales', 'typical mammals', 'apes'))
d$u[which(d$u == .1)] <- 'low(10%)'
d$u[which(d$u == .5)] <- 'medium(50%)'
d$u[which(d$u == .9)] <- 'high(90%)'
d$u <- factor(d$u, levels = c('low(10%)', 'medium(50%)', 'high(90%)'))
d$uni_var <- factor(d$uni_var, levels = c('uni_s', 'uni_l', 'var_s', 'var_l'))

FA1 <- ggplot(d) + 
  geom_line(aes(x = age * .1, y = r, group = interaction(uni_var, case, sexes, u), color = uni_var, linetype = sexes), linewidth = 1.25) +
  labs(x = expression(paste('Age (relative to population mean lifespan)')), y = expression(paste(bar(italic('r'))))) +
  scale_y_continuous(breaks = c(0, .06, .12)) +
  scale_colour_manual(values = c4[c(3, 4, 1, 2)]) +
  facet_grid(u ~ case) +
  the
FA1 <- ggdraw(FA1) + draw_label(expression(bold(paste('\u03bc'[f], '=', '\u03bc'[m], '=0.1', sep = ''))), x = .1, y = .95, size = 17)
ggsave(plot = FA1, './supp/Fig_A1.pdf', width = 22, height = 16, units = 'cm', device = cairo_pdf)

################################################################################ SELECTIVE PRESSURES
# ON BEHAVIOURS WITH FECUNDITY OUTCOMES
FA2A <- ggplot(d) + 
  geom_line(aes(x = age * .1, y = c / critical_b, group = interaction(uni_var, case, sexes, u), color = uni_var, linetype = sexes), linewidth = 1.25) +
  labs(x = expression(paste('Age (relative to population mean lifespan)')), y = expression(paste(italic('c/b'^''['*']), ' (fecundity)')), tag = '(A)') + 
  scale_y_continuous(breaks = c(0, .06, .12)) +
  scale_colour_manual(values = c4[c(3, 4, 1, 2)]) +
  facet_grid(u ~ case) +
  geom_hline(yintercept = 0, linetype = 'solid', color = "grey0", linewidth = .2) + 
  the +
  theme(plot.tag.position = c(.02, .95), text = element_text(size = 20, face = "bold"))
FA2A <- ggdraw(FA2A) + draw_label(expression(bold(paste('\u03bc'[f], '=', '\u03bc'[m], '=0.1', sep = ''))), x = .1, y = .95, size = 17)

################################################################################ SELECTIVE PRESSURES
# ON BEHAVIOURS WITH MORTALITY OUTCOMES
df_0 <- res[which(res$uni_var == 'uni_s' & res$fec_or_mor == 'MORTALITY'), ]
df_0 <- rbind(df_0, res[which(res$uni_var == 'uni_l' & res$fec_or_mor == 'MORTALITY'), ])
df_1 <- res[which(res$uni_var == 'var_s' & res$fec_or_mor == 'MORTALITY'), ]
df_1 <- rbind(df_1, res[which(res$uni_var == 'var_l' & res$fec_or_mor == 'MORTALITY'), ])
d <- rbind(df_0, df_1)

# complement data for the counterpart homogeneous population
temp <- res[which(res$uni_var == 'uni_s' & res$fec_or_mor == 'MORTALITY'), ]
temp$u <- .1
d <- rbind(d, temp)
temp <- res[which(res$uni_var == 'uni_s' & res$fec_or_mor == 'MORTALITY'), ]
temp$u <- .9
d <- rbind(d, temp)
#
temp <- res[which(res$uni_var == 'uni_l' & res$fec_or_mor == 'MORTALITY'), ]
temp$u <- .1
d <- rbind(d, temp)
temp <- res[which(res$uni_var == 'uni_l' & res$fec_or_mor == 'MORTALITY'), ]
temp$u <- .9
d <- rbind(d, temp)

d$sexes[which(d$sexes == 'F')] <- 'females'
d$sexes[which(d$sexes == 'M')] <- 'males'
d$case[which(d$case == 'ML')] <- 'typical mammals'
d$case[which(d$case == 'FL')] <- 'apes'
d$case[which(d$case == 'BN')] <- 'whales'
d$case <- factor(d$case, levels = c('whales', 'typical mammals', 'apes'))
d$u[which(d$u == .1)] <- 'low(10%)'
d$u[which(d$u == .5)] <- 'medium(50%)'
d$u[which(d$u == .9)] <- 'high(90%)'
d$u <- factor(d$u, levels = c('low(10%)', 'medium(50%)', 'high(90%)'))
d$uni_var <- factor(d$uni_var, levels = c('uni_s', 'uni_l', 'var_s', 'var_l'))

FA2B <- ggplot(d) + 
  geom_line(aes(x = age * .1, y = c / critical_b, group = interaction(uni_var, case, sexes, u), color = uni_var, linetype = sexes), linewidth = 1.25) + 
  labs(x = expression(paste('Age (relative to population mean lifespan)')), y = expression(paste(italic('c/b'^''['*']), ' (mortality)')), tag = '(B)') + 
  scale_y_continuous(breaks = c(-0.04, 0, .04)) +
  scale_colour_manual(values = c4[c(3, 4, 1, 2)]) +
  facet_grid(u ~ case) +
  geom_hline(yintercept = 0, linetype = 'solid', color = "grey0", linewidth = .2) +
  the +
  theme(plot.tag.position = c(.02, .95), text = element_text(size = 20, face = "bold"))
FA2B <- ggdraw(FA2B) + draw_label(expression(bold(paste('\u03bc'[f], '=', '\u03bc'[m], '=0.1', sep = ''))), x = .1, y = .95, size = 17)
twoplots <- align_plots(FA2A, FA2B, align = "hv", axis = "tblr")

ggdraw(twoplots[[1]])
ggsave('./supp/Fig_A2A.pdf', width = 24, height = 16, units = 'cm', device = cairo_pdf)
ggdraw(twoplots[[2]])
ggsave('./supp/Fig_A2B.pdf', width = 24, height = 16, units = 'cm', device = cairo_pdf)

if(file.exists("./Rplots.pdf")) {file.remove("./Rplots.pdf")}
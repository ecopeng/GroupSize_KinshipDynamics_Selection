rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source("s.R")
if(!dir.exists(paste(arg[1], '/main/', sep = ''))) {dir.create(paste(arg[1], '/main/', sep = ''))}

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
small <- paste('small(', paste(S, '+', S, sep = ''), ')', sep = '')
large <- paste('large(', paste(L, '+', L, sep = ''), ')', sep = '')
df_1 <- res[which(res$uni_var == 'var_s' & res$fec_or_mor == 'FECUNDITY'), ]
df_2 <- rbind(df_1, res[which(res$uni_var == 'var_l' & res$fec_or_mor == 'FECUNDITY'), ])
d <- rbind(df_1, df_2)
d$uni_var[which(d$uni_var == 'var_s')] <- small
d$uni_var[which(d$uni_var == 'var_l')] <- large
d$uni_var <- factor(d$uni_var, levels = c(small, large))
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
F1 <- ggplot(d) + 
  geom_line(aes(x = age * .1, y = r, group = interaction(uni_var, case, sexes, u), color = sexes, linetype = uni_var), linewidth = 1) + 
  labs(x = expression(paste('Age (relative to population mean lifespan)')), y = expression(paste(bar(italic('r'))))) + 
  scale_y_continuous(breaks = c(.02, .06, .1)) +
  scale_colour_manual(values = c4[2:1]) +
  facet_grid(u ~ case) +
  the + 
  theme(text = element_text(size = 20, face = "bold"))
F1 <- ggdraw(F1) + draw_label(expression(bold(paste('\u03bc'[f], '=', '\u03bc'[m], '=0.1', sep = ''))), x = .1, y = .95, size = 17)
ggsave(plot = F1, './main/Fig_1.pdf', width = 22, height = 16, units = 'cm', device = cairo_pdf)

################################################################################ SELECTIVE PRESSURES
# ON BEHAVIOURS WITH FECUNDITY OUTCOMES
x_seq <- seq(.1, 3, by = 1e-05)
yic <- 0
intersections <- d %>%
  group_by(uni_var, case, sexes, u) %>%
  summarise(interpolated = approx(x = age * .1, y = c / critical_b, xout = x_seq)$y) %>%
  mutate(x_seq = x_seq) %>%
  slice_min(abs(interpolated - yic))
summary(abs(intersections$interpolated))
cp <- intersections[which(abs(intersections$interpolated) <= 1e-05), ]
ycoods <- rep(.06, nrow(cp))
F2A <- ggplot(d) + 
  geom_hline(yintercept = 0, linetype = 'dotted', color = "black", linewidth = .5) +
  geom_vline(data = cp, aes(xintercept = x_seq), linetype = 'dotted', color = "black", linewidth = .5) + 
  geom_richtext(data = cp, mapping = aes(x = x_seq, y = ycoods, label = sub("\\(.*", "", uni_var)), inherit.aes = F, vjust = .5, angle = 90, size = 4.5, label.size = NA) + 
  geom_line(aes(x = age * .1, y = c / critical_b, group = interaction(uni_var, case, sexes, u), color = sexes, linetype = uni_var), linewidth = 1) + 
  labs(x = expression(paste('Age (relative to population mean lifespan)')), y = expression(paste(italic('c/b'^''['*']), ' (fecundity)')), tag = '(A)') + 
  scale_y_continuous(breaks = c(0, .05, .1)) +
  scale_colour_manual(values = c4[2:1]) +
  facet_grid(u ~ case) +
  the +
  theme(plot.tag.position = c(.02, .95), text = element_text(size = 20, face = "bold")) +
  geom_point(data = cp, aes(x_seq, interpolated), size = 3.5, color = 'black', shape = 20)
F2A <- ggdraw(F2A) + draw_label(expression(bold(paste('\u03bc'[f], '=', '\u03bc'[m], '=0.1', sep = ''))), x = .08, y = .95, size = 17)

################################################################################ SELECTIVE PRESSURES
# ON BEHAVIOURS WITH MORTALITY OUTCOMES
df_1 <- res[which(res$uni_var == 'var_s' & res$fec_or_mor == 'MORTALITY'), ]
df_2 <- rbind(df_1, res[which(res$uni_var == 'var_l' & res$fec_or_mor == 'MORTALITY'), ])
d <- rbind(df_1, df_2)
d$uni_var[which(d$uni_var == 'var_s')] <- small
d$uni_var[which(d$uni_var == 'var_l')] <- large
d$uni_var <- factor(d$uni_var, levels = c(small, large))
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
x_seq <- seq(.1, 3, by = 1e-05)
yic <- 0
intersections <- d %>%
  group_by(uni_var, case, sexes, u) %>%
  summarise(interpolated = approx(x = age * .1, y = c / critical_b, xout = x_seq)$y) %>%
  mutate(x_seq = x_seq) %>%
  slice_min(abs(interpolated - yic))
cp <- intersections[which(abs(intersections$interpolated) <= 1e-05), ]
F2B <- ggplot(d) + 
  geom_hline(yintercept = 0, linetype = 'dotted', color = "black", linewidth = .5) +
  geom_vline(data = cp, aes(xintercept = x_seq), linetype = 'dotted', color = "black", linewidth = .5) + 
  geom_line(aes(x = age * .1, y = c / critical_b, group = interaction(uni_var, case, sexes, u), color = sexes, linetype = uni_var), linewidth = 1) +
  labs(x = expression(paste('Age (relative to population mean lifespan)')), y = expression(paste(italic('c/b'^''['*']), ' (mortality)')), tag = '(B)') +
  scale_y_continuous(breaks = c(-.04, 0, .04)) +
  scale_colour_manual(values = c4[2:1]) +
  facet_grid(u ~ case) +
  the +
  theme(plot.tag.position = c(.02, .95), text = element_text(size = 20, face = "bold")) +
  geom_point(data = cp, aes(x_seq, interpolated), size = 3, color = 'black', shape = 21, stroke = .85) +
  geom_point(data = cp, aes(x_seq, interpolated), size = 1, color = 'black', shape = 20)
F2B <- ggdraw(F2B) + draw_label(expression(bold(paste('\u03bc'[f], '=', '\u03bc'[m], '=0.1', sep = ''))), x = .08, y = .95, size = 17)
twoplots <- align_plots(F2A, F2B, align = "hv", axis = "tblr")
ggsave('./main/Fig_2A.pdf', plot = ggdraw(twoplots[[1]]), width = 22, height = 16, units = 'cm', device = cairo_pdf)
ggsave('./main/Fig_2B.pdf', plot = ggdraw(twoplots[[2]]), width = 22, height = 16, units = 'cm', device = cairo_pdf)

################################################################################ BIRTH RATE FEMALE KILLER WHALES
rm(list = ls())
source('s.R')

female_data <- read.csv('../data/births.csv', stringsAsFactors = F)
model <- glm(birth ~ age_class * pod, data = female_data, family = binomial)
pred_data <- expand.grid(age_class = c("young", "old"), pod = c("K","J","L"))
pod_size <- read.csv('../data/pod_size.csv', stringsAsFactors = F)
birth_rate <- predict(model, newdata = pred_data, se.fit = T, type = 'response')
pred_data$birth_rate <- birth_rate$fit
pred_data$se_rate <- birth_rate$se.fit

pod <- reshape2::melt(pod_size)
p_pod_size <- ggplot(data = pod, aes(x = variable, y = value, group = variable)) + 
  stat_boxplot(geom = 'errorbar', width = .25) +
  labs(x = 'Pod', y = 'Size') +
  geom_boxplot(fill = 'grey80', color = 'black', outlier.size = .75, fatten = 4, width = .5) + 
  the +
  scale_y_continuous(breaks = c(15, 30, 45)) +
  theme(plot.tag.position = c(.5, .925),
        text = element_text(size = 25, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 30, margin = margin(0, 5, 0, 1, 'pt')),
        axis.ticks.y = element_line(color = "black", linewidth = (1)),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(0, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, margin = margin(0, 5, 0, 1, 'pt')))

p_birth_rate <- ggplot(pred_data, aes(x = pod, y = birth_rate, color = age_class)) +
  geom_errorbar(aes(ymin = birth_rate - se_rate, ymax = birth_rate + se_rate), width = 0, position = position_dodge(.25), linewidth = 2, show.legend = T) +
  geom_point(size = 8, shape = 16, fill = 'white', position = position_dodge(.25)) +
  the + 
  labs(y = 'Birth rate', x = 'Pod') +
  scale_y_continuous(breaks = c(.05, .1, .15)) +
  theme(plot.tag.position = c(.25, .925),
        text = element_text(size = 25, face = "bold"),
        legend.direction = "vertical",
        legend.position = c(.875, .875),
        legend.text = element_text(size = 16, margin = margin(0, 8, 0, 0, 'pt')),
        legend.key.size = unit(.1, "lines"),
        axis.text.x = element_text(size = 30, margin = margin(5, 0, 5, 0, 'pt')),
        axis.text.y = element_text(size = 30, margin = margin(0, 5, 0, 1, 'pt')),
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks = element_line(color = "black", linewidth = (1)),
        axis.title.x = element_text(size = 40, margin = margin(5, 0, 5, 0, 'pt')),
        axis.title.y = element_text(size = 40, margin = margin(0, 5, 0, 1, 'pt')),
        panel.border = element_rect(linewidth = 1),
        plot.margin = unit(c(3, .5, -2.5, 1), "lines"))
F3 <- egg::ggarrange(tag_facet(p_pod_size, tag_pool = "A", size = 8), tag_facet(p_birth_rate, tag_pool = "B", size = 8), ncol = 1)
ggsave(filename = './main/Fig_3.pdf', plot = F3, width = 22, height = 25, units = 'cm', device = cairo_pdf)

if(file.exists("./Rplots.pdf")) {file.remove("./Rplots.pdf")}
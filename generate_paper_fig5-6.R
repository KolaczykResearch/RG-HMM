library('ggplot2')
library('gridExtra')

#####################################################
# FIGURE 5
#####################################################
# ---------------
# EXP2_vary_B - VERSION2 - all together
# ---------------
data = read.csv('Simulations/Estimation/EXP2_vary_B/data_exp2.csv')
data$N = as.factor(data$N)
data$M = as.factor(data$M)
data$B = as.factor(data$B)
# remove B = 150000
data = data[data$B!=150000,] 
pd = position_dodge(0)

# ====== p, q
data_new = data[c('B', 'mean_p', 'se_p')]
colnames(data_new) = c('B', 'mean', 'se')
data_new2 = data[c('B', 'mean_q', 'se_q')]
colnames(data_new2) = c('B', 'mean', 'se')
data_new = rbind(data_new, data_new2)
data_new['Parameter'] = rep(c('p', 'q'), each=4)
# data_new = data_new %>% 
#   mutate(Parameter = recode_factor(Parameter, 'p' = 'alpha', 'q' = 'beta'))
plt1 = ggplot(data_new, aes(x=B, y=mean, group=Parameter)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape=Parameter), size=4) + 
  geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.3, linetype="dashed", color = "red") + 
  scale_shape_manual(values=c(4,6), labels = parse_format()) + ylim(0.1, 0.8) +
  labs(y='Estimates') + theme_bw()
plt1

# ====== gamma
data_new = data[c('B', 'mean_gamma', 'se_gamma')]
data_new['Parameter'] = rep(c('gamma'), 4)
plt3 = ggplot(data_new, aes(x=B, y=mean_gamma, group=Parameter)) +
  geom_errorbar(aes(ymin=mean_gamma-se_gamma, ymax=mean_gamma+se_gamma), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape=Parameter), size=4) + 
  scale_shape_manual(values=c(7), labels = parse_format()) +
  geom_hline(yintercept=2, linetype="dashed", color = "red") + 
  labs(y='Estimates') + theme_bw()
plt3

# =======alpha
data_new = data[c('B', 'mean_alpha', 'se_alpha')]
data_new['Parameter'] = rep(c('alpha'), 4)
plt4 = ggplot(data_new, aes(x=B, y=mean_alpha, group=Parameter)) +
  geom_errorbar(aes(ymin=mean_alpha-se_alpha, ymax=mean_alpha+se_alpha), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape=Parameter), size=3) + 
  scale_shape_manual(values=c(8), labels = parse_format()) + ylim(0.03, NA) + 
  geom_hline(yintercept=0.03, linetype="dashed", color = "red") + 
  labs(y='Estimates') + theme_bw()
plt4

# ===== beta
data_new = data[c('B', 'mean_beta', 'se_beta')]
data_new['Parameter'] = rep(c('beta'), 4)
plt5 = ggplot(data_new, aes(x=B, y=mean_beta, group=Parameter)) +
  geom_errorbar(aes(ymin=mean_beta-se_beta, ymax=mean_beta+se_beta), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(aes(shape=Parameter), size=3) + 
  scale_shape_manual(values=c(9), labels = parse_format()) + ylim(0.01, 0.11) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "red") + 
  labs(y='Estimates') + theme_bw()
  # labs(y=expression('Type-II Error Rate'~beta)) + theme_bw()
plt5

# 10.13*5.18
grid.arrange(plt1, plt3, plt4, plt5, layout_matrix = 
               cbind(rep(1, 14), c(2,2,2,2,2,2,2,2,3,3,3,4,4,4)))

#####################################################
# FIGURE 6
#####################################################
# ---------------
# EXP3_vary_NM - VERSION2 - all together
# ---------------
data = read.csv('Simulations/Estimation/EXP3_vary_NM/data_exp3.csv')
data$N = as.factor(data$N)
data$M = as.factor(data$M)
data = data[data$N!=100,]
pd = position_dodge(0.05)

# --- for p, q
data_new = data[c('M', 'N', 'mean_p', 'se_p')]
colnames(data_new) = c('M', 'N', 'mean', 'se')
data_new2 = data[c('M', 'N', 'mean_q', 'se_q')]
colnames(data_new2) = c('M', 'N', 'mean', 'se')
data_new = rbind(data_new, data_new2)
data_new['Parameter'] = rep(c('p', 'q'), each=12)
plt2_1 = ggplot(data_new, aes(x=M, y=mean, group=interaction(N, Parameter), 
                              color=N, shape=Parameter, linetype=N)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(size=4, position=pd) + 
  geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.3, linetype="dashed", color = "red") + 
  scale_shape_manual(values=c(4,6), labels = parse_format()) + ylim(0.1, 0.8) +
  labs(y='Estimates') + theme_bw() + 
  guides(color = guide_legend(override.aes = list(size = 1))) + scale_linetype_manual(values=c(1,2,10))
plt2_1

# ====== gamma
data_new = data[c('M', 'N', 'mean_gamma', 'se_gamma')]
colnames(data_new) = c('M', 'N', 'mean', 'se')
data_new['Parameter'] = rep(c('gamma'), 12)
plt2_3 = ggplot(data_new, aes(x=M, y=mean, group=interaction(Parameter, N), 
                              color=N, shape=Parameter, linetype=N)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(size=4, position=pd) + 
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  scale_shape_manual(values=c(7), labels = parse_format()) +
  labs(y='Estimates') + theme_bw() + 
  scale_linetype_manual(values=c(1,2,10)) + 
  guides(color = FALSE, linetype=FALSE)
plt2_3

# ====== alpha
data_new = data[c('M', 'N', 'mean_alpha', 'se_alpha')]
colnames(data_new) = c('M', 'N', 'mean', 'se')
data_new['Parameter'] = rep(c('alpha'), 12)
plt2_4 = ggplot(data_new, aes(x=M, y=mean, group=interaction(Parameter, N), 
                              color=N, shape=Parameter, linetype=N)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(size=4, position=pd) + 
  geom_hline(yintercept=0.03, linetype="dashed", color = "red") +
  scale_shape_manual(values=c(8), labels = parse_format()) +
  scale_linetype_manual(values=c(1,2,10)) + 
  labs(y='Estimates') + theme_bw() + 
  guides(color = FALSE, linetype=FALSE)
plt2_4

# ====== beta
data_new = data[c('M', 'N', 'mean_beta', 'se_beta')]
colnames(data_new) = c('M', 'N', 'mean', 'se')
data_new['Parameter'] = rep(c('beta'), 12)
plt2_5 = ggplot(data_new, aes(x=M, y=mean, group=interaction(Parameter, N), 
                              color=N, shape=Parameter, linetype=N)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(size=4, position=pd) + 
  geom_hline(yintercept=0.01, linetype="dashed", color = "red") +
  scale_shape_manual(values=c(9), labels = parse_format()) +
  labs(y='Estimates') + theme_bw() + 
  scale_linetype_manual(values=c(1,2,10)) + 
  guides(color = FALSE, linetype=FALSE)
plt2_5

# 10.13*5.18
grid.arrange(plt2_1, plt2_3, plt2_4, plt2_5, layout_matrix = 
               cbind(rep(1, 14), c(2,2,2,2,2,2,2,2,3,3,3,4,4,4)))





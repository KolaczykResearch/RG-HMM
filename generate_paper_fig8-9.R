library('ggplot2')
# ------------------------------------------
# -----STUDY 1------------------------------
# ------------------------------------------
######################################
# Figure 8 (left)
#####################################
data10_p = read.csv('Simulations/Testing/N_10/data10.csv', header=T)
data20_p = read.csv('Simulations/Testing/N_20/data20.csv', header=T)
data30_p = read.csv('Simulations/Testing/N_30/data30.csv', header=T)

data = rbind(data10_p, data20_p)
data = rbind(data, data30_p)
data$scaled_t = data$scaled_t_obs*2

data$B = as.factor(data$B)
data$scaled_t = as.factor(data$scaled_t)
data$N = as.factor(data$N)

plt = ggplot(data, aes(x=scaled_t, y=Rate_of_Detection, group=B, color=B,
shape=B)) +
geom_errorbar(aes(ymin=Rate_of_Detection-se, ymax=Rate_of_Detection+se),
width=.1, position=position_dodge(width=0.4)) +
geom_point(position=position_dodge(width=0.4), size=2) +
labs(y='Rate of Detection') +
scale_shape_manual(values=c(1,2,3)) + theme_bw() +
scale_color_manual(values=c(rgb(227, 60, 61, maxColorValue = 255),
rgb(61, 127, 182, maxColorValue = 255), rgb(81, 173, 79, maxColorValue = 255))) +
theme(text = element_text(size = 14))
plt = plt + facet_grid(N ~ process)
plt
# 6.89*6.99

######################################
# Figure 8 (right), Figure 9
#####################################
data10 = read.csv('Simulations/Testing/N_10/data10_all.csv', header=T)
data20 = read.csv('Simulations/Testing/N_20/data20_all.csv', header=T)
data30 = read.csv('Simulations/Testing/N_30/data30_all.csv', header=T)
data20$seed = data20$seed -1000
data30$seed = data30$seed -100

data_all = rbind(data10, data20)
data_all = rbind(data_all, data30)
data_all$scaled_t = data_all$scaled_t_obs*2

data_all$B = as.factor(data_all$B)
data_all$seed = as.factor(data_all$seed)
data_all$N = as.factor(data_all$N)
data_all$scaled_t = as.factor(data_all$scaled_t)

model = glmer(success ~ N + scaled_t + process + B + (1|seed), data=data_all, family=binomial(logit))
summary(model)

# Figure 8 (right)
# theme_set(theme_sjplot())
set_theme(
axis.textcolor = "black",
base = theme_bw()
)
plot_model(model, type = "pred", terms = c("scaled_t", 'B', 'process'),
title='Fitted Probabilities of Successful Detection from GLMM',
axis.title = 'Probability of Successful Detection',
ci.lvl=0.95, dodge=0.3) +
aes(shape = group) + scale_shape_manual(values = c(1, 2, 3)) +
labs(shape = "B")
# 6.89*6.99

# Figure 9
# ---- coefficients
plot_model(model, vline.color = "grey", order.terms = c(3, 4, 5, 7, 8, 6, 1, 2),
show.values = TRUE, value.offset = .4, title='') +
theme(text = element_text(size = 14))

######################################
# Table 4
#####################################
Anova(model)
xtable(Anova(model))

# ------------------------------------------
# -----STUDY 2------------------------------
# ------------------------------------------
###############
# 2.1 Summary plots (mean and estimated sd) for rate of detection
###############
data_low_r = read.csv('N_10_20_30_vary_r/data10_20_30_low_r.csv', header=T)
data_low_r$scaled_t = data_low_r$scaled_t_obs*2
data_low_r$N = as.factor(data_low_r$N)
data_low_r$rate_obs = as.factor(data_low_r$rate_obs)

plt = ggplot(data_low_r, aes(x=rate_obs, y=Rate_of_Detection, group=process, color=process)) +
geom_errorbar(aes(ymin=Rate_of_Detection-se, ymax=Rate_of_Detection+se),
width=.1, position=position_dodge(width=0.2)) +
geom_point(position=position_dodge(width=0.2)) +
labs(y='Rate of Detection', x='observation rate')
plt = plt + facet_wrap(~N, ncol=3) + theme(text = element_text(size = 14))
plt
# 6.89*6.99

###############
# 2.2 Fit GLMM
###############
data_low_r = read.csv('N_10_20_30_vary_r/data10_20_30_low_r_all.csv', header=T)
data_low_r[data_low_r$N==20, ]$seed = data_low_r[data_low_r$N==20, ]$seed  -1000
data_low_r[data_low_r$N==30, ]$seed = data_low_r[data_low_r$N==30, ]$seed  -100
data_all = data_low_r # should also consider rate_obs=1.5 from other three folders
#
data_all$B = as.factor(data_all$B)
data_all$scaled_t_obs = as.factor(data_all$scaled_t_obs)
data_all$seed = as.factor(data_all$seed)
data_all$N = as.factor(data_all$N)
data_all$rate_obs = as.factor(data_all$rate_obs)
data_all$observation_rate = data_all$rate_obs

model2 = glmer(success ~ N + process+observation_rate+ (1|seed), data=data_all, family=binomial(logit))
summary(model2)

Anova(model2)
xtable(Anova(model2))

plot_model(model2, type = "pred", terms = c("observation_rate", 'process'),
title='Fitted Probabilities of Successful Detection from GLMM', axis.title = 'Probability of Successful Detection')
# 6.89*6.99
plot_model(model2, vline.color = "grey",
show.values = TRUE, value.offset = .4, title='')





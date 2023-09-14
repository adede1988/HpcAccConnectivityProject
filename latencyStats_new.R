library(tidyverse)
library(ggpubr)
library(rstatix)
library(performance)
library(lmerTest)

#standard power comparison across conditions in different times and frequencies independently ####
data = read.csv('R:\\MSS\\Johnson_Lab\\dtf8829\\latencyAOVdat.csv')
data = read.csv('R:\\MSS\\Johnson_Lab\\dtf8829\\latencyAOVdatNEW.csv')
data <- data %>% filter(cond == 'subMiss' | cond == 'subHit' | cond=='hit_on' | cond=='miss_on')
data$encRet[data$cond == 'subMiss' | data$cond == 'subHit'] <- 'enc'
data$encRet[data$cond=='hit_on' | data$cond=='miss_on'] <- 'ret'
data$hitMiss = 'hit'
data$hitMiss[data$cond=='miss_on' | data$cond == 'subMiss'] = 'miss'
data$realID = paste(data$subID, data$chi)
data$adjTime = data$centerOfMass / data$RT

data <- data %>% filter(RT<3000)
# data <- data %>% filter(reg == 'hip' | reg == 'phg')

#regress out RT first



data.rt = lm(centerOfMass ~ RT, data = data)
data$rt_res = data.rt$residuals

data.lmer <- lmer(rt_res ~ reg*hitMiss*encRet + (1 | chi/subID), REML = TRUE, 
                  control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE), data = data)

Anova(data.lmer)
# summary(data.lmer)
performance(data.lmer)
performance(data.rt)


#show that there's a lot of variance with RT
data %>% 
  ggplot(aes(x=RT, y=centerOfMass, color = reg)) + 
  geom_point()
  



library(ggeffects)

model_cond <- ggemmeans(data.lmer,
                        terms = c("reg","hitMiss", "encRet"),
                        ci.lvl = 0.83) %>% 
  rename(encRet = facet, reg = x, hitMiss = group, rt_res = predicted)


#overall plot

model_cond %>%
  ggplot(aes(x = hitMiss, y = rt_res, group = reg)) + 
  geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
  geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(1)) +
  geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(1)) + 
  facet_wrap(~encRet) + 
  geom_point(aes(fill = reg, color = reg, group = reg, x = hitMiss, y = rt_res), 
             size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.1), data = data) + 
  ylim(c(-400, 400))


model_cond %>%
  ggplot(aes(x = encRet, y = rt_res, group = reg)) + 
  geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
  geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(1)) +
  geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(1)) + 
  facet_wrap(~hitMiss) + 
  geom_point(aes(fill = reg, color = reg, group = reg, x = encRet, y = rt_res), 
             size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.1), data = data) + 
  ylim(c(-400, 400))



#plot of region effect alone
model_cond <- ggemmeans(data.lmer,
                        terms = c("reg"),
                        ci.lvl = 0.83) %>% 
  rename(reg = x, rt_res = predicted)

model_cond %>%
  ggplot(aes(x = reg, y = rt_res)) + 
  # geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
  geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(.2)) +
  geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(.2)) + 
  geom_point(aes(fill = reg, color = reg, group = reg, x = reg, y = rt_res), 
             size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.2), data = data) + 
  geom_boxplot(aes(x = reg, y = rt_res, fill = reg), outlier.shape = NA, alpha = .5, width = .1, colour = "black", data = data) +
  ylim(c(-500, 500))


  
  #plot of hitMiss effect alone
  model_cond <- ggemmeans(data.lmer,
                          terms = c("hitMiss"),
                          ci.lvl = 0.83) %>% 
    rename(hitMiss = x, rt_res = predicted)
  data$temp = paste(data$realID, data$encRet)
  model_cond %>%
    ggplot(aes(x = hitMiss, y = rt_res)) + 
    # geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
    geom_point(aes(fill = hitMiss, group = hitMiss, color = hitMiss), size = 1, shape = 19, alpha = 0.2, position = position_dodge(.2)) +
    geom_errorbar(aes(colour = hitMiss, ymin = conf.low, ymax = conf.high), 
                  width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(.2)) + 
    geom_point(aes(fill = hitMiss, color = hitMiss, group = hitMiss, x = hitMiss, y = rt_res), 
               size = 3, shape = 19, alpha = .2, position = position_dodge(.2), data = data) + 
    geom_boxplot(aes(x = hitMiss, y = rt_res, fill = hitMiss), 
                 outlier.shape = NA, alpha = .5, width = .1, colour = "black", data = data) +
    geom_line(aes(group = temp, x = hitMiss, y = rt_res), alpha = .05, linewidth = 2, data = data) +
  ylim(c(-500, 500))
  
  
  #plot of interaction between reg X encRet
  model_cond <- ggemmeans(data.lmer,
                          terms = c("reg", "encRet"),
                          ci.lvl = 0.83) %>% 
    rename(reg = x, rt_res = predicted, encRet = group)
  
  model_cond %>%
    ggplot(aes(x = encRet, y = rt_res, group = reg)) + 
    geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
    geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(1)) +
    geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                  width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(1)) + 
    geom_point(aes(fill = reg, color = reg, group = reg, x = encRet, y = rt_res), 
               size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.1), data = data) + 
    ylim(c(-400, 400))
  
  
  #####################################################################################################################
  ######################## analaysis of proportion of time ############################################################
  
  data$adjTime = data$centerOfMass / data$RT
  data.lmer <- lmer(adjTime ~ reg*hitMiss*encRet + (1|realID), REML = TRUE, 
                    control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE), data = data)
  
  Anova(data.lmer)
  summary(data.lmer)
  performance(data.lmer)

  

  model_cond <- ggemmeans(data.lmer,
                          terms = c("reg","hitMiss", "encRet"),
                          ci.lvl = 0.83) %>% 
    rename(encRet = facet, reg = x, hitMiss = group, adjTime = predicted)
  
  
  #overall plot
  
  model_cond %>%
    ggplot(aes(x = hitMiss, y = adjTime, group = reg)) + 
    geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
    geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(1)) +
    geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                  width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(1)) + 
    facet_wrap(~encRet) + 
    geom_point(aes(fill = reg, color = reg, group = reg, x = hitMiss, y = adjTime), 
               size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.1), data = data) + 
    ylim(c(.25, .75))
  
  
  
  
  
  #plot of region effect alone
  model_cond <- ggemmeans(data.lmer,
                          terms = c("reg"),
                          ci.lvl = 0.83) %>% 
    rename(reg = x, adjTime = predicted)
  
  model_cond %>%
    ggplot(aes(x = reg, y = adjTime)) + 
    # geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
    geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(.2)) +
    geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                  width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(.2)) + 
    geom_point(aes(fill = reg, color = reg, group = reg, x = reg, y = adjTime), 
               size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.2), data = data) + 
    geom_boxplot(aes(x = reg, y = adjTime, fill = reg), outlier.shape = NA, alpha = .5, width = .1, colour = "black", data = data)
    ylim(c(.25, .75))
  
  
  
    #plot of hitMiss effect alone
    model_cond <- ggemmeans(data.lmer,
                            terms = c("hitMiss"),
                            ci.lvl = 0.83) %>% 
      rename(hitMiss = x, adjTime = predicted)
    data$temp = paste(data$realID, data$encRet)
    model_cond %>%
      ggplot(aes(x = hitMiss, y = adjTime)) + 
      # geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
      geom_point(aes(fill = hitMiss, group = hitMiss, color = hitMiss), size = 1, shape = 19, alpha = 0.2, position = position_dodge(.2)) +
      geom_errorbar(aes(colour = hitMiss, ymin = conf.low, ymax = conf.high), 
                    width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(.2)) + 
      geom_point(aes(fill = hitMiss, color = hitMiss, group = hitMiss, x = hitMiss, y = adjTime), 
                 size = 1.5, shape = 19, alpha = .2, position = position_dodge(.2), data = data) + 
      geom_boxplot(aes(x = hitMiss, y = adjTime, fill = hitMiss), 
                   outlier.shape = NA, alpha = .5, width = .1, colour = "black", data = data) +
      geom_line(aes(group = temp, x = hitMiss, y = adjTime), alpha = .05, linewidth = 2, data = data) +
      ylim(c(.25, .75))
    
    
    # plot of encRet effect alone
  
    
    model_cond <- ggemmeans(data.lmer,
                            terms = c("encRet"),
                            ci.lvl = 0.83) %>% 
      rename(encRet = x, adjTime = predicted)
    data$temp = paste(data$realID, data$hitMiss)
    model_cond %>%
      ggplot(aes(x = encRet, y = adjTime)) + 
      # geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
      geom_point(aes(fill = encRet, group = encRet, color = encRet), size = 1, shape = 19, alpha = 0.2, position = position_dodge(.2)) +
      geom_errorbar(aes(colour = encRet, ymin = conf.low, ymax = conf.high), 
                    width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(.2)) + 
      geom_point(aes(fill = encRet, color = encRet, group = encRet, x = encRet, y = adjTime), 
                 size = 1.5, shape = 19, alpha = .2, position = position_dodge(.2), data = data) + 
      geom_boxplot(aes(x = encRet, y = adjTime, fill = encRet), 
                   outlier.shape = NA, alpha = .5, width = .1, colour = "black", data = data) +
      geom_line(aes(group = temp, x = encRet, y = adjTime), alpha = .05, linewidth = 2, data = data) +
      ylim(c(.25, .75))
    
  
  
  
    #plot of interaction between reg X encRet
    model_cond <- ggemmeans(data.lmer,
                            terms = c("reg", "encRet"),
                            ci.lvl = 0.83) %>% 
      rename(reg = x, adjTime = predicted, encRet = group)
    
    model_cond %>%
      ggplot(aes(x = encRet, y = adjTime, group = reg)) + 
      geom_line(aes(group = reg), alpha = 1 , colour = "gray48", linewidth = .5, position = position_dodge(1)) +
      geom_point(aes(fill = reg, group = reg, color = reg), size = 1, shape = 19, alpha = 0.2, position = position_dodge(1)) +
      geom_errorbar(aes(colour = reg, ymin = conf.low, ymax = conf.high), 
                    width = 0.4, alpha = 0.95, size = 1.3, position = position_dodge(1)) + 
      geom_point(aes(fill = reg, color = reg, group = reg, x = encRet, y = adjTime), 
                 size = 1.5, shape = 19, alpha = 1, position = position_jitterdodge(.1), data = data) + 
      ylim(c(.25, .75))
  
  
  
  
  
  ### model broken down to encoding only
  
    encDat = data %>% filter(encRet == 'enc' & (reg == 'hip' | reg == 'phg'))
    data.lmer <- lmer(rt_res ~ reg*hitMiss + (1|realID), REML = TRUE, 
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE), data = encDat)
    
    Anova(data.lmer)
    summary(data.lmer)
    performance(data.lmer)
  
  ### model broken down to retrieval only
  
    retDat = data %>% filter(encRet == 'ret' & (reg == 'hip' | reg == 'phg'))
    data.lmer <- lmer(rt_res ~ reg*hitMiss + (1|realID), REML = TRUE, 
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE), data = retDat)
    # data.lmer <- lm(adjTime ~ reg*hitMiss , data = retDat)
    
    Anova(data.lmer)
    summary(data.lmer)
    performance(data.lmer)
  
  
  ###########################################################################################

pd_2 <- position_dodge(width = .5)
ggplot(model_cond, aes(x=facet,y=predicted)) + 
  geom_point(aes(colour = group), position = pd_2) + 
  geom_errorbar(aes(colour = group, ymin = conf.low, ymax = conf.high), width = 0.4, alpha = 0.95, size = 1.3, position = pd_2) + 
  labs(colour = "region") +
  # scale_color_manual(values=c("turquoise4", "tan1","black")) + #spec num colors for num regions
  # geom_line(aes(group = facet, color = x), size = .8, position = pd_2) +
  # scale_y_reverse() +
  # facet_grid(x~facet)
  ylab("residual latency after regressing out RT") +
  geom_point(aes(x = hitMiss, y = rt_res, colour = reg), alpha = 0.3, size = 3, position = pd_2, data = data) + #raw data frame here
  xlab("Handedness") +
  theme_light() +
  theme(legend.position="right",
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        panel.border = element_rect(color="black", fill=NA, size = 1), 
        strip.background = element_rect(fill="gray", color="black"), 
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16))






















ggplot(data, aes(y = rt_res, x = hitMiss, fill = encRet)) +
  # geom_flat_violin(aes(fill = encRet), position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = F, alpha = .5, colour = NA) +
  geom_line(aes(group = realID), alpha = 0.1 , colour = "gray48", linewidth = .5, position = position_dodge(0.3)) +
  geom_point(aes(fill = encRet, group = realID, color = encRet), size = 1.5, shape = 19, alpha = 0.2, position = position_dodge(0.3)) +
  geom_boxplot(aes(x = hitMiss, y = rt_res, fill = encRet), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  #scale_y_continuous(limits = c(-3.5, -1.5), breaks = seq(-3.5, -1.5, by = 0.5)) +
  xlab("hit v. miss") + ylab ("latency residual") +
  scale_color_manual(values=c("turquoise4", "tan1"), name = "Task") +
  scale_fill_manual(values=c("turquoise4", "tan1"), name = "Task") +
  theme_bw() +
  facet_grid(encRet~reg) +
  theme(legend.position="NONE",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=alpha('blue', 0)),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold"))


#follow effects and plot the significant effects
ggplot(data, aes(y = rt_res, x = hitMiss)) +
  # geom_flat_violin(aes(fill = encRet), position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = F, alpha = .5, colour = NA) +
  geom_line(aes(group = realID), 
            alpha = 0.1 , colour = "gray48", 
            linewidth = .5, position = position_dodge(0.3)) +
  geom_point(aes(fill = hitMiss, group = realID, color = hitMiss), 
             size = 1.5, shape = 19, alpha = 0.2, position = position_dodge(0.3)) +
  geom_boxplot(aes(x = hitMiss, y = rt_res, fill = hitMiss), 
               outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  #scale_y_continuous(limits = c(-3.5, -1.5), breaks = seq(-3.5, -1.5, by = 0.5)) +
  xlab("hit v. miss") + ylab ("latency residual") +
  scale_color_manual(values=c("turquoise4", "tan1"), name = "Task") +
  scale_fill_manual(values=c("turquoise4", "tan1"), name = "Task") +
  theme_bw() +
  theme(legend.position="NONE",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=alpha('blue', 0)),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold"))




#subhit v. submiss

encDat = data %>% filter(cond == 'subMiss' | cond == 'subHit')



encDat.lmer <- lmer(centerOfMass ~ RT + reg * hitMiss + (1+chi|subID), REML = TRUE, 
              control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE), data = encDat)
encDat.lm = lm(centerOfMass ~ RT + reg * hitMiss,  data = encDat)

# test = aov(centerOfMass ~ cond*reg + Error(1+chi), data = encDat)
summary(encDat.lmer)
summary(encDat.lm)

data %>%
  ggplot(aes(x=RT, fill = hitMiss, y = centerOfMass)) +
  geom_line(aes(group = realID), alpha = .8 , colour = "gray48", linewidth = 1, position = position_dodge(0.3)) +
  geom_point(aes(fill = hitMiss, color = hitMiss), size = 1.5, shape = 19, alpha = 1, position = position_dodge(0.3)) + 
  facet_grid(encRet~reg)




encDat %>%
ggplot(aes(x = reg, y = centerOfMass, fill = cond)) +
   # geom_violin(draw_quantiles = c(.25, .5, .75), alpha = .5) + 
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.9)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4,
              position = position_jitterdodge(dodge.width = .9)) +

  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) +
  ggtitle('encoding data') +
  guides(color = guide_legend(override.aes = list(size=10)),
         shape = guide_legend(override.aes = list(size=10)))  

#hit v. miss retrieval data

retDat = data %>% filter(cond=='hit_on' | cond=='miss_on')

test = aov(centerOfMass ~ cond*reg + Error(chi), data = retDat)

library(lme4)
library(lmerTest)
library(car)

test2 <- lmer(centerOfMass ~ encRet * reg * hitMiss + (1+chi|subID), REML = TRUE, 
              control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE), data = data)
Anova(test2)
summary(test2)

retDat %>%
  ggplot(aes(x = reg, y = centerOfMass, fill = cond)) +
  # geom_violin(draw_quantiles = c(.25, .5, .75), alpha = .5) + 
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.9)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4,
              position = position_jitterdodge(dodge.width = .9)) +
  
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) +
  ggtitle('retrieval data') +
  guides(color = guide_legend(override.aes = list(size=10)),
         shape = guide_legend(override.aes = list(size=10)))  


#
data %>%
  ggplot(aes(x = encRet, y = centerOfMass, color = hitMiss)) + 
  geom_jitter()

ggplot(data, aes(y = centerOfMass, x = encRet, fill = hitMiss)) +
  # geom_flat_violin(aes(fill = encRet), position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = F, alpha = .5, colour = NA) +
  geom_line(aes(group = chi), alpha = 0.1 , colour = "gray48", linewidth = .5, position = position_dodge(0.3)) +
  geom_point(aes(fill = hitMiss, group = chi, color = hitMiss), size = 1.5, shape = 19, alpha = 0.2, position = position_dodge(0.3)) +
  geom_boxplot(aes(x = encRet, y = centerOfMass, fill = hitMiss), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  #scale_y_continuous(limits = c(-3.5, -1.5), breaks = seq(-3.5, -1.5, by = 0.5)) +
  xlab("Task [visual memory vs rest]") + ylab ("Aperiodic Slope") +
  scale_color_manual(values=c("turquoise4", "tan1"), name = "Task") +
  scale_fill_manual(values=c("turquoise4", "tan1"), name = "Task") +
  theme_bw() +
  facet_wrap(~reg) +
  theme(legend.position="NONE",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=alpha('blue', 0)),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold"))














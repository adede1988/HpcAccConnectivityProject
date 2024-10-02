ggplot(filter(long_df, encRet=='enc'), aes(y = Char, x = hitMiss, fill = hitMiss)) +
  geom_flat_violin(aes(fill =  hitMiss), position = position_nudge(x = .3, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  #geom_line(aes(group = realID), alpha = 0.2 , colour = "gray48", size = .8, position = position_dodge(0.3)) +
  geom_jitter(aes(fill = hitMiss, group = realID, color = hitMiss), size = 3, 
              shape = 19, alpha = 0.1,  
              position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = .1)) +
  geom_boxplot(aes(x = hitMiss, y = Char, fill = hitMiss), outlier.shape = NA, alpha = .5, width = .2,position = position_nudge(x = 0, y = 0), colour = "black") +
  ylab("Char path len") + xlab("hit miss") +
  scale_color_manual(values=c("turquoise4", "tan1"), name = "memory") +
  scale_fill_manual(values=c("turquoise4", "tan1"), name = "memory") +
  theme_bw() +
  theme(legend.position="right",
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


ggplot(filter(long_df, encRet=='enc'), aes(y = ST, x = hitMiss, fill = hitMiss)) +
  geom_flat_violin(aes(fill =  hitMiss), position = position_nudge(x = .3, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  #geom_line(aes(group = realID), alpha = 0.2 , colour = "gray48", size = .8, position = position_dodge(0.3)) +
  geom_jitter(aes(fill = hitMiss, group = realID, color = hitMiss), size = 3, 
              shape = 19, alpha = 0.1,  
              position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = .1)) +
  geom_boxplot(aes(x = hitMiss, y = ST, fill = hitMiss), outlier.shape = NA, alpha = .5, width = .2,position = position_nudge(x = 0, y = 0), colour = "black") +
  ylab("node strength (weighted)") + xlab("hit miss") +
  scale_color_manual(values=c("turquoise4", "tan1"), name = "memory") +
  scale_fill_manual(values=c("turquoise4", "tan1"), name = "memory") +
  theme_bw() +
  theme(legend.position="right",
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


ggplot(filter(long_df, encRet=='enc'), aes(y = Saturation, x = hitMiss, fill = hitMiss)) +
  geom_flat_violin(aes(fill =  hitMiss), position = position_nudge(x = .3, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  #geom_line(aes(group = realID), alpha = 0.2 , colour = "gray48", size = .8, position = position_dodge(0.3)) +
  geom_jitter(aes(fill = hitMiss, group = realID, color = hitMiss), size = 3, 
              shape = 19, alpha = 0.1,  
              position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = .1)) +
  geom_boxplot(aes(x = hitMiss, y = Saturation, fill = hitMiss), outlier.shape = NA, alpha = .5, width = .2,position = position_nudge(x = 0, y = 0), colour = "black") +
  ylab("node strength (binarized)") + xlab("hit miss") +
  scale_color_manual(values=c("turquoise4", "tan1"), name = "memory") +
  scale_fill_manual(values=c("turquoise4", "tan1"), name = "memory") +
  theme_bw() +
  theme(legend.position="right",
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

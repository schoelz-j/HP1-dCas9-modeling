library(ggridges)

# Ridge plots of simulations results
expression <- read.csv(here('datasets/S2_binding_expression.csv'))
expression <- expression %>%
  group_by(FlyBase_FBgn) %>%
  summarise(Avg_TPM = mean(S2)) %>%
  mutate(log_TPM = log(Avg_TPM)) %>%
  select(FlyBase_FBgn, log_TPM) %>%
  filter(is.finite(log_TPM))
colnames(expression)[1] <- 'FBgn'
expression$log_TPM[expression$FBgn==dcas_fbgns[1]]
baselines <- c(expression$log_TPM[expression$FBgn==dcas_fbgns[1]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[2]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[3]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[4]])
# Cpr11B
cpr11b_plus <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/cpr11b_hp1_plus_simulated.csv')
str(cpr11b_plus)

cpr11b_ridge <- data.frame(values = c(cpr11b_plus$HP1a, cpr11b_plus$HP1B, cpr11b_plus$HP1C),
                           model = c(rep('HP1a', 1000), rep('HP1B', 1000), rep('HP1C', 1000)))
cpr11b_plus_pred <- ggplot(cpr11b_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[1], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

cpr11b_only <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/cpr11b_hp1_only_simulated.csv')
cpr11b_only_ridge <- data.frame(values = c(cpr11b_only$HP1a, cpr11b_only$HP1B, cpr11b_only$HP1C),
                                model = cpr11b_ridge$model)
cpr11b_only_pred <- ggplot(cpr11b_only_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[1], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

# Dgt3
dgt3_plus <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/dgt3_hp1_plus_simulated.csv')
dgt3_ridge <- data.frame(values = c(dgt3_plus$HP1a, dgt3_plus$HP1B, dgt3_plus$HP1C),
                         model = cpr11b_ridge$model)
dgt3_plus_pred <- ggplot(dgt3_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[2], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

dgt3_only <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/dgt3_hp1_only_simulated.csv')
dgt3_only_ridge <- data.frame(values = c(dgt3_only$HP1a, dgt3_only$HP1B, dgt3_only$HP1C),
                              model = dgt3_ridge$model)
dgt3_only_pred <- ggplot(dgt3_only_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[2], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

# CecA1
ceca1_plus <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/ceca1_hp1_plus_simulated.csv')
ceca1_ridge <- data.frame(values = c(ceca1_plus$HP1a, ceca1_plus$HP1B, ceca1_plus$HP1C),
                          model = dgt3_ridge$model)

ceca1_plus_pred <- ggplot(ceca1_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[3], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

ceca1_only <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/ceca1_hp1_only_simulated.csv')
ceca1_only_ridge <- data.frame(values = c(ceca1_only$HP1a, ceca1_only$HP1B, ceca1_only$HP1C),
                               model = ceca1_ridge$model)
ceca1_only_pred <- ggplot(ceca1_only_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[3], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

# Mtk
mtk_plus <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/mtk_hp1_plus_simulated.csv')
mtk_ridge <- data.frame(values = c(mtk_plus$HP1a, mtk_plus$HP1B, mtk_plus$HP1C),
                        model = dgt3_ridge$model)
mtk_plus_pred <- ggplot(mtk_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[3], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

mtk_only <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_outputs/mtk_hp1_only_simulated.csv')
mtk_only_ridge <- data.frame(values = c(mtk_only$HP1a, mtk_only$HP1B, mtk_only$HP1C),
                             model = mtk_ridge$model)
mtk_only_pred <- ggplot(mtk_only_ridge, aes(x = values, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[3], lwd = 1.5, color = 'grey')+
  scale_fill_manual(values = c('#473E92', '#2D8775', '#8F3086'))+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')

pdf('Figure4_2.pdf', height = 10, width = 9)
ggarrange(cpr11b_results, cpr11b_only_pred, cpr11b_plus_pred,
          dgt3_results, dgt3_only_pred, dgt3_plus_pred,
          ceca1_results, ceca1_only_pred, ceca1_plus_pred,
          mtk_results, mtk_only_pred, mtk_plus_pred,
          nrow = 4, ncol = 3, labels = c('A.', 'B.', 'C.', 'D.', 'E.', 'F.', 'G.', 'H.', 'I.',
                                         'J.', 'K.', 'L.'))
dev.off()
getwd()


# Ridge plots of predictions for remaining genes in other figures
dcas_fbgns <- c('FBgn0033483', 'FBgn0033649', 'FBgn0038965')
baselines <- c(expression$log_TPM[expression$FBgn==dcas_fbgns[1]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[2]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[3]])

ggplot(expression, aes(x = log_TPM))+
  geom_histogram()+
  geom_vline(xintercept = baselines[3])

hp1a_plus_preds <- read.csv('hp1a_hp1_plus_simulated.csv')
hp1a_only_preds <- read.csv('hp1a_hp1_only_simulated.csv')
hp1a_plus_preds$model <- rep('HP1 Plus', nrow(hp1a_plus_preds))
hp1a_only_preds$model <- rep('HP1 Only', nrow(hp1a_only_preds))
hp1a <- data.frame(rbind(hp1a_plus_preds, hp1a_only_preds))


hp1a_egr <- ggplot(hp1a, aes(x = egr, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[1], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#3d5a80', '#98c1d9'))

hp1a_pyr <- ggplot(hp1a, aes(x = pyr, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[2], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#3d5a80', '#98c1d9'))

hp1a_mats <- ggplot(hp1a, aes(x = mats, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[3], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#3d5a80', '#98c1d9'))

dcas_fbgns <- c('FBgn0030613', 'FBgn0003386', 'FBgn0035844')
baselines <- c(expression$log_TPM[expression$FBgn==dcas_fbgns[1]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[2]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[3]])

ggplot(expression, aes(x = log_TPM))+
  geom_histogram()+
  geom_vline(xintercept = baselines[3])
hp1b_plus_preds <- read.csv('hp1b_hp1_plus_simulated.csv')
hp1b_only_preds <- read.csv('hp1b_hp1_only_simulated.csv')
hp1b_plus_preds$model <- rep('HP1 Plus', nrow(hp1b_plus_preds))
hp1b_only_preds$model <- rep('HP1 Only', nrow(hp1b_only_preds))
hp1b <- data.frame(rbind(hp1b_plus_preds, hp1b_only_preds))

hp1b_rab3 <- ggplot(hp1b, aes(x = rab3, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[1], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#86B049', '#DFF5CE'))

hp1b_shaw <- ggplot(hp1b, aes(x = shaw, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[2], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#86B049', '#DFF5CE'))

hp1b_cg76 <- ggplot(hp1b, aes(x = cg76, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[3], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#86B049', '#DFF5CE'))

dcas_fbgns <- c('FBgn0003884', 'FBgn0033794')
baselines <- c(expression$log_TPM[expression$FBgn==dcas_fbgns[1]],
               expression$log_TPM[expression$FBgn==dcas_fbgns[2]])

hp1c_plus_preds <- read.csv('hp1c_hp1_plus_simulated.csv')
hp1c_only_preds <- read.csv('hp1c_hp1_only_simulated.csv')
hp1c_plus_preds$model <- rep('HP1 Plus', nrow(hp1c_plus_preds))
hp1c_only_preds$model <- rep('HP1 Only', nrow(hp1c_only_preds))
hp1c <- data.frame(rbind(hp1c_plus_preds, hp1c_only_preds))

hp1c_alpha <- ggplot(hp1c, aes(x = alpha, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[1], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#89609e', '#dfc0eb'))

hp1c_cg26 <- ggplot(hp1c, aes(x = cg26, y = model, fill = model))+
  geom_density_ridges(scale = 20, alpha = 0.8, color = 'white')+
  geom_vline(xintercept = baselines[1], lwd = 1.5, color = 'grey')+
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  ylab('')+
  xlab('Predicted Expression')+
  scale_fill_manual(values = c('#89609e', '#dfc0eb'))

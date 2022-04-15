setwd("Mohammad/")

library(tidyverse)
library(ggplot2)
library(ggtree)
library(reshape2)

library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)





tree <- read.tree(paste0("all_nuc.treefile"))
info <- read.csv2("order_family_species.csv")



df <- info[, 3:5]
rownames(df) <- info$full_tip_label



g <- left_join(tibble(full_tip_label = tree$tip.label), info[, 2:5])



p <-
  ggtree(tree) +
    geom_tiplab(size = 2, align = T, linetype = "blank", hjust = .5) +
    # geom_text(aes(label = node)) +
    theme_tree2()

p


p1_a <- 
  p %<+% 
    g + 
    geom_tippoint(aes(color = GC_skew)) +
    scale_color_continuous(low = 'red', high = 'green') +
    theme(legend.position = "right")

p1_a



p1_b <- 
  p %<+% 
    g + 
    geom_tippoint(aes(color = Total_GC)) +
    scale_color_continuous(low = 'red', high = 'green') +
    theme(legend.position = "right")

p1_b



p1_c <- 
  p %<+% 
    g + 
    geom_tippoint(aes(color = Seq_length)) +
    scale_color_continuous(low = 'red', high = 'green') +
    theme(legend.position = "right")

p1_c



# svg(filename = "images/GC_skew_1.svg", width = 10, height = 10)
# p1_a
# dev.off()
# 
# svg(filename = "images/Total_GC_1.svg", width = 10, height = 10)
# p1_b
# dev.off()
# 
# svg(filename = "images/Seq_length_1.svg", width = 10, height = 10)
# p1_c
# dev.off()





####################################
#color_scale <- c("Seq_length" = "red", "Total_GC" = "blue", "GC_skew" = "purple")

# p2 <- 
#   df %>% 
#     mutate(ord = 1:n()) %>% 
#     melt(id.vars = c("ord")) %>% 
#     group_by(variable) %>%
#     mutate(col_val = scale(value)) %>%
#     ungroup() %>%
#     arrange(variable, ord) %>% 
# 
#     ggplot(aes(x = variable, y = ord, fill = variable, alpha = col_val)) + 
#       geom_tile() +
#       scale_fill_manual(values = color_scale) +
#       scale_alpha(range = c(0.1, 0.9)) +
#       theme_tree2() +
#       theme(legend.position = "none") +
# 
# p2
####################################



library(aplot)
library(ggnewscale)


p2 <- 
  df %>% 
    mutate(ord = 1:n()) %>%
    
    ggplot() +
      
      geom_tile(aes(x = rep("GC_skew", nrow(df)), y = ord, fill = GC_skew)) +
      # scale_fill_gradient2("GC_skew", limits = range(df$GC_skew), 
      #                      low = "#9ecae1", high = "#3182bd") +
      
      new_scale("fill") +
      
      geom_tile(aes(x = rep("Seq_length", nrow(df)), y = ord, fill = Seq_length)) +
      scale_fill_gradient2("Seq_length", limits = range(df$Seq_length),
                           low = "#762A83", high = "#1B7837") +
      
      new_scale("fill") +
      
      geom_tile(aes(x = rep("Total_GC", nrow(df)), y = ord, fill = Total_GC)) +
      scale_fill_gradient2("Total_GC", limits = range(df$Total_GC), 
                           low = "#1B7837", high = "#762A83") +
      
      theme_tree2()

p2



p3 <- p %>% insert_right(p2, width = 0.5)
p3



# svg(filename = "images/all_variables.svg", width = 15, height = 10)
# p3
# dev.off()



####################################



df <- info[, 3:5]
rownames(df) <- info$full_tip_label



p4d <- phylo4d(tree, df)
barplot(p4d)

# svg(filename = "images/phylo4d_barplot.svg", width = 10, height = 16)
# barplot(p4d)
# dev.off()





ps <- phyloSignal(p4d = p4d, method = "all")
ps

# write_csv(ps$stat, "phyloSignal_statistics.csv")
# write_csv(ps$pvalue, "phyloSignal_pvalues.csv")



# this will take a while but it's not too bad
phylosim <- phyloSim(tree = tree, method = "all", nsim = 100, reps = 99) 



plot(phylosim, stacked.methods = F, quantiles = c(0.05, 0.95))

# svg(filename = "images/phylosim.svg", width = 20, height = 10)
# plot(phylosim, stacked.methods = F, quantiles = c(0.05, 0.95))
# dev.off()



plot.phylosim(phylosim, what = "pval", stacked.methods = T)

# svg(filename = "images/phylosim_stacked.svg", width = 20, height = 10)
# plot.phylosim(phylosim, what = "pval", stacked.methods = T)
# dev.off()



GC_skew <- phyloCorrelogram(p4d, trait = "GC_skew")
Seq_length <- phyloCorrelogram(p4d, trait = "Seq_length")
Total_GC <- phyloCorrelogram(p4d, trait = "Total_GC")

plot(GC_skew)
plot(Seq_length)
plot(Total_GC)

svg(filename = "images/GC_skew_Correlogram.svg", width = 10, height = 10)
plot(GC_skew)
dev.off()
# 
svg(filename = "images/Seq_length_Correlogram.svg", width = 10, height = 10)
plot(Seq_length)
dev.off()
# 
svg(filename = "images/Total_GC_Correlogram.svg", width = 10, height = 10)
plot(Total_GC)
dev.off()





collembola <- lipaMoran(p4d)
collembola.p4d <- lipaMoran(p4d, as.p4d = T)



barplot.phylo4d(p4d, bar.col = (collembola$p.value < 0.05) + 1,
                center = F , scale = F)

# svg(filename = "images/collembola_signal_1.svg", width = 10, height = 16)
# barplot.phylo4d(p4d, bar.col = (collembola$p.value < 0.05) + 1,
#                 center = F , scale = F)
# dev.off()



barplot.phylo4d(collembola.p4d, bar.col = (collembola$p.value < 0.05) + 1, 
                center = F, scale = F)

# svg(filename = "images/collembola_signal_2.svg", width = 10, height = 16)
# barplot.phylo4d(collembola.p4d, bar.col = (collembola$p.value < 0.05) + 1, 
#                 center = F, scale = F)
# dev.off()
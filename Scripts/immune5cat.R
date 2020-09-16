# 5 category immune
# Re-stratify
slide_out = read.csv('/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/Test_slide.csv')
tile_out = read.csv('/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/Test_tile.csv')
ref = read.csv('/Users/rh2740/documents/GBM/Results/5cat/GBM_immune_subtype_score.txt', sep='\t')
ref$immune = gsub('im', '', ref$immune)

colnames(ref) = gsub('case', 'slide', colnames(ref))
new_slide = merge(slide_out, ref[,1:2])
new_tile = merge(tile_out, ref[,1:2])

new_slide$stratify1 = (new_slide$immune < 3)
new_slide$stratify2 = (new_slide$immune != 3)
new_slide$stratify1 = gsub('TRUE', 'high', new_slide$stratify1)
new_slide$stratify2 = gsub('FALSE', 'low', new_slide$stratify2)
new_slide$stratify2 = gsub('TRUE', 'high', new_slide$stratify2)
new_slide$stratify1 = gsub('FALSE', 'low', new_slide$stratify1)

new_tile$stratify1 = (new_tile$immune < 3)
new_tile$stratify2 = (new_tile$immune != 3)
new_tile$stratify1 = gsub('TRUE', 'high', new_tile$stratify1)
new_tile$stratify2 = gsub('FALSE', 'low', new_tile$stratify2)
new_tile$stratify2 = gsub('TRUE', 'high', new_tile$stratify2)
new_tile$stratify1 = gsub('FALSE', 'low', new_tile$stratify1)

write.csv(new_slide, '/Users/rh2740/documents/GBM/Results/5cat/slide.csv')
write.csv(new_tile, '/Users/rh2740/documents/GBM/Results/5cat/tile.csv')

# Calculate new AUROC
slide = read.csv('/Users/rh2740/documents/GBM/Results/5cat/slide.csv')
tile = read.csv('/Users/rh2740/documents/GBM/Results/5cat/tile.csv')
library(pROC)
library(MLmetrics)
slide_roc.1 = roc(factor(slide$stratify1), slide$high_score, levels=c('high', 'low'))
slide_roc.2 = roc(factor(slide$stratify2), slide$high_score, levels=c('high', 'low'))
tile_roc.1 = roc(factor(tile$stratify1), tile$high_score, levels=c('high', 'low'))
tile_roc.2 = roc(factor(tile$stratify2), tile$high_score, levels=c('high', 'low'))
slide_roc.old = roc(factor(slide$True_label), slide$high_score, levels=c('high', 'low'))
tile_roc.old = roc(factor(tile$True_label), tile$high_score, levels=c('high', 'low'))

slide_PRAUC.1 = PRAUC(slide$high_score, factor(slide$stratify1))
slide_PRAUC.2 = PRAUC(slide$high_score, factor(slide$stratify2))
tile_PRAUC.1 = PRAUC(tile$high_score, factor(tile$stratify1))
tile_PRAUC.2 = PRAUC(tile$high_score, factor(tile$stratify2))
slide_PRAUC.old = PRAUC(slide$high_score, factor(slide$True_label))
tile_PRAUC.old = PRAUC(tile$high_score, factor(tile$True_label))

# New plots
input_file=read.csv('/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/tSNE_Figure.csv')
ref = read.csv('/Users/rh2740/documents/GBM/Results/5cat/GBM_immune_subtype_score.txt', sep='\t')
colnames(ref) = gsub('case', 'slide', colnames(ref))

tsne = merge(input_file, ref[,1:2])

slide = read.csv('/Users/rh2740/documents/GBM/Results/5cat/slide.csv')
tsne = merge(tsne, slide$stratify2)
tsne$immune_label = factor(tsne$y)
write.csv(tsne, '/Users/rh2740/documents/GBM/Results/5cat/tsne.csv')
tsne['immune_label'] <- lapply(tsne['immune_label'] , factor)

library(ggplot2)
library(gridExtra)

tsne = read.csv('/Users/rh2740/documents/GBM/Results/5cat/tsne.csv')
pa=ggplot(data=tsne,aes_string(x='tsne1',y='tsne2',col='immune'))+
  scale_color_manual(values = c("#8E063B", "#D95F02", "#7570B3", "#F6C971", "#949E15"))+
  geom_point(alpha=1, size=1)+ scale_shape(solid = TRUE)+
  #theme(legend.position='bottom')+
  xlim(-60,60)+
  ylim(-60,60)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


pdf(file='/Users/rh2740/documents/GBM/Results/5cat/5catimmune.pdf',
    width=10,height=10,useDingbats=FALSE)

grid.arrange(pa, nrow=1)

dev.off()


#4cat immune
slide_out = read.csv('/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/Test_slide.csv')
tile_out = read.csv('/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/Test_tile.csv')
ref = read.csv('/Users/rh2740/documents/GBM/Results/4cat/gbm_all_subtype_collections.v5.1.tsv', sep='\t')
ref$immune = gsub('im', '', ref$immune)

colnames(ref) = gsub('case', 'slide', colnames(ref))
new_slide = merge(slide_out, ref[,c(1,15)])
new_tile = merge(tile_out, ref[,c(1,15)])

new_slide$stratify1 = (new_slide$immune < 4)
new_slide$stratify1 = gsub('TRUE', 'high', new_slide$stratify1)
new_slide$stratify1 = gsub('FALSE', 'low', new_slide$stratify1)

new_tile$stratify1 = (new_tile$immune < 3)
new_tile$stratify1 = gsub('TRUE', 'high', new_tile$stratify1)
new_tile$stratify1 = gsub('FALSE', 'low', new_tile$stratify1)

write.csv(new_slide, '/Users/rh2740/documents/GBM/Results/4cat/slide.csv')
write.csv(new_tile, '/Users/rh2740/documents/GBM/Results/4cat/tile.csv')

# Calculate new AUROC
slide = read.csv('/Users/rh2740/documents/GBM/Results/4cat/slide.csv')
tile = read.csv('/Users/rh2740/documents/GBM/Results/4cat/tile.csv')

slide.trunc = slide[slide$immune != 3,]
tile.trunc = tile[tile$immune != 3,]
library(pROC)
library(MLmetrics)
slide_roc.1 = roc(factor(slide$stratify1), slide$high_score, levels=c('high', 'low'))
slide.trunc_roc = roc(factor(slide.trunc$stratify1), slide.trunc$high_score, levels=c('high', 'low'))
tile_roc.1 = roc(factor(tile$stratify1), tile$high_score, levels=c('high', 'low'))
tile.trunc_roc = roc(factor(tile.trunc$stratify1), tile.trunc$high_score, levels=c('high', 'low'))
slide_roc.old = roc(factor(slide$True_label), slide$high_score, levels=c('high', 'low'))
tile_roc.old = roc(factor(tile$True_label), tile$high_score, levels=c('high', 'low'))

slide_PRAUC.1 = PRAUC(slide$high_score, factor(slide$stratify1))
slide_PRAUC.1 = PRAUC(slide$high_score, factor(slide$stratify1))
tile_PRAUC.1 = PRAUC(tile$high_score, factor(tile$stratify1))
slide_PRAUC.1 = PRAUC(slide$high_score, factor(slide$stratify1))
slide_PRAUC.old = PRAUC(slide$high_score, factor(slide$True_label))
tile_PRAUC.old = PRAUC(tile$high_score, factor(tile$True_label))

# New plots
input_file=read.csv('/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/tSNE_Figure.csv')
ref = read.csv('/Users/rh2740/documents/GBM/Results/4cat/GBM_immune_subtype_score.txt', sep='\t')
colnames(ref) = gsub('case', 'slide', colnames(ref))

tsne = merge(input_file, ref[,1:2])

slide = read.csv('/Users/rh2740/documents/GBM/Results/4cat/slide.csv')
tsne = merge(tsne, slide$stratify1)
tsne$immune_label = factor(tsne$y)
write.csv(tsne, '/Users/rh2740/documents/GBM/Results/4cat/tsne.csv')
tsne['immune_label'] <- lapply(tsne['immune_label'] , factor)

library(ggplot2)
library(gridExtra)

tsne = read.csv('/Users/rh2740/documents/GBM/Results/4cat/tsne.csv')
pa=ggplot(data=tsne,aes_string(x='tsne1',y='tsne2',col='immune'))+
  scale_color_manual(values = c("#8E063B", "#D95F02", "#F6C971", "#7570B3"))+
  geom_point(alpha=1, size=1)+ scale_shape(solid = TRUE)+
  #theme(legend.position='bottom')+
  xlim(-60,60)+
  ylim(-60,60)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


pdf(file='/Users/rh2740/documents/GBM/Results/4cat/4catimmune.pdf',
    width=10,height=10,useDingbats=FALSE)

grid.arrange(pa, nrow=1)

dev.off()




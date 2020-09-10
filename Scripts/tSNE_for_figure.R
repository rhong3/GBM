## reduce dimensionality of the acitvation layer of model
## visualize the manifold


input_file='/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/For_tSNE.csv'
output_file='/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/tSNE_Figure.csv'
out_fig='/Users/rh2740/documents/GBM/Results/1109/FS1immune0/out/Figure.pdf'
start=14   # telomere 15, others 14
bins=50
POS_score=c('high_score')
MDP = 0.5 # 0.5 for binary; 1/length(POS_score)

library(Rtsne)
ori_dat = read.table(file=input_file,header=T,sep=',')
# P = ori_dat[which(ori_dat$Prediction==1),]
# N = ori_dat[which(ori_dat$Prediction==0),]
# N = ori_dat[sample(nrow(N), 10000), ]
# sp_ori_dat = rbind(P, N)
sp_ori_dat=ori_dat[sample(nrow(ori_dat), 20000), ]
# sp_ori_dat=ori_dat
X = as.matrix(sp_ori_dat[,start:dim(sp_ori_dat)[2]])
res = Rtsne(X, initial_dims=100, check_duplicates = FALSE)
Y=res$Y
out_dat = cbind(sp_ori_dat[,1:(start-1)],Y)

dat = cbind(out_dat,x_bin=cut(out_dat[,start],bins),
            y_bin=cut(out_dat[,(start+1)],bins))

dat = cbind(dat, x_int = as.numeric(dat$x_bin),
            y_int = as.numeric(dat$y_bin))

colnames(dat)[start:(start+1)]=c('tsne1','tsne2')

dat$True_label=as.factor(dat$True_label)
dat$slide=as.factor(dat$slide)

write.table(dat, file=output_file, row.names = F, sep=',')

## plot the manifold with probability
library(ggplot2)
library(gridExtra)
palist <- list()
pblist <- list()
for(i in 1:length(POS_score)){
  palist[[i]]=ggplot(data=dat,aes_string(x='tsne1',y='tsne2',col=POS_score[i]))+
    scale_color_gradient2(high='red',mid='gray',low='steelblue',midpoint=MDP)+
    geom_point(alpha=1, size=1)+ scale_shape(solid = TRUE)+
    #theme(legend.position='bottom')+
    xlim(-60,60)+
    ylim(-60,60)+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}

pdf(file=out_fig,
    width=10,height=10,useDingbats=FALSE)

grid.arrange(palist[[1]], nrow=1)

dev.off()

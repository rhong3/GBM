## reduce dimensionality of the acitvation layer of model
## visualize the manifold


input_file='/Users/rh2740/documents/GBM/Results/1109/F1telomere0/out/For_tSNE.csv'
output_file='/Users/rh2740/documents/GBM/Results/1109/F1telomere0/out/tSNE_P_N.csv'
out_fig='/Users/rh2740/documents/GBM/Results/1109/F1telomere0/out/P_N.pdf'
start=15   # telomere 15, others 14
bins=50
POS_score=c('short_score','normal_score','long_score')
MDP = 1/length(POS_score) # 0.5 for binary; 1/length(POS_score)

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
    scale_color_gradient2(high='darkorange',mid='white',low='steelblue',midpoint=MDP)+
    geom_point(alpha=0.2)+
    #theme(legend.position='bottom')+
    xlim(-60,60)+
    ylim(-60,60)
  
  pblist[[i]]=ggplot(data=dat,aes_string(x='tsne1',y='tsne2'))+
    geom_point(aes(col=True_label),alpha=0.2)+
    scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    #theme(legend.position='bottom')+
    xlim(-60,60)+
    ylim(-60,60)
}

p3=ggplot(data=dat,aes_string(x='tsne1',y='tsne2'))+
  geom_point(aes(col=slide),alpha=0.2)+
  theme(legend.position='none')+
  xlim(-60,60)+
  ylim(-60,60)

p4=ggplot(data=subset(dat,True_label==1),aes_string(x='tsne1',y='tsne2'))+
  geom_point(aes(col=slide),alpha=0.2)+
  theme(legend.position='none')+
  xlim(-60,60)+
  ylim(-60,60)

pdf(file=out_fig,
    width=14,height=7)

for(i in 1:length(palist)){
  grid.arrange(palist[[i]],pblist[[i]],nrow=1)
}
grid.arrange(p3,p4,nrow=1)

dev.off()

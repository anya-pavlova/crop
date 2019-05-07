options(stringsAsFactors = F)
d = read.csv('/home/anna/Documents/Agro_ms/sunflover_artecls/fatty_acids_results_no_blank_normal_plus.csv')
dim(d)
d[1:10,1:7]
hist(d[,-1:-3])
d = list(peaks=d[,1:4],intensity=as.matrix(d[,-1:-4]))


mean.int=apply(d$intensity,1,mean,na.rm=T)
mean.int=mean.int-min(mean.int)
plot(d$peaks$rt,d$peaks$mz,pch=19,cex=mean.int/4+0.1)

colnames(d$intensity)[order(apply(is.na(d$intensity),2,mean))]

cor = cor(d$intensity,u='p')
dim(cor)
mds = cmdscale(1-cor,k=2)
plot(mds,col=grepl('QC',colnames(cor))+1)

rn = gsub('S|new','',sapply(strsplit(colnames(d$intensity),'_'),'[',4))
qc = rn %in% c('ACN','QC')
rn = as.numeric(rn)
table(is.na(rn),qc)
batch = rn
batch[] = 1
batch[rn > 97] = 2
batch[rn > 193] = 3
batch[rn > 276] = 4
batch[rn > 346] = 5
batch[rn > 441] = 6
batch[rn > 536] = 7
batch[qc] = 8

qc.batch = as.numeric(substr(sapply(strsplit(colnames(d$intensity),'_'),'[',5),1,1))
batch[qc] = qc.batch[qc] 
table(batch)

table(batch)
range(rn,na.rm=T)
hist(rn[batch==1])

plot(mds,col=batch)
legend('topleft',col=1:8,legend=1:8,pch=19)

f = batch %in% c(1:2,5:7)
f=T
mds = cmdscale(1-cor[f,f],k=2)
plot(mds,col=batch[f],xlim=c(-0.1,0.1),ylim=c(-0.1,0.1))
plot(mds,col=grepl('QC',colnames(cor))+1)

####### look on GC
gc = read.csv('/home/anna/Documents/Agro_ms/sunflover_artecls/need_sunf_samples.csv',header = T,skip = 1,check.names = F,dec = ',',fileEncoding = 'utf8')

rownames(gc) = paste0('S',gc$extraction_number)
gc = t(as.matrix(gc[,-1:-2]))

lc = d$intensity
rownames(lc)=d$peaks$fatty_acid
colnames(lc) = paste0('S',rn)
rownames(gc) = gsub('С','C',rownames(gc))

cmn.sid = intersect(colnames(lc),colnames(gc))
cmn.fa =  intersect(rownames(lc),rownames(gc))

colnames(gc)[!(colnames(gc) %in% cmn.sid)]

gc = gc[cmn.fa,cmn.sid]
lc = lc[cmn.fa,cmn.sid]

lc = 2^lc
lc = sweep(lc,2,apply(lc,2,sum),'/')*100
#gc = log2(gc)
par(mfrow=c(3,4))
for(i in 1:nrow(gc)){
  #f = !is.na(gc[i,]) & !is.na(lc[i,])
  plot(gc[i,],lc[i,],xlab='GC',ylab='LC',main=rownames(gc)[i])
  abline(lm(lc[i,] ~ gc[i,]),col='red')
  legend('topleft',legend=paste0('rho=',round(cor(gc[i,],lc[i,],m='sp'),3)))
}

plot(apply(gc,1,mean),apply(lc,1,mean),log='xy')
text(apply(gc,1,mean),apply(lc,1,mean),rownames(gc),adj=c(0,1))
plot(d$peaks$rt,d$peaks$mz)

col = sample(rainbow(12))
plot(gc,lc,log='xy',col=rep(col,times=49),pch=1:12)
legend('bottomright',col=col,pch=1:12,legend=rownames(gc))
abline(a=0,b=1,col='red')

plot(d$peaks$rt,d$peaks$mz,cex=mean.int/4+0.1,pch=19)
text(d$peaks$rt,d$peaks$mz,d$peaks$fatty_acid)
text(d$peaks$rt,d$peaks$mz,d$peaks$fatty_acid,adj=c(0,1))

boxplot(mean.int ~ I(d$peaks$fatty_acid %in% rownames(gc)))

names(mean.int) = d$peaks$fatty_acid
gc.mean=log2(apply(gc,1,mean)[names(mean.int)])
mi = cbind(mean.int-mean(mean.int),gc.mean-mean(gc.mean,na.rm = T))
mi = mi[order(mi[,1]),]
image(1:2,1:nrow(mi),t(mi),col=topo.colors(100),xaxt = 'n',yaxt = 'n',xlab='Method',ylab='FA')
axis(1,1:2,labels = c('LC','GC'))
axis(2,1:nrow(mi),labels = rownames(mi),las=2)

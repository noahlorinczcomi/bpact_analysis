# shared count SIMEX function ####
library(data.table);library(dplyr);library(magrittr)
rm(list=ls(all=TRUE))
shared_count_simex=function(
    M=8148, # number of jointly tested genes
    # delta1=0.98, # proportion of non-causal trait 1 genes
    # delta2=0.99, # proportion of non-causal trait 2 genes
    Mcausal1=10, # estimated number of causal genes for trait 1
    Mcausal2=10, # estimated number of causal genes for trait 1
    ngwas1=450000, # sample size of trait 1 GWAS
    ngwas2=16000, # sample size of trait 2 GWAS
    estimated_overlap_count=10, # data-estimated number of overlapping causal genes
    upsilon_overlap=0.2, # n01/sqrt(n1*n2)*Corr(t1,t2)
    # assumed simulation parameters
    h21=0.2, # SNP heritability trait 1
    h22=0.2, # SNP heritability trait 1
    niter=1000, # number of iterations
    m=50, # number of tested SNPs per gene
    mcausal_snps_per_gene=3, # number of causal SNPs per causal gene
    LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
    LD_rho=0.5, # correlation parameter
    nlambdas=10, # number of lambdas to evaluate in SIMEX
    doplot=TRUE, # should a plot be created?
    verbose=FALSE,
    do.extrapolatation=TRUE
) {
  # this function is used to estimate the number of genes which are inferred to be shared
  # due only to sample overlap
  library(mvnfast,lib='/home/lorincn/Rpkgs')
  tr=function(x) sum(diag(x))
  cs=function(n,rho) {a=matrix(rho,n,n);diag(n)=1}
  # genome-wide parameters
  ixs=1:M # gene indices
  delta1=1-Mcausal1/18257
  delta2=1-Mcausal2/18257
  # mcausalgenes1=round((1-delta1)*M) # number causal genes for trait 1
  # mcausalgenes2=round((1-delta2)*M) # number causal genes for trait 2
  mcausalgenes1=Mcausal1
  mcausalgenes2=Mcausal2
  ix1=sample(ixs,mcausalgenes1,replace=FALSE) # indices of causal trait 1 genes
  if(estimated_overlap_count==0) {
    ix2=sample(ixs[-ix1],mcausalgenes2,replace=FALSE) # indices of causal trait 2 genes
  } else {
    ix2_a=sample(ix1,estimated_overlap_count,replace=FALSE)
    ix2_b=sample(ixs[-c(ix1,ix2_a)],mcausalgenes2-estimated_overlap_count,replace=FALSE)
    ix2=c(ix2_a,ix2_b)
  }
  # SNP-level parameters
  causalsnpix=round(seq(1,m,length.out=mcausal_snps_per_gene)) # indices of causal SNPs for each gene
  mcausalsnps1=mcausalgenes1*mcausal_snps_per_gene # total number of causal trait 1 SNPs
  mcausalsnps2=mcausalgenes2*mcausal_snps_per_gene # total number of causal trait 2 SNPs
  R=diag(m)
  if(tolower(LD_type)=='ar') R=LD_rho^toeplitz(0:(m-1))
  if(tolower(LD_type)=='cs') R=cs(LD_rho)
  Z0_1=rep(0,m); Z0_1[causalsnpix]=sqrt(ngwas1*h21/mcausalsnps1) # true Z-stats for trait 1 for each causal gene
  Z0_2=rep(0,m); Z0_2[causalsnpix]=sqrt(ngwas2*h22/mcausalsnps2) # true Z-stats for trait 2 for each causal gene
  # moments of null and non-null test statistic distributions
  e0_1=m # null mean of trait 1
  v0_1=2*tr(R%*%R) # null variance of trait 1
  e1_1=m+c(t(Z0_1)%*%Z0_1) # non-null mean of trait 1
  v1_1=v0_1+4*c(t(Z0_1)%*%R%*%Z0_1) # non-null variance of trait 1
  e0_2=m # null mean of trait 2
  v0_2=2*tr(R%*%R) # null variance of trait 2
  e1_2=m+c(t(Z0_2)%*%Z0_2) # non-null mean of trait 2
  v1_2=v0_2+4*c(t(Z0_2)%*%R%*%Z0_2) # non-null variance of trait 2
  # distribution parameters
  beta0_1=e0_1/v0_1 # rate parameter of null distribution for trait 1
  alpha0_1=e0_1*beta0_1 # shape parameter of null distribution for trait 1
  beta1_1=e1_1/v1_1 # rate parameter of non-null distribution for trait 1
  alpha1_1=e1_1*beta1_1 # shape parameter of non-null distribution for trait 1
  beta0_2=e0_2/v0_2 # rate parameter of null distribution for trait 2
  alpha0_2=e0_2*beta0_2 # shape parameter of null distribution for trait 2
  beta1_2=e1_2/v1_2 # rate parameter of non-null distribution for trait 2
  alpha1_2=e1_2*beta1_2 # shape parameter of non-null distribution for trait 2
  # function to perform simulation - defined down here so I don't need to give it any arguments - they come from above
  simres=function(overlap_corr,verbose.=verbose) {
    SigmaOverlap=matrix(c(1,overlap_corr,overlap_corr,1),2,2) # sample overlap correlation matrix
    K=kronecker(SigmaOverlap,R) # variance-covariance matrix of Z-stats for each trait and gene
    RES1=RES2=matrix(nr=niter,nc=M) # store results in these
    numsig1=numsig2=c()
    for(i in 1:M) {
      Z1=Z2=rep(0,m)
      if(i %in% ix1) Z1[causalsnpix]=sqrt(ngwas1*h21/mcausalsnps1)
      if(i %in% ix2) Z2[causalsnpix]=sqrt(ngwas2*h22/mcausalsnps2)
      z=rmvn(niter,c(Z1,Z2),K)
      T1=rowSums(z[,1:m]^2)
      T2=rowSums(z[,-c(1:m)]^2)
      numsig1[i]=mean(pgamma(T1,shape=alpha0_1,rate=beta0_1,lower.tail=FALSE)<(0.05/12727))
      numsig2[i]=mean(pgamma(T2,shape=alpha0_2,rate=beta0_2,lower.tail=FALSE)<(0.05/12727))
      # log versions to avoid numerical instability
      ## trait 1
      log_f0_1=dgamma(T1,shape=alpha0_1,rate=beta0_1,log=TRUE)
      log_f1_1=dgamma(T1,shape=alpha1_1,rate=beta1_1,log=TRUE)
      log_ratio_1=log_f0_1-log_f1_1
      delta_ratio_1=delta1/(1-delta1)
      p1=1/(1+exp(log_ratio_1)*delta_ratio_1)
      ## trait 2
      log_f0_2=dgamma(T2,shape=alpha0_2,rate=beta0_2,log=TRUE)
      log_f1_2=dgamma(T2,shape=alpha1_2,rate=beta1_2,log=TRUE)
      log_ratio_2=log_f0_2-log_f1_2
      delta_ratio_2=delta2/(1-delta2)
      p2=1/(1+exp(log_ratio_2)*delta_ratio_2)
      RES1[,i]=p1
      RES2[,i]=p2
      # number of significant genes
      if(i %% floor(0.1*M)==0 & verbose.) cat(round(i/M*100),'% complete\n',sep='')
    }
    out=list(shared_counts=rowSums(RES1*RES2), # rows are simulations, columns are genes
             counts1=rowSums(RES1), # rows are simulations, columns are genes
             counts2=rowSums(RES2), # rows are simulations, columns are genes
             numsig1=numsig1,
             numsig2=numsig2) 
    return(out) # niter-length vector of estimated shared counts
  }
  # perform SIMEX
  corrs=seq(upsilon_overlap,sign(upsilon_overlap)*0.99,length.out=nlambdas)
  RESI=matrix(nr=niter,nc=length(corrs))
  counts1=counts2=matrix(nr=niter,nc=length(corrs))
  numsig1=numsig2=c()
  for(i in 1:length(corrs)) {
    if(verbose) cat(i,'\n')
    fitt=simres(corrs[i],verbose.=verbose)
    RESI[,i]=fitt$shared_counts
    counts1[,i]=fitt$counts1
    counts2[,i]=fitt$counts2
    numsig1[i]=sum(fitt$numsig1)
    numsig2[i]=sum(fitt$numsig2)
  }
  # end
  if(do.extrapolatation) {
    dat=cbind(corrs,colMeans(RESI))
    colnames(dat)=c('sample_overlap_corr','estimated_shared_count')
    x=dat[,1]
    y=dat[,2]
    upsilon=x[1]
    initshared=y[1]
    X=cbind(x,x^2,x^3)
    net=glmnet::cv.glmnet(X,y,alpha=1/2)
    fit=glmnet::glmnet(X,y,alpha=1/2,lambda=net$lambda.min)
    est=as.matrix(coef(fit))
    X=cbind(1,X)
    # extrapolate
    s=seq(0,1*sign(upsilon),length.out=20)
    S=cbind(1,s,s^2,s^3)
    ypred=S%*%est
    yl=sort(range(c(y,ypred)))
    dyl=diff(yl)
    if(doplot) {
      plot(x,y,type='b',
           xlab=expression(upsilon[tt]*'*'),
           ylab='est. shared genes',
           xlim=sort(c(0,sign(upsilon)*1)),ylim=yl,pch=19,col='gray90')
      points(upsilon,y[1],col='#FA7F84',pch=19)
      text(upsilon,y[1]+dyl/9,col='#FC494F',cex=4/5,srt=90,label=round(y[1]))
      text(0,ypred[1]+dyl/9,col='#5846FA',cex=4/5,srt=90,label=round(ypred[1]))
      points(0,ypred[1],col='#9C91FF',pch=19)
      lines(s,ypred,lwd=1/2)
      legend('top',legend=c('original','adjusted','simulated'),pch=c(19,19,19),
             col=c('#FA7F84','#9C91FF','gray90'),cex=4/5)
    }
    adj=ypred[1]/dat[1,2]
    out=list(
      adj_shared_count=estimated_overlap_count*adj,
      init_shared_count_est=estimated_overlap_count,
      adj=adj,
      lambda_res=x,
      RES=RESI,
      dat=dat,
      counts1=counts1, # estimated total causal gene counts for trait 1
      counts2=counts2,
      numsig1=numsig1,
      numsig2=numsig2)
  } else {
    out=list(
      counts1=counts1,
      counts2=counts2,
      numsig1=numsig1,
      numsig2=numsig2
    )
  }
  return(out)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# actually running for all trait pairs
trait_info=list(
  AD=c('AD_EUR.Rds','EUR',487511),
  AFIB=c('AFIB_EUR.Rds','EUR',2339188),
  ALS=c('ALS_EUR.Rds','EUR',138086),
  CAD=c('CAD_EUR.Rds','EUR',1165690),
  SCZ=c('SCZ_EUR.Rds','EUR',320404),
  CKD=c('CKD_EUR.Rds','EUR',625219),
  INT=c('INT_EUR.Rds','EUR',216381),
  MDD=c('MDD_EUR.Rds','EUR',807553),
  MS=c('MS_EUR.Rds','EUR',456348),
  PD=c('PD_EUR.Rds','EUR',611485),
  SLEEP=c('SLEEP_EUR.Rds','EUR',446118),
  STROKE=c('STROKE_EUR.Rds','EUR',1296908),
  T2D=c('T2D_ALL.Rds','EUR',2535601),
  ADHD=c('ADHD_EUR.Rds','EUR',225534),
  BIP=c('BIPOLAR_EUR.Rds','EUR',413466),
  BMI=c('BMI_EUR.Rds','EUR',694649),
  DBP=c('DBP_EUR.Rds','EUR',757601),
  EDU=c('EDU_EUR.Rds','EUR',3037499),
  HDL=c('HDL_EUR.Rds','EUR',1300000),
  LBD=c('LBD_EUR.Rds','EUR',16516),
  LDL=c('LDL_EUR.Rds','EUR',1300000),
  LUPUS=c('LUPUS_EUR.Rds','EUR',324698),
  PP=c('PP_EUR.Rds','EUR',757601),
  RA=c('RA_EUR.Rds','EUR',456348),
  SBP=c('SBP_EUR.Rds','EUR',757601),
  TG=c('TG_EUR.Rds','EUR',1300000),
  agingR1=c('SurrealGANR1_EUR.Rds','EUR',49482),
  agingR2=c('SurrealGANR2_EUR.Rds','EUR',49482),
  agingR3=c('SurrealGANR3_EUR.Rds','EUR',49482),
  agingR4=c('SurrealGANR4_EUR.Rds','EUR',49482),
  agingR5=c('SurrealGANR5_EUR.Rds','EUR',49482),
  mvAGING=c('mvAGING_EUR.Rds','EUR',1900000)
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# calculating adjusted shared counts all trait pairs ####
## heritabilities
setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/')
h2s=readRDS('rgs_many_traits_w_h2.Rds') %>% 
  as_tibble() %>%
  filter(!(trait1 %in% c('FCV','FCF','agingAvg','agingJoint','PTSD','COVID19')),
         !(trait2 %in% c('FCV','FCF','agingAvg','agingJoint','PTSD','COVID19'))) %>%
  mutate(trait1=ifelse(trait1=='BIPOLAR','BIP',trait1),
         trait2=ifelse(trait2=='BIPOLAR','BIP',trait2))
h2=data.frame(
  trait=c(h2s$trait1,h2s$trait2),
  h2=c(h2s$h2_trait1,h2s$h2_trait2),
  h2se=c(h2s$h2se_trait1,h2s$h2se_trait2)) %>%
  distinct(trait,.keep_all=TRUE)
## estimated counts
yin=readRDS('/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/propshared_causal_genes_all_pairs.Rds')
counts=rbind(
  yin %>% select(trait=trait1,post_count=count_trait1),
  yin %>% select(trait=trait2,post_count=count_trait2)
) %>%
  group_by(trait) %>%
  summarise(post_count=mean(post_count)) %>%
  ungroup()
# data=readRDS(
#   # '/Users/lorincn/Documents/working_space/data_while_at_ASHG/num_causal_genes_posterior_tau_delta_results.Rds'
#   '/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/num_causal_genes_posterior_tau_delta_results.Rds'
# )
# posterior_parameters=data$posterior_parameters
# posterior_hyperparameters=data$posterior_hyperparameters
# pp=posterior_parameters %>%
#   group_by(trait) %>%
#   summarise(delta=mean(mean_delta),
#             tau=mean(mean_tau))
## estimated shared gene counts
res=readRDS('/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/shared_causal_gene_counts_matrix_32_traits.Rds')
shared_est=res$shared_est
## sample overlap correlations
setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/')
Upsilon=readRDS('GWAS_sample_overlap_correlation_matrix.Rds')
# start SIMEX for each trait pair 
gg=t(combn(names(trait_info),2))
resdf=data.frame()
dats=list()
for(i in 1:nrow(gg)) {
  trait1=gg[i,1]
  trait2=gg[i,2]
  cat(trait1,', ',trait2,'\n',sep='')
  count1=counts %>% filter(trait==trait1) %>% pull(post_count)
  count2=counts %>% filter(trait==trait2) %>% pull(post_count)
  # delta1=pp %>% filter(trait==trait1) %>% pull(delta)
  # delta2=pp %>% filter(trait==trait2) %>% pull(delta)
  h2_1=h2 %>% filter(trait==trait1) %>% pull(h2)
  h2_2=h2 %>% filter(trait==trait2) %>% pull(h2)
  n1=as.numeric(trait_info[[gg[i,1]]][3])
  n2=as.numeric(trait_info[[gg[i,2]]][3])
  ctt=shared_est[which(rownames(shared_est)==trait1),which(colnames(shared_est)==trait2)]
  ctt=round(ctt)
  upsilon=Upsilon[which(rownames(Upsilon)==trait1),which(colnames(Upsilon)==trait2)]
  # fit SIMEX
  setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/results/SIMEX_plots')
  png(paste0(trait1,'_',trait2,'_sharing_from_sample_overlap.png'),width=4,height=3,units='in',res=300)
  simex=shared_count_simex(
    M=8148, # number of jointly tested genes
    Mcausal1=count1, # estimated number of causal genes for trait 1
    Mcausal2=count2, # estimated number of causal genes for trait 2
    # delta1=delta1, # proportion of non-causal trait 1 genes
    # delta2=delta2, # proportion of non-causal trait 2 genes
    ngwas1=n1, # sample size of trait 1 GWAS
    ngwas2=n2, # sample size of trait 2 GWAS
    estimated_overlap_count=ctt, # data-estimated number of overlapping causal genes
    upsilon_overlap=upsilon, # n01/sqrt(n1*n2)*Corr(t1,t2)
    # assumed simulation parameters
    h21=h2_1, # SNP heritability trait 1
    h22=h2_2, # SNP heritability trait 1
    niter=50, # number of iterations
    m=50, # number of tested SNPs per gene
    mcausal_snps_per_gene=3, # number of causal SNPs per causal gene
    LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
    LD_rho=0.5, # correlation parameter
    nlambdas=20, # number of lambdas to evaluate in SIMEX
    doplot=TRUE, # should a plot be created?
    verbose=FALSE
  )
  dev.off()
  dats[[i]]=simex$dat
  names(dats)[i]=paste0(trait1,'|',trait2)
  saveRDS(dats,'/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/results/simex_dats.Rds')
  # save results
  toadd=data.frame(
    trait1=trait1,
    trait2=trait2,
    adj_factor=simex$adj,
    adj_shared_count=simex$adj_shared_count,
    init_shared_count=simex$init_shared_count_est
  )
  resdf=rbind(resdf,toadd)
}
setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/number_causal_genes/results')
saveRDS(resdf,'shared_counts_SIMEX_adj_all_pairs.Rds')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# plotting adjusted shared counts ####
library(dplyr);library(ggplot2);library(ggrepel)
setwd('/Volumes/chengflab/Cheng-Noah/manuscripts/number_causal_genes/results')
init_res=readRDS('/Volumes/chengflab/Cheng-Noah/manuscripts/number_causal_genes/shared_causal_gene_counts_matrix_32_traits.Rds')
init_res=init_res$shared_est
adj_res=readRDS('/Volumes/chengflab/Cheng-Noah/manuscripts/number_causal_genes/adjusted_shared_causal_gene_counts_matrix_32_traits.Rds')
adj_res=adj_res[rownames(init_res),colnames(init_res)]
adj_res=adj_res*init_res
# make upper triangles NAs so I don't get duplicates
init_res[upper.tri(init_res,diag=TRUE)]=NA
adj_res[upper.tri(adj_res,diag=TRUE)]=NA
resdf=bind_cols(
  data.frame(trait1=rep(colnames(init_res),each=nrow(init_res)),
             trait2=rep(rownames(init_res),ncol(init_res)),
             init_count=c(init_res)),
  data.frame(adj_count=c(adj_res))
) %>%
  filter(trait1!=trait2) %>%
  na.omit()
# if adjusted count > initial count, just use initial count
table(resdf$init_count<resdf$adj_count)
resdf=resdf %>% mutate(adj_count=ifelse(adj_count>init_count,init_count,adj_count))
# plot differences between initial and adjusted estimates of shared counts
p1=ggplot(resdf,aes(x=init_count,y=adj_count)) +
  geom_point(color='#DBC4FF',alpha=0.5) +
  geom_point(color='#7F3BF5',
             data=resdf %>%
               filter(init_count>100) %>%
               filter(adj_count/init_count<0.975)) +
  geom_abline(intercept=0,slope=1,linetype='dashed',lwd=1/3) +
  geom_text_repel(aes(label=paste0(trait1,',',trait2)),
                  max.overlaps=100,
                  force=10,
                  size=3,
                  data=resdf %>%
                    filter(init_count>100) %>%
                    filter(adj_count/init_count<0.975)) +
  theme_bw() +
  theme(strip.background=element_rect(fill="white")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='initial shared gene count estimate',
       y='adjusted shared gene count estimate')
# plot distribution of adjustment factor
p2=resdf %>%
  filter(init_count>0) %>%
  ggplot(aes(x=adj_count/init_count)) +
  geom_histogram(fill='#DBC4FF',color='#7F3BF5',bins=50) +
  # scale_x_continuous(breaks=c(0,0.25,0.5,0.75,0.9,1),
  #                    labels=c(0,0.25,0.5,0.75,0.9,1),
  #                    limits=c(0,1)) +
  theme_bw() +
  theme(strip.background=element_rect(fill="white")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='adjusted shared count / original shared count')
p2
ggpubr::ggarrange(p1,p2,nrow=1,ncol=2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# simulation of false positive rate ####
mcausals=c(100,250,500,750)
h2s=c(0.4,0.3,0.2,0.1)
ngwas=c(1e5,2.5e5,5e5)
COUNTS1=COUNTS2=SD1=SD2=array(dim=c(length(mcausals),length(h2s),length(ngwas)))
NUMSIG1=NUMSIG2=COUNTS1
sss=0
for(i in 1:length(mcausals)) {
  for(j in 1:length(h2s)) {
    for(k in 1:length(ngwas)) {
      sss=sss+1
      cat(round(sss/prod(dim(COUNTS1))*100),'%\n',sep='')
      simex=shared_count_simex(
        M=8148, # number of jointly tested genes
        Mcausal1=mcausals[i], # estimated number of causal genes for trait 1
        Mcausal2=mcausals[length(mcausals)-i+1], # estimated number of causal genes for trait 2
        ngwas1=ngwas[k], # sample size of trait 1 GWAS
        ngwas2=ngwas[length(ngwas)-k+1], # sample size of trait 2 GWAS
        estimated_overlap_count=0, # data-estimated number of overlapping causal genes
        upsilon_overlap=0.5, # n01/sqrt(n1*n2)*Corr(t1,t2)
        # assumed simulation parameters
        h21=h2s[j], # SNP heritability trait 1
        h22=h2s[length(h2s)-j+1], # SNP heritability trait 1
        niter=20, # number of iterations
        m=50, # number of tested SNPs per gene
        mcausal_snps_per_gene=3, # number of causal SNPs per causal gene
        LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
        LD_rho=0.5, # correlation parameter
        nlambdas=2, # number of lambdas to evaluate in SIMEX
        doplot=TRUE, # should a plot be created?
        verbose=FALSE,
        do.extrapolatation=FALSE)
        # store results
        counts1=simex$counts1
        counts2=simex$counts2
        COUNTS1[i,j,k]=mean(simex$counts1)
        COUNTS2[i,j,k]=mean(simex$counts2)
        SD1[i,j,k]=sd(simex$counts1)
        SD2[i,j,k]=sd(simex$counts2)
        NUMSIG1[i,j,k]=mean(simex$numsig1)
        NUMSIG2[i,j,k]=mean(simex$numsig2)
    }
  }
}
todf=function(arr,mcausals,h2s,ngwas) {
  data.frame(
    x=c(arr),
    mcausal=rep(rep(mcausals,length(h2s)),length(ngwas)),
    h2=rep(rep(h2s,each=length(mcausals)),length(ngwas)),
    ngwas=rep(ngwas,each=length(mcausals)*length(h2s)))
}
estdf=rbind(
  todf(COUNTS1,mcausals,h2s,ngwas),
  todf(COUNTS2,rev(mcausals),rev(h2s),rev(ngwas))
)
sddf=rbind(
  todf(SD1,mcausals,h2s,ngwas),
  todf(SD2,rev(mcausals),rev(h2s),rev(ngwas))
)
numdf=rbind(
  todf(NUMSIG1,mcausals,h2s,ngwas),
  todf(NUMSIG2,rev(mcausals),rev(h2s),rev(ngwas))
)
library(dplyr);library(ggplot2)
plotdf=bind_cols(
  estdf,
  sddf %>% dplyr::select(se=x),
  numdf %>% dplyr::select(nsig=x))
plotdf$ngwas=ifelse(plotdf$ngwas<1e6,
  paste0(round(plotdf$ngwas/1e3),'K'),
  paste0(round(plotdf$ngwas/1e6),'M'))
plotdf$h2=paste0('SNP h2=',round(plotdf$h2,1))
p1=ggplot(plotdf,aes(x=mcausal,y=x,color=factor(ngwas))) +
  geom_abline(intercept=0,slope=1,lwd=1/3,color='gray80') +
  geom_line(aes(group=factor(ngwas))) +
  geom_line(aes(x=mcausal,y=x+2*se,color=factor(ngwas),group=factor(ngwas)),lty=2) +
  geom_line(aes(x=mcausal,y=x-2*se,color=factor(ngwas),group=factor(ngwas)),lty=2) +
  geom_point() +
  scale_color_manual('GWAS n',
                     values=ggpubfigs::friendly_pal('contrast_three')) +
  facet_wrap(~h2,nrow=1) +
  theme_bw() +
  scale_x_continuous(breaks=mcausals,labels=mcausals) +
  scale_y_continuous(breaks=mcausals,labels=mcausals) +
  theme(strip.background=element_rect(fill="gray90")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.position='bottom') +
  labs(x='true number of causal genes',
       y='estimated number of \ncausal genes')

p2=ggplot(plotdf,aes(x=mcausal,y=nsig,color=factor(ngwas))) +
  geom_abline(intercept=0,slope=1,lwd=1/3,color='gray80') +
  geom_line(aes(group=factor(ngwas))) +
  geom_point() +
  scale_color_manual('GWAS n',
                     values=ggpubfigs::friendly_pal('contrast_three')) +
  facet_wrap(~h2,nrow=1) +
  theme_bw() +
  scale_x_continuous(breaks=mcausals,labels=mcausals) +
  scale_y_continuous(breaks=mcausals,labels=mcausals) +
  theme(strip.background=element_rect(fill="gray90")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.position='bottom') +
  labs(x='true number of causal genes',
       y='number of Bonferroni-\nsignificant genes')
ggpubr::ggarrange(p1,p2,nrow=2,ncol=1,common.legend=TRUE)
setwd('/Volumes/chengflab/Cheng-Noah/manuscripts/number_causal_genes/results')
ggsave('estimated_numcausal_simulation_performance.pdf',
       width=7,height=5,dpi=300)

rm(list = ls())  
library(dplyr)
library(MRPRESSO)
library(TwoSampleMR)
library(dplyr)
library(ggsci)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggview)
library(plotly)
library(TSMRhelper)
library(forestploter)
library(grid)
library(stringr)
library(grDevices)
library(parallel)
library(vroom)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')


# 读取暴露
exp_clumped<- vroom("338_all_data.csv")

outcome_dat<-extract_outcome_data(snps = exp_clumped$SNP,outcomes = "ieu-b-7",proxies = F)

# 合并数据
dat <- harmonise_data(exposure_dat = exp_clumped,outcome_dat = outcome_dat)
dat <- distinct(dat, SNP, .keep_all = TRUE)
dat<-dat[dat$mr_keep,]
write.csv(dat,file = "338_PD_dat.csv",row.names = F)

# 分组
dat<-split(dat,dat$id.exposure)
# 跑mr
library(parallel)
cl <- makeCluster(detectCores())
res<-parLapply(cl,dat,mr)
stopCluster(cl)
res<-do.call(rbind,res)
res<-generate_odds_ratios(res)
write.csv(res,file = "338_PD_res.csv",row.names = F)

# 异质性检验
dat<-do.call(rbind,dat)
heterogeneity<-mr_heterogeneity(dat)
write.csv(heterogeneity,file = "338_PD_heterogeneity.csv",row.names = F)

# 多效性检验
pleiotropy_test<-mr_pleiotropy_test(dat)
write.csv(pleiotropy_test,file = "338_PD_pleiotropy_test.csv",row.names = F)

# 绘图
res_id<-read.csv("./338_PD_res.csv")%>%get_sbeta_res()%>%dplyr::filter(method=="Inverse variance weighted",pval<0.05)%>%pull(id.exposure)
dat_total<-read.csv("./338_PD_dat.csv")
for (idr in res_id) {
  dat<-subset(dat_total,id.exposure==idr)
  scatter_plot<-mr_scatter_plot(mr(dat),dat)[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme(axis.title.y = element_text(size = 20))+
    theme_bw(base_size = 16)+ theme(
      plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")
    )
  ggsave(plot=scatter_plot,filename = paste0(idr,"_scatter_plot_plot.pdf"),device = "pdf",width = 10,height = 10)
  funnel_plot<-mr_funnel_plot(mr_singlesnp(dat,all_method=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_simple_mode","mr_weighted_mode")))[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme(axis.title.y = element_text(size = 20))+
    theme_bw(base_size = 16)+ theme(
      plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")
    )
  ggsave(plot=funnel_plot,filename = paste0(idr,"_funnel_plot_plot.pdf"),device = "pdf",width = 10,height = 10)
  forest_plot<-mr_forest_plot(mr_singlesnp(dat))[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme_bw()+
    theme(legend.position = 'none',plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"))
  ggsave(plot=forest_plot,filename = paste0(idr,"_forest_plot_plot.pdf"),device = "pdf",width = 10,height = 10)
  
  leaveoneout_plot<-mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme_bw()+
    theme(legend.position = 'none',plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),
          axis.title.x = element_text(size = 12))
  ggsave(plot=leaveoneout_plot,filename = paste0(idr,"_leaveoneout_plot_plot.pdf"),device = "pdf",width = 10,height = 10)  
}


# 用空格调整列宽
dt$` ` <- paste(rep(" ", 30), collapse = " ")

dt$pval<-round(dt$pval,digits = 3)
dt$pval<-ifelse(dt$pval<0.001,"<0.001",dt$pval)

colnames(dt)[5]<-"italic(P)*-Value"

dt<-dt%>%dplyr::rename(Outcome=outcome,
                       `Trait`=exposure,
                       Method=method,
                       nSNP=nsnp,)

dt$Method<-ifelse(dt$Method=="Inverse variance weighted","IVW",dt$Method)
dt<-subset(dt, Method=="IVW")
tm <- forest_theme(base_size = 8,colhead=list(fg_params = list(parse=TRUE)))
p=forest(dt[,c(2:5,10,9)],#按顺序选择数据1-4列、12-13列为绘图区域和OR列的位置（图放中间）、8-11列作为森林图元素（Q、Q_pval和多效性P值）
         est = dt$or,
         lower = dt$or_lci95,
         upper = dt$or_uci95,
         sizes = 0.6,
         ci_column = 5,
         ref_line = 1,
         xlim = c(0,3),
         ticks_at = c(0,1,2,3),
         theme = tm)

ggsave("338_PD_森林图.pdf",plot=p,device = "pdf",width = 20,height = 20)

#处理好之后结果放到一个文件夹里面
id_473_Uc<-as.data.frame(res_id)   
setnames(id_473_Uc, old = names(id_473_Uc)[1], new = "id")
write.csv(id_473_Uc,"id_338_PD.csv",row.names = F)

####### 反向mr######
exp_clumped<-extract_instruments("ieu-b-7",p1=5e-8)#5e-6 1e-5
table(!((((exp_clumped$beta.exposure)^2)/((exp_clumped$se.exposure)^2))<10))
exp_clumped <- exp_clumped[!((((exp_clumped$beta.exposure)^2)/((exp_clumped$se.exposure)^2))<10),]


# 结局
# 读取Trait.csv文件
idTotrait <- read.csv("./338_Trait.csv")

# 读取id_473_Uc.csv文件
total_data <- read.csv("id_338_PD.csv")
id_list <- total_data$id

# 设置文件夹路径 原始数据路径
folder_path <- "F:/脑脊液转SNP原始数据"

# 构建文件路径并添加后缀
file_paths <- file.path(folder_path, paste0(id_list, "_buildGRCh37.tsv.csv"))

dat <- list()

for (i in 1:length(id_list)) {
  print(id_list[i])
  file_path <- file_paths[i]
  outcome_dat <- vroom(file_path)
  outcome_dat$p_value <- as.numeric(outcome_dat$p_value)
  # 筛选snp
  #colnames(outcome_dat)
  outcome_dat <- outcome_dat %>% filter(SNP %in% unique(exp_clumped$SNP))
  outcome_dat$id <- str_split(id_list[i], "_buildGRCh37.tsv.csv", simplify = TRUE)[, 1]
  outcome_dat$phe <- idTotrait[idTotrait$id == outcome_dat$id[1], ]$Trait
  # 转换格式
  outcome_dat <- format_data(outcome_dat,
                             type = "outcome",
                             snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "standard_error",
                             eaf_col = "effect_allele_frequency",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "p_value",
                             chr_col = "CHR",
                             pos_col = "BP",
                             id_col = "id",
                             phenotype_col = "phe")
  dat[[i]] <- harmonise_data(exposure_dat = exp_clumped, outcome_dat = outcome_dat)
}
dat<-do.call(rbind,dat)

dat<-dat[dat$mr_keep,]
write.csv(dat,file = "PD_338_dat.csv",row.names = F)

re_res<-mr(dat)
re_res<-generate_odds_ratios(re_res)
write.csv(re_res,file = "PD_338_res.csv",row.names = F)

# 异质性检验
heterogeneity<-mr_heterogeneity(dat)
write.csv(heterogeneity,file = "PD_338_heterogeneity.csv",row.names = F)

# 多效性检验
pleiotropy_test<-mr_pleiotropy_test(dat)
write.csv(pleiotropy_test,file = "PD_338_pleiotropy_test.csv",row.names = F)

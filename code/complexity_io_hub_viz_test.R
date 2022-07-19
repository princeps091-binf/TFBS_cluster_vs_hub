library(vroom)
library(GenomicRanges)
library(tidyverse)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
########################################################################################################
# Utils. Fn
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

########################################################################################################
# Files
TFBS_cluster_file<-"./data/wgEncodeUWDukeDnaseHMEC.fdr01.TFCluster.hg19.peak"
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/H1_union_trans_res_dagger_tbl.Rda"
spec_res_folder<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"
########################################################################################################
TFBS_cluster_tbl<-vroom(TFBS_cluster_file,col_names = F) %>% 
  mutate(ID=paste(X1,X2,X3,sep="_"))

TFBS_cluster_tbl %>% 
  ggplot(.,aes(X4))+
  geom_density()+
  scale_x_log10()
###################################
hub_tbl<-get_tbl_in_fn(hub_file)

hub_tbl<-do.call(bind_rows,map(unique(hub_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_folder,chromo,"_spec_res.Rda"))
  tmp_tbl<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>% 
    mutate(bins=map(bins,as.numeric)) 
  
})) %>% 
  dplyr::select(chr,node,res,bins)

plan(multisession,workers=4)
hub_tbl<-hub_tbl %>% 
  mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    inter_cl_Grange<-   GRanges(seqnames=chr,
                                ranges = IRanges(start=as.numeric(bins),
                                                 end=as.numeric(bins) + res_num[res]-1
                                ))
    inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
    return(inter_cl_Grange)
    
  }))
plan(sequential)

tmp_l<-hub_tbl %>% 
  #  filter(res%in% c("1Mb","500kb","100kb")) %>% 
  #  filter(res%in% c("10kb","50kb","5kb")) %>% 
  dplyr::select(GRange) %>% as.list

cl_GRange<-IRanges::reduce(do.call("c",tmp_l$GRange))
TFBS_cluster_Grange<-   GRanges(seqnames=TFBS_cluster_tbl$X1,
                            ranges = IRanges(start=TFBS_cluster_tbl$X2,
                                             end=TFBS_cluster_tbl$X3
                            ))
TFBS_cluster_Grange<-TFBS_cluster_Grange[seqnames(TFBS_cluster_Grange) %in% unique(seqnames(cl_GRange))]

findOverlaps(TFBS_cluster_Grange,cl_GRange) %>% 
  as_tibble %>% 
  distinct(queryHits)

in_cl_peak<-TFBS_cluster_tbl$ID[unique(queryHits(findOverlaps(TFBS_cluster_Grange,cl_GRange)))]

gg_tmp<-TFBS_cluster_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  ggplot(.,aes(X4,color=hub.io))+
  scale_color_brewer(palette="Set1")+
  theme_minimal()+
  geom_density()+
  scale_x_log10()
gg_tmp

in_vec<-TFBS_cluster_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  filter(hub.io=="in") %>% 
  dplyr::select(X4) %>% 
  unlist

out_vec<-TFBS_cluster_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  filter(hub.io=="out") %>% 
  dplyr::select(X4) %>% 
  unlist

wilcox.test(in_vec,out_vec,alternative = "greater")

library(vroom)
library(GenomicRanges)
library(tidyverse)
library(furrr)
library(parallel)
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

cage_tbl_coord_build_fn<-function(cage_tbl,ID_col){
  
  cage_coord<-cage_tbl %>% dplyr::select(ID_col) %>% unlist
  cage_start<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split=','),'[',1)),split = '\\.'),'[',1)))
  cage_end<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split=','),'[',1)),split = '\\.'),'[',3)))
  cage_chr<-unlist(lapply(strsplit(cage_coord,split = ':'),'[',1))
  return(cage_tbl%>%mutate(chr=cage_chr,start=cage_start,end=cage_end))
  
}

cage_mean_compute_fn<-function(cage_tbl){
  print("compute m")
  cl<-makeCluster(5)
  tmp_m<-parallel::parApply(cl,X = as.matrix(cage_tbl[,-1]),MARGIN = 1,function(x){
    mean(x)
  })
  stopCluster(cl)
  rm(cl)
  cage_tbl<-cage_tbl%>%
    mutate(score=tmp_m)%>%
    filter(score>0)
  return(cage_tbl)
  
}

produce_GRange_fn<-function(tmp_tbl){
  cage_Grange<-GRanges(seqnames=tmp_tbl$chrom,
                       ranges = IRanges(start=tmp_tbl$start,
                                        end=tmp_tbl$end
                       ))
  
  return(cage_Grange)
}

########################################################################################################
# Files
TFBS_cluster_file<-"./data/wgEncodeUWDukeDnaseHMEC.fdr01.TFCluster.hg19.peak"
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/HMEC_union_trans_res_dagger_tbl.Rda"
spec_res_folder<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

## Cell-line sample ID
H1<-c("CNhs14067","CNhs14068","CNhs13964")
GM12878<-c('CNhs12331','CNhs12332','CNhs12333')
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')
#TSS
cage_tss_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#',col_select = contains(c("Annotation",H1,HMEC,GM12878)))
#Enhancers
cage_enh_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",comment = '#',col_select = contains(c("Id",H1,HMEC,GM12878)))
########################################################################################################
#Enhancers
tmp_enh_tbl<-cage_enh_tbl %>%
  dplyr::select(c("Id",all_of(HMEC))) %>% 
  filter(if_any(where(is.numeric), ~ .x > 0)) %>% 
  cage_mean_compute_fn(.) %>% 
  mutate(chrom=str_split_fixed(Id,":|-",3)[,1],
         start=str_split_fixed(Id,":|-",3)[,2],
         end=str_split_fixed(Id,":|-",3)[,3]) %>% 
  filter(!(is.na(chrom))) %>% 
  dplyr::select(chrom,start,end,score) %>% 
  mutate(start=as.integer(start),
         end=as.integer(end))
#Tss
tmp_tss_tbl<-cage_tss_tbl %>%
  dplyr::select(contains(c("Annotation",all_of(HMEC)))) %>% 
  filter(if_any(where(is.numeric), ~ .x > 0)) %>% 
  cage_mean_compute_fn(.) %>% 
  mutate(chrom=str_split_fixed(`00Annotation`,":|,|\\.\\.",4)[,1],
         start=str_split_fixed(`00Annotation`,":|,|\\.\\.",4)[,2],
         end=str_split_fixed(`00Annotation`,":|,|\\.\\.",4)[,3]) %>% 
  filter(!(grepl("STAT",`00Annotation`))) %>% 
  dplyr::select(chrom,start,end,score) %>% 
  mutate(start=as.integer(start),
         end=as.integer(end))

TFBS_cluster_tbl<-vroom(TFBS_cluster_file,col_names = F) %>% 
  mutate(ID=paste(X1,X2,X3,sep="_")) %>% 
  dplyr::rename(chrom=X1,start=X2,end=X3,score=X4)


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



TFBS_cluster_GRange<-TFBS_cluster_tbl %>% produce_GRange_fn
tmp_tss_GRange<-tmp_tss_tbl %>% produce_GRange_fn
tmp_enh_GRange<-tmp_enh_tbl %>% produce_GRange_fn

findOverlaps(tmp_tss_GRange,TFBS_cluster_GRange) %>% 
  as_tibble %>% distinct(queryHits) %>% summarise(n=n()/length(tmp_tss_GRange))

findOverlaps(tmp_enh_GRange,TFBS_cluster_GRange)%>% 
  as_tibble %>% distinct(queryHits) %>% summarise(n=n()/length(tmp_enh_GRange))

findOverlaps(tmp_tss_GRange,cl_GRange) %>% 
  as_tibble %>% distinct(queryHits) %>% summarise(n=n()/length(tmp_tss_GRange))

findOverlaps(tmp_enh_GRange,cl_GRange)%>% 
  as_tibble %>% distinct(queryHits) %>% summarise(n=n()/length(tmp_enh_GRange))


###  load necessary packages   ################################
#library(phyloseq)
#library('taxize')
#library(BiocInstaller)
#library(seqinr)
#library(ape)
library(RColorBrewer)
library(rdrop2)
#library(DECIPHER)
library(vegan)
#library(stats)
#library(corrplot)
library(fields)
library(psych)
library(tidyverse) 
###############################################################################
# function to remove rows with n number of NA's
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
###############################################################################
###############################################################################
###set directory and variables for specific computer Dropbox
dir<-"/Users/geraldn/Dropbox/"
out_file<-"Documents/KAUST/eDNA/R/pipe_summary"
stud_pat<-"Dammam"  # matches study specific title from pipe (begining of files).
## set some other universal variables
primers<-c('18s_stoeck','18smini',"euka02","co1","vert")
minboots<-c(50,70,90)
min_samples_include<-2  # used duing filtering SV's in < this number will be removed
bl_prop_cut<-0.001  ## remove SV's with greater than this proportion in max of blanks
sam_sum_prop_cut<-0.001 ## remove SV's with less than this propotion of summed reads of samples
###############################################################################
#  get sample data    and   out put from dada2
# sample data need "type" column with..... "sample" "extraction blank" "pcr blank" "positive control"
###############################################################################
# sample data
sam<-openxlsx::read.xlsx(paste(dir,"Documents/KAUST/eDNA/Samples_Data/Extraction test aug18/EXPERIMENT DESIGN 2018 PCR TEMPLATE.xlsx", sep=""),
                            startRow=1, sheet=1, check.names=T)
#  unique(sam$type)
# summary, seqtab and taxass    list.files()
setwd(paste(dir,out_file,sep=""))
f_tax<-list.files(pattern = paste(stud_pat,".*taxass.*\\.csv",sep="")) ## 3 x number of primers
f_sum<-list.files(pattern = paste(stud_pat,".*summary\\.csv", sep=""))  # one per primer
f_rds<-list.files(pattern = paste(stud_pat,".*seqtab\\.rds", sep=""))
##   import all data into lists            names(sum)
tax_list<-setNames(lapply(f_tax, read.csv), tools::file_path_sans_ext(basename(f_tax)))  # mes1<-tax_list[[1]]
##  !! make sure first colomn is sequences names seq  !!!!
sum<-setNames(lapply(f_sum, read.csv), tools::file_path_sans_ext(basename(f_sum)))
sv_list<-setNames(lapply(f_rds, readRDS), tools::file_path_sans_ext(basename(f_rds))) # names(sv_list) mes<-sv_list[[1]]

#######################################
# tidy     !  summary tables   !    table from dada2 pipe
## add column for primer and sampple_id
for (i in seq_along(sum)){
  sum[[i]]["primer"] <- primers[i]
  sum[[i]]["sample_id"] <- rownames(sv_list[[i]]) }  #  
sum2<-bind_rows(sum, .id = 'names(sum)')

#######################################
# tidy sample dataframe  combine with summary
#   names(sam)  unique(sam$primer)  unique(sum2$primer)
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##    !!!!  TODO  !!???????    , join multiple sample sheets, make sample/primer unique column
##     issue becaue column name is X_number and this needs to match later on  ?????????????
###  solution???  make list for each sample sheet/primer combo
sam2<-sam %>% 
  mutate(sample_id=as.character(sample_id)) %>% 
  left_join(sum2[,c(8,9,7)]) %>% 
  dplyr::rename(reads_per_sample=nonchim)   ### reads after DADA2 filters
# limit samples if has more than one primer per miseq run  !!!!!!!!!!!!!!!!!!!!!
#    +++++++++++    will need to change for other studies
# set sample list for each primer
sam_prim1<-sam2[sam2$primer==primers[1],] 
sam_prim2<-sam2[sam2$primer==primers[2],] 
sam_list<-list(sam_prim1=sam_prim1,sam_prim2=sam_prim2) # names(sam_list)
#  get index to call column names of sample types   mes<-sam_list2[[1]]
sam_list<- lapply(sam_list, function(x) { 
  row_number<-c(1:length(as.numeric(row.names(x)))) # make index for choose columns of SVs
  x<-cbind(x,row_number)    })
# replicate to match list with primers and nboots.
 sam_list<-rep(sam_list, each=length(minboots))
## make lists for each sample type, needed for filtering  x<-sample_index_list[[1]]
sample_index_list<- lapply(sam_list, function(x) { 
  x<-x$row_number[x$type=="sample"]     })
#
blank_index_list<- lapply(sam_list, function(x) { 
  x<-x$row_number[x$type=="extraction blank" | x$type=="pcr blank"]     })
#
positive_index_list<- lapply(sam_list, function(x) { 
  x<-x$row_number[x$type=="positive control" | x$type=="mock"]     })


### merge otu and matadata ######################################
#  make elements into dataframe
dat<-purrr::map(sv_list, data.frame)
# make rownames a column
dat<-purrr::map(dat, function(x) rownames_to_column(x, var="sample_id"))  #  mes<- names(dat[[1]]) names(dat)

################################################################################################################
##     transpose seq table and merge to tax table     ### used for heatmap
dat1<- lapply(dat, function(x) { 
        s<-names(x)
        id<-x[,1]
        x<-t(x[,-1])
        x<-unname(x)
        x<-data.frame(x)
        #colnames(x)<-paste("X",id, sep="") ;return(x)
        x$seq<-s[-1]  ;return(x)
        })                                                      #   x<-dat1[[1]]  x2<-dat1[[2]] names(dat1)
# assign colnames - which are sample_id
for (i in seq_along(dat1)){
  colnames(dat1[[i]]) <- c(paste("X_",dat[[i]][,1],sep=""),"seq")
}
########  join sample/otu with taxon  ##############
##   mes5<-tax_list[[1]]  names(dat3)   names(tax_list)     mes<-dat3[[2]]  names(mes)
## duplicate dat1 to match tax_list    names(dat3)
dat2<-rep(dat1, each=3)
minboot_seq<-rep(minboots, times=length(primers))
primers_seq<-rep(primers, each=length(minboots))
############  add with tax_list- for each taxonomy
dat3<-mapply(cbind,dat2,tax_list, SIMPLIFY=F)   # names(dat3)
## add column for minboot
for (i in seq_along(dat3)){
    dat3[[i]]["primer"] <- primers_seq[i] 
    dat3[[i]]["minboot"] <- minboot_seq[i] 
    dat3[[i]]["seq.1"] <- c(1:length(row.names(dat3[[i]])))
}

####################################################################################
####    begin filter SV based on # samples and blanks
####################################################################################
    #  mes1<-dat3[[2]]  mes2<-dat3[[3]]   # names(dat3)  mes0<-tax_list[[2]]

#### remove SV's in < number of samples  - lapply through all lists
      #  hist(apply(mes1[,samples_index],1,function(x) length(which(x>0)))) # sum SV's per sample
      #  mes2<-mes1[apply(mes1[,samples_index],1,function(x) (length(which(x>0)))>=min_samples_include),]
dat4<-dat3
for (i in seq_along(dat4)){   
        mes1<-dat4[[i]]
        mes2<-mes1[apply(mes1[,sample_index_list[[i]]],1,function(x) (length(which(x>0)))>=min_samples_include),]

        dat4[[i]]<-mes2
      }     #  mes2<-dat4[[3]]  mes1<-dat3[[3]]
#### remove SV's based on blanks - lapply through all lists 
 
 dat_filt<-dat4   # make copy for next step   names(dat4)
 for (i in seq_along(dat_filt)){    #  i<-5
     mes1<-dat_filt[[i]] # mes1<-dat4[[5]]     mes1<-dat_filt[[5]]
     bl_max<- apply(mes1[,blank_index_list[[i]]], 1, max)  # get max reads for each SV from all blanks
     mes1$bl_max_prop<-bl_max/sum(bl_max)  #     sum(bl_max_prop)       hist(log(bl_max))  order(-bl_max_percent)
     #bl_mean<- apply(mes1[,blank_index], 1, mean)  # get max reads for each SV from all blanks
     #sam_max<- apply(mes1[,sample_index], 1, max) # get max reads for each SV from samples
     #sam_mean<- apply(mes1[,sample_index], 1, mean) # get max reads for each SV from samples
     #bl_per_max_sam<-bl_max/sam_max  # percent of max blank to max read of single samples
     #bl_per_sum<-bl_max/sam_sum   # percent of max blank to sum of smaples
     
     sam_sum<- apply(mes1[,sample_index_list[[i]]], 1, sum) # get sum reads for each SV from samples
     mes1$sam_sum_prop<-sam_sum/sum(sam_sum)  # % of above   hist(log(sam_sum))
     
     mes2<-mes1[mes1$bl_max_prop < bl_prop_cut,] # remove SV's with this propotion in max blanks
     
     mes3<-mes2[mes2$sam_sum_prop > sam_sum_prop_cut,] # remove based on % reads in samples- consider 0.0001
     mes3<-select(mes3,-sam_sum_prop,-bl_max_prop)
     dat_filt[[i]]<-mes3
 }  
#  names(dat_filt)   mes<-dat_filt[[2]]  q<-mes1[is.na(mes1$primer),]
###################################################################################
###################################################################################
####################################################################################
# sum all reads for each unique taxon and !   gather !!
 dat_filt_taxa<-dat_filt
for (i in seq_along(dat_filt_taxa)){
  dat_filt_taxa[[i]]<-dat_filt_taxa[[i]] %>% 
   # mes1<-mes %>% 
    select(-seq.1,-seq) %>% 
    group_by_at(vars(Superkingdom:minboot)) %>% 
    summarise_all(funs(sum)) %>% 
    gather(key="sample_id",value="reads", starts_with("X"))  %>% 
    mutate(sample_id=gsub("X_","",sample_id))
    }   # names(dat3)    mes<-dat_filt[[5]] names(dat_filt)  mes<-x$seq[x$seq.1 %in% c(1,8,2)]

### simplify list into one dataframe and tidy    ##################################
 dat_filt_taxa_df<-do.call(rbind,dat_filt_taxa)
 sam3<-select(sam2,-primer)# get rid of primer column in sample
### merge with sample data
 dat_filt_taxa_df<-dat_filt_taxa_df %>% 
  ungroup() %>% 
  left_join(sam3) %>% 
  mutate(sample_id_u=paste(primer,sample_id,minboot,sep="_")) %>% 
  mutate(miseq_id=paste(primer,minboot,sep="_")) %>% 
  unite(lineage, Superkingdom:Species, sep=";",remove=F) %>% 
  mutate(lineage2=make.names(lineage)) %>% 
   group_by(sample_id_u) %>% 
   mutate(sample_sum_after_filt=sum(reads)) %>%  # get sum reads per sample after filter
   ungroup() 
 #  unique(dat_filt_taxa_df$primer)
#   mes<-dat_filt_taxa_df[dat_filt_taxa_df$miseq_id=="NA_NA",]  # unique(dat_filt_taxa_df$miseq_id)

#####################################################################################################
#####################################################################################################
###  rarify  and simpson diversity ???    names(rar)    x<-rar[[2]]     x2<-rar2[[2]]   
rar<-dat_filt_taxa_df %>% 
    select(miseq_id,sample_id,lineage2,reads) %>% 
    split(list(.$miseq_id))
## rarify
rar2<-lapply(rar, function(x) { 
  x1<-select(x,sample_id,lineage2, reads)
  x1<-spread(x1, key=lineage2,value=reads) 
  
  nn<-cbind(x1$sample_id,rowSums(x1[,-1])) ## use only sample not blanks to get miniumum smaple size for rarefy !!!!!!
  nn<-data.frame(nn)
  colnames(nn) <- c("sample_id", "r_sum")
  nn$r_sum<-as.numeric(as.character(nn$r_sum))
  nn<-left_join(nn,sam3)
  nn<-nn[nn$type=="sample",]
  s<-min(nn$r_sum)  
  
  x2<-data.frame(vegan::rrarefy(x1[,-1],sample=s))
  x3<-bind_cols(x1[,1],x2)
  x4<-gather(x3, key="lineage2",value="rare_reads", -sample_id) # 
  x<-x4
})
rar3<-bind_rows(rar2, .id = "miseq_id")  # combine rarified lists
## diversity
div<-lapply(rar2, function(x) { 
  #  x<-rar2[[1]]
  x1<-select(x,sample_id,lineage2, rare_reads)
  x1<-spread(x1, key=lineage2,value=rare_reads) 
  s<-min(rowSums(x1[,-1]))
  x2<-data.frame(vegan::diversity(x1[,-1],index = "shannon"))
  x3<-bind_cols(x1[,1],x2)
  names(x3)[2]<-"shan_div_rare" # 
  x<-x3
})
div2<-bind_rows(div, .id = "miseq_id")

##  combine rare data     names(dat6)
dat_filt_taxa_df<-dat_filt_taxa_df %>% 
  left_join(rar3) %>%  ## unique(dat5$miseq_id)  unique(rar3$miseq_id)
  left_join(div2) %>%
  select(-lineage2)

### richness rarefied richness
rich<-dat_filt_taxa_df %>% 
  select(miseq_id,sample_id_u,lineage,reads,rare_reads) %>% 
  group_by(miseq_id,sample_id_u) %>% 
  summarise(rich=length(rare_reads[reads>0]),rich_rare=length(reads[rare_reads>0]))

dat_filt_taxa_df<-dat_filt_taxa_df %>% 
  left_join(rich)
# remove taxonomy with <= n number of na's, bigger number more kept   ???
#dat5<-delete.na(dat4, 4)      unique(dat)
##########################################
# export  to project folder
out_file<-"Documents/KAUST/eDNA/R/csv/"
data.table::fwrite(dat_filt_taxa_df,paste(dir,out_file,"extr_test_all_reads.csv",sep=""),row.names=F, sep=",")










#####################################################################################################
#####################################################################################################
#####################################################################################################
#  old code from original quant paper
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
## remove SVs in only one sample and how many   #########################################################
d<-c(9:104)
hist(apply(ts5[,d],1,function(x) length(which(x>0))))
ts5g1<-ts5[apply(ts5[,d],1,function(x) (length(which(x>0)))>1),] ##92 of 133 more than 1 sample
ts5g2<-ts5[apply(ts5[,d],1,function(x) (length(which(x>0)))>2),] ##82 of 133 more than 2 sample
ts5g3<-ts5[apply(ts5[,d],1,function(x) (length(which(x>0)))>3),] ##73 of 133 more than 3 sample

ukki<-length(which(ts5less1$Kingdom=="unknown"))  #unknown meta
ukph<-length(which(ts5less1$Phylum=="unknown"))-ukki  #unknown meta
ukcl<-length(which(ts5less1$Class=="unknown"))-ukph #unknown cordates
ukk<-length(rownames(ts5less1))-sum(ukki,ukph,ukcl)
ukkk<-c(ukki,ukph,ukcl,ukk)
rdf<-data.frame(cbind(c("unknown","unknown meta","unknown cordates","fish"),ukkk))
#
ukki2<-length(which(ts5less2$Kingdom=="unknown"))  #unknown meta
ukph2<-length(which(ts5less2$Phylum=="unknown"))-ukki  #unknown meta
ukcl2<-length(which(ts5less2$Class=="unknown"))-ukph #unknown cordates
ukk2<-length(rownames(ts5less2))-sum(ukki2,ukph2,ukcl2)
ukkk2<-c(ukki2,ukph2,ukcl2,ukk2)
#
ukki3<-length(which(ts5less3$Kingdom=="unknown"))  #unknown meta
ukph3<-length(which(ts5less3$Phylum=="unknown"))-ukki  #unknown meta
ukcl3<-length(which(ts5less3$Class=="unknown"))-ukph #unknown cordates
ukk3<-length(rownames(ts5less3))-sum(ukki2,ukph2,ukcl2)
ukkk3<-c(ukki3,ukph3,ukcl3,ukk3)
##   table of how many sv wre removed based on number of samples they are in - 1,2 and 3 
rdf<-cbind(rdf,ukkk2,ukkk3)
###################################################################################
## Determine what SV's to remove based on blank and total reads #########################
### use ts5g1 92 SV's - removed SV's with only 1 occurance
gc<-c(9:104)  # colnames(ts5g1)  all samples
gs<-c(9:99,103:104) ##  non blank samples
bl<-c(100:102)   ##Blanks
s<-colSums(ts5g1[,gc])
s<-s[order(s)]    ##   hist(s)
tot_reads<-sum(ts5g1[,gc])
reads_in_tbl<-sum(ts5g1[,c(100,101)]) ##total reads in blanks
co<-tot_reads*.0001  ## 0.01 percent of reads
co<-100
blrow<-rep(1:length(rownames(ts5g1)),times=length(bl))
blcol<-rep(bl,each=length(rownames(ts5g1)))
svtotr<-rowSums(ts5g1[,gc]) ### total reads for each SV
ch<-ts5g1[rowSums(ts5g1[,gc])<co,]   ## SV table with SV removed with less than co


j<-102
i<-14
## set up balnk tables
numgt<-seq(1:(length(rownames(ts5g1))*length(bl)))
bl_reads<-seq(1:(length(rownames(ts5g1))*length(bl)))
per_reads<-seq(1:(length(rownames(ts5g1))*length(bl)))
numz<-seq(1:(length(rownames(ts5g1))*length(bl)))
max_read<-seq(1:(length(rownames(ts5g1))*length(bl)))

## for each blank go through each sample get samples info based on blank and summarize --xx3
k<-1
for (j in bl) {
  for (i in 1:length(rownames(ts5g1))) {
    numgt[k]<-length(which(ts5g1[i,gs]>ts5g1[i,j]))  #num samples greater than blank  (gs nonblank smaples)
    bl_reads[k]<-ts5g1[i,j]
    per_reads[k]<-ts5g1[i,j]/sum(ts5g1[i,gc])
    numz[k]<-length(which(ts5g1[i,gs]==0))  ## number of zeros
    max_read[k]<-max(ts5g1[i,gs])  ## max read per sample
      k<-k+1
  }
}
numgt<-length(gs)-numgt  # now big numbers are bad - blank greater than most smaples
  xx3<-data.frame(cbind(numgt,bl_reads,per_reads,max_read,numz,blrow,blcol))
  xx3$numofbl_grsam_no_zero<-xx3$numgt-xx3$numz  ## hist(xx3$numofbl_grsam_no_zero)
  xx3$bl_per_max<-bl_reads/max_read ## percent of max num read that is blank
  xx3$SVtotread<-rep(svtotr, times=3)
  xx3$bl_per_tot<-bl_reads/xx3$SVtotread
  xx44<-xx3[which(xx3$bl_per_max > 0.1)  ,]
  # removal criteria    ##
 # remov<-xx3[which(xx3$bl_per_max > 0.1 | xx3$bl_per_tot > 0.01 | xx3$SVtotread < co)  ,]
  remov<-xx3[which(xx3$bl_per_max > 0.1 | xx3$bl_per_tot > 0.01)  ,] ## find SV's that blank is greater >10% of max read or blank is >1% of total reads
  remov2<-remov[!duplicated(remov$blrow),]
  ## final new #######################################
  ts5good<-ts5g1[-remov2$blrow,]    #############   50 SV,s #################
  
  #########################################################################################################
  ####################    inserted 19/ mar from below #############################################################################
  ## convert NA to unknown for taxonomy assignment ##########################
  #
  ts5good[] <- lapply(ts5good, function(x){
    # check if you have a factor first:
    if(!is.factor(x)) return(x)
    # otherwise include NAs into factor levels and change factor levels:
    x <- factor(x, exclude=NULL)
    levels(x)[is.na(levels(x))] <- "unknown"
    return(x)
  })
  ts8[] <- lapply(ts8, function(x){
    # check if you have a factor first:
    if(!is.factor(x)) return(x)
    # otherwise include NAs into factor levels and change factor levels:
    x <- factor(x, exclude=NULL)
    levels(x)[is.na(levels(x))] <- "unknown"
    return(x)
  })
  
  
# ts5g1[,101]
# ts5g1[21,]
#   hist(xx)
#################   remove by orders that are not marine     class(ts5good$Class)
 ts5good1<-ts5good[!ts5good$Class=="Amphibia",]   #removed 9 now 41, when no co - removed 15 now 71
ts5good1<-ts5good1[!ts5good1$Order=="Primates",] # removed 4 now 37, when no co - removed 6 now 65
ts5good1<-ts5good1[!ts5good1$Kingdom=="unknown",] # removed 8 now 29, when no co - removed 15 now 50
ts5good1<-ts5good1[!ts5good1$Phylum=="unknown",] # removed 6 now 23, when no co - removed 7 now 43
ts5good1<-ts5good1[!ts5good1$Family=="Geoemydidae",] # removed     when no co - removed 2 now 41
#ts5good1<-ts5good1[!ts5good1$Class=="unknown",] # removed 6 now 23, when no co - removed 22 now 19
  
#########################################################################################################
###########   updat summary table with new removals -- read_sum1  ################################################
qq<-ts5g1[,gc]  ## colnames(ts5g1)
read_sum1[6,]<-cbind(sum(qq),mean(colSums(qq)),sd(colSums(qq)),length(row.names(qq)), mean(colSums(qq>0,na.rm = TRUE)),sd(colSums(qq>0,na.rm = TRUE)))
qq<-ts5good[,gc]
read_sum1[7,]<-cbind(sum(qq),mean(colSums(qq)),sd(colSums(qq)),length(row.names(qq)), mean(colSums(qq>0,na.rm = TRUE)),sd(colSums(qq>0,na.rm = TRUE)))
qq<-ts5good1[,gc]  ##  colnames(ts5good1)
read_sum1[8,]<-cbind(sum(qq),mean(colSums(qq)),sd(colSums(qq)),length(row.names(qq)), mean(colSums(qq>0,na.rm = TRUE)),sd(colSums(qq>0,na.rm = TRUE)))
row.names(read_sum1)[6]<-"SV_1_sample" 
row.names(read_sum1)[7]<-"blank_based_removals" 
row.names(read_sum1)[8]<-"non_marine_order_removal" 
#    write.table(read_sum1, paste(dir,"Documents/KAUST/eDNA/R/export/dada/read_sum.csv",sep=""),row.names=T, sep=",")
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#################    reset depending on what you want filtered   #################################################
##
ts5good<-ts5good1
##
#########################################################################################################

############################################################################################
##make all with reduced taxa,   remove taxonomy for taxa accumulation 
ts5tran<-t(ts5good[,gc])
colnames(ts5tran)<-ts5good$Row.names
allg<-merge(ts5tran,sam, by="row.names",all.x=TRUE)  #  colnames(all2)  ###used for taxa acculmulati
all2<-allg  ### all2 was caluclated previously with all taxa/reads
############################################################################################
################################################################################################################
#############################################################################################
####              begin taxa accumulation
############################################################################################################
#####    2 plots        ####################################################
#dev.off()
par(mfrow=c(2,1),mar=c(1.5,1.5,1,1),oma=c(2,2,1,1))
####################################### 16 subsamples
all3<-all2[all2$subsamples=="1" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat, all3$dna_extraction),]
svv<-c(2:42)
dat<-all3[,svv]   ##  colnames(all3)     z<-all3[,135:151]
cc<-c("red","orange","green","blue")
ll<-c(2,1,2,1)
df <- lapply(c(1,17,33,49),function(i) specaccum(dat[seq(i,i+15),], method="random",permutations = 1000))
dfs<-df
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.12), ci.lty=0, lwd=3, lty=ll[1],
     xlim=c(1.5,16), ylim=c(0,30), ylab="Number of sequence variants", xlab="Number of subsamples")
for (i in 2:4) {
  plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0,lty=ll[i], 
       lwd=3)}
tit<-c("Coral Reef, extracellular DNA","Coral Reef, intra- and extracellular DNA","Seagrass meadow, extracellular DNA","Seagrass meadow, intra- and extracellular DNA")
legend(6,13, tit, col=cc, lty=ll, lwd=2,
       bty='n', inset=.1, y.intersp=1.5, merge=T)
text(1.2,30,labels="A")
####################################### though time plto 2
all3<-all2[all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat),]
dat<-all3[,svv]  ##   zz<-all3[,135:151]
cc<-c("red","green")
ll<-c(2)
df <- lapply(c(1,12),function(i) specaccum(dat[seq(i,i+9),], method="random",permutations = 1000))
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.1), ci.lty=0, lwd=3, lty=ll[1],
     xlim=c(1.5,10), ylim=c(0,30), ylab="Number of OTUs", xlab="Number of samples (taken every 2 weeks)")
i=2
plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0, lty=ll, 
     lwd=3)
tit<-c("Coral reef, extracellular DNA","Seagrass meadow, extracellular DNA")
legend(5,7, tit, col=cc, lty=ll, lwd=2,
       bty='n', inset=.1, y.intersp=1.5, merge=T)
text(1.3,30,labels="B")
ex<-1.1
mtext("Number of samples", SOUTH<-1,line=0.4, cex=ex, at=.5 ,outer=TRUE)
mtext("Number of sequence variants", WEST<-2,las=0,line=0.55, cex=ex, at=.5 ,outer=TRUE)
##################################################################################################################
##################################################################################################################
####################################### 1 plot of species sumation    ####################################
## spatial samples and blanks
all3<-all2[all2$subsamples=="1" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat, all3$dna_extraction),]
svv<-c(2:42)
dat<-all3[,svv]   ## 
cc<-c("red","dark red","green","dark green")
ll<-c(4,1,4,1)
###
df <- lapply(c(1,17,33,49),function(i) specaccum(dat[seq(i,i+15),], method="random",permutations = 1000))
dfs<-df
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.12), ci.lty=0, lwd=4, lty=1,
     xlim=c(1.5,16), ylim=c(0,30), ylab="Number of sequence variants", xlab="Number of subsamples")
for (i in 2:4) {
  plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0,lty=1, 
       lwd=4)}
tit<-c("Coral Reef, extracellular DNA","Coral Reef, intra- and extracellular DNA","Seagrass meadow, extracellular DNA","Seagrass meadow, intra- and extracellular DNA")
titn<-c("Spatial replicates","Temporal replicates") ##for legend 2
legend(6.5,6.5, tit, col=cc, lty=1, lwd=2,
       bty='n', inset=.1, y.intersp=1.5, merge=T)

## add temporal samples
all3<-all2[all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat),]
dat<-all3[,svv]  ##   zz<-all3[,135:151]
cc<-c("red","green")
ll<-c(1)
df <- lapply(c(1,12),function(i) specaccum(dat[seq(i,i+9),], method="random",permutations = 1000))
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.1), ci.lty=0, lwd=1.3, lty=ll[1],
      ylab="Number of OTUs", xlab="Number of samples (taken every 2 weeks)", add=T)
i=2
plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0, lty=ll, 
     lwd=1.3)

legend(6.5,10.5, titn, col="black", lty=c(1), lwd=c(4,1.1),
       bty='n', inset=.1, y.intersp=1.5, merge=T)
###########################################################################################################
###########################################################################################################
###########################################################################################################
######### keep only abundant SV's   ############################################################################
all3<-all2[all2$subsamples=="1" | all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat, all3$dna_extraction),]   ## names(all3)
svv<-c(2:42)
dat<-all3[,svv]   ## 
##################################################
##    get number individuals in 1,5,10 smaples
q<-rep(0,times=1000)
for (i in 1:1000) {
 #g<-sample(1:16, 10, replace=F)
 g<-sample(17:32, 5, replace=F)
m<-dat[g,]
sm<-sum(m)
q[i]<-sm
}
describe(q)
# smaples       5                10
#coral ex   6354 , 2701      12611  2822
# croal all 18794,10926      36126  10999
################################################
dat2<-dat[,order(colSums(dat),decreasing=T)]
#  hist(colSums(dat2))
t<-sum(dat2)  ## 109014
t4<-sum(dat2[,1:4])  ## 
t4/t    ## .62
t6<-sum(dat2[,1:6])  ## 73750
t6/t    ## 5 had .68  and 6 had 72
t10<-sum(dat2[,1:10])  ## 102156
t10/t    ## .88  11 had 90%
t20<-sum(dat2[,1:20])  ## 180098
t20/t    ## .97
df <- lapply(c(1,17,33,49),function(i) specaccum(dat[seq(i,i+15),], method="random",permutations = 1000))
df[[4]] $sites
###########################################################################################################
###########################################################################################################
###    species sum plots with most abudnant 10 and 20 SV's
#  dat10<-dat2[,1:10]  ##
#  dat20<-dat2[,1:20]  ##

par(mfrow=c(2,1),mar=c(1.5,1.5,1,1),oma=c(2,2,1,1))
####################################### top 10   ####################################
## spatial samples and blanks
all3<-all2[all2$subsamples=="1" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat, all3$dna_extraction),]
svv<-c(2:42)
dat<-all3[,svv]   ## 
dat2<-dat[,order(colSums(dat),decreasing=T)]
dat2<-dat2[,1:10]  ##
dat<-dat2
cc<-c("red","dark red","green","dark green")
ll<-c(4,1,4,1)
###
df <- lapply(c(1,17,33,49),function(i) specaccum(dat[seq(i,i+15),], method="random",permutations = 1000))
dfs<-df
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.12), ci.lty=0, lwd=4, lty=1,
     xlim=c(1.5,16), ylim=c(0,12), ylab="Number of sequence variants", xlab="Number of subsamples")
for (i in 2:4) {
  plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0,lty=1, 
       lwd=4)}
tit<-c("Coral Reef, extracellular DNA","Coral Reef, intra- and extracellular DNA","Seagrass meadow, extracellular DNA","Seagrass meadow, intra- and extracellular DNA")
titn<-c("Spatial replicates","Temporal replicates") ##for legend 2
legend(6.5,3, tit, col=cc, lty=1, lwd=2,
       bty='n', inset=.1, y.intersp=1.5, merge=T)
## add temporal samples
all3<-all2[all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat),]
dat<-all3[,svv]  ##   zz<-all3[,135:151]
dat2<-dat[,order(colSums(dat),decreasing=T)]
dat2<-dat2[,1:10]  ##
dat<-dat2
cc<-c("red","green")
ll<-c(1)
df <- lapply(c(1,12),function(i) specaccum(dat[seq(i,i+9),], method="random",permutations = 1000))
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.1), ci.lty=0, lwd=1.3, lty=ll[1],
     ylab="Number of OTUs", xlab="Number of samples (taken every 2 weeks)", add=T)
i=2
plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0, lty=ll, 
     lwd=1.3)
legend(6.5,4.5, titn, col="black", lty=c(1), lwd=c(4,1.1),
       bty='n', inset=.1, y.intersp=1.5, merge=T)
text(1.2,12,labels="A")
####################################### top 20   ####################################
## spatial samples and blanks
all3<-all2[all2$subsamples=="1" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat, all3$dna_extraction),]
svv<-c(2:42)
dat<-all3[,svv]   ## 
dat2<-dat[,order(colSums(dat),decreasing=T)]
dat2<-dat2[,1:20]  ##
dat<-dat2
cc<-c("red","dark red","green","dark green")
ll<-c(4,1,4,1)
###
df <- lapply(c(1,17,33,49),function(i) specaccum(dat[seq(i,i+15),], method="random",permutations = 1000))
dfs<-df
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.12), ci.lty=0, lwd=4, lty=1,
     xlim=c(1.5,16), ylim=c(0,25), ylab="Number of sequence variants", xlab="Number of subsamples")
for (i in 2:4) {
  plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0,lty=1, 
       lwd=4)}
## add temporal samples
all3<-all2[all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat),]
dat<-all3[,svv]  ##   zz<-all3[,135:151]
dat2<-dat[,order(colSums(dat),decreasing=T)]
dat2<-dat2[,1:20]  ##
dat<-dat2
cc<-c("red","green")
ll<-c(1)
df <- lapply(c(1,12),function(i) specaccum(dat[seq(i,i+9),], method="random",permutations = 1000))
plot(df[[1]], col=cc[1] , ci.type="polygon", ci.col=adjustcolor(cc[1],alpha.f=0.1), ci.lty=0, lwd=1.3, lty=ll[1],
     ylab="Number of OTUs", xlab="Number of samples (taken every 2 weeks)", add=T)
i=2
plot(df[[i]], add=T, col=cc[i] , ci.type="polygon", ci.col=adjustcolor(cc[i], alpha.f=0.12), ci.lty=0, lty=ll, 
     lwd=1.3)
text(1.2,25,labels="B")
ex<-1
mtext("Number of samples", SOUTH<-1,line=0.4, cex=ex, at=.5 ,outer=TRUE)
mtext("Number of sequence variants", WEST<-2,las=0,line=0.55, cex=ex, at=.5 ,outer=TRUE)
##################################################################################################################
#################################################################################################################
##################################################################################################################
######    figure % unique and overlap of SV's  ##############################################################################
all3<-all2[all2$subsamples=="1" | all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat, all3$dna_extraction),]   ## names(all3)    1:46 reef
svv<-c(2:42)
reef<-all3[all3$habitat=="coral",svv]
sg<-all3[all3$habitat=="seagrass",svv]
reef<-reef[,colSums(reef)>0] ## 33 SV
sg<-sg[,colSums(sg)>0] ## 33 SV
share<-sg[,-which(colnames(sg) %in% colnames(reef))]  ## 25 in common  8 unique
uniq<-reef[,-which(colnames(reef) %in% colnames(sg))] 
33/41
##################################################################################################################
#################################################################################################################
####################################################################################################################
##  time series ???????
all3<-all2[all2$subsamples=="5" & all2$sample_type=="sediment",]#only single subsamples
all3<-all3[order(all3$habitat),]
dat<-all3[,svv]  ##   zz<-all3[,135:151]
cc<-c("red","green")
all3$richnes<-rowSums(dat>0,na.rm = TRUE)
mess<-all3[,c(1,46,60)]  ## names(all3)

##################################################################################################################
#################################################################################################################
####################################################################################################################
###############################################################################################################
#######################################################################################################
########  ##            ordination    ###################################################################
all3<-all2[!all2$sample_type=="trap",]#remove trap
all3<-all3[all3$subsamples=="1" | all3$subsamples=="5",]# 1 5 and blank

all3<-all3[order(all3$habitat, all3$dna_extraction),]
all3$habitat<-factor(all3$habitat,levels=c("coral","seagrass","blank"))
all3$samp_cat<-factor(all3$samp_cat,levels=c("spatial","temporal")) ##  levels(all3$samp_cat)
all3$hab_dna<-factor(all3$hab_dna,levels=c("coral extracellular","coral total","seagrass extracellular","seagrass total")) #  #  levels(all3$hab_dna)
all3$for_mer<-factor(all3$for_mer,levels=c("coral temporal extracellular","coral spatial extracellular","coral spatial total",            
                                      "seagrass temporal extracellular","seagrass spatial extracellular","seagrass spatial total")) #  
## unique(all3$for_mer)
tit<-c("Coral reef, extracellular DNA","Coral reef, intra- and extracellular DNA","Seagrass meadow, extracellular DNA","Seagrass meadow, intra- and extracellular DNA", "Blanks and control")
tit2<-c("Spatial replicates","Temporal replicates","Positive control","Extraction blank","PCR blank") #
titn<-c("Spatial replicates","Temporal replicates")
svv<-c(2:42)
dat<-all3[,svv]   ##  colnames(all3)     z<-all3[,135:150]    colnames(z)
dat<-sqrt(dat)    
#dat<-log(dat+1)  
#dat<-dat
all3$habdna<-paste(all3$habitat, all3$dna_extraction)
all3$habdna<-as.factor(all3$habdna)
levels(all3$habdna)
dat.env<-all3[,c(1,43:59)]
data(dat.env)
cc<-c("red","dark red","green","dark green") ##c ex, c all , g ex, g all, blanks 
cc2<-c("red","red","dark red","green","green","dark green") ##c ex, c all , g ex, g all, blanks 
ll<-c(1,4,4,1,4,4)
pp<-c(1,2,3,4,8)   #21,23
ss<-c(2,1,1)
ord<-metaMDS(dat,trace=F,autotransform = F)
#ord<-rda(dat,trace=F)
plot(ord, disp="sites", type='n', xlim=c(-2,2), ylim=c(-1,1))
points(ord, col=cc[dat.env$hab_dna], cex=1, pch = pp[dat.env$samp_cat])
ordiellipse(ord, dat.env$for_mer, kind = "se", conf = 0.95, col=cc2,lwd=ll,lty=1 )    ## levels(dat.env$for_mer)
# add legend
xx<--2.13
legend(xx,-0.95, col=cc, legend=tit, pch = 15, cex = 1,bty='n', y.intersp=1.5)
#legend(xx,0.5, col="black", legend=tit2, pch = pp, cex = 1,bty='n', y.intersp=1.5)
legend(xx,-0.6, titn, col="black", pch = pp[c(1,2.5)], lty=c(1), lwd=c(3.2,1.1), bty='n', inset=.1, y.intersp=1.5, merge=T)
#############################################################################################
#############################################################################################
#######################################################################################################
#############################################################################################
#############################################################################################
#######################################################################################################
########              heat maps        ##################################################
###################################################################################################
######     combine all SV by Family and sum number of reads
#######
ts5<-ts5good1    ###   change depending on filtering wanted       colnames(ts5)
######
fam.5sum<-aggregate(ts5[,-c(1:8)],ts5[,2:6], FUN=sum)    ## sum reads for al Seq Var
fam.5len<-aggregate(ts5[,-c(1:8)],ts5[,2:6], function(x) length(which(x>0)))  
###   sort
fam.5len1<-fam.5len[order(fam.5len[,1],fam.5len[,2],fam.5len[,1],fam.5len[,3],fam.5len[,4],fam.5len[,5]),] 
fam.5sum1<-fam.5sum[order(fam.5sum[,1],fam.5sum[,2],fam.5sum[,1],fam.5sum[,3],fam.5sum[,4],fam.5sum[,5]),]
#fam.8len1<-fam.8len[order(fam.8len[,1],fam.8len[,2],fam.8len[,1],fam.8len[,3],fam.8len[,4],fam.8len[,5]),]
##  make separte pylogeny table to join pre heat map
pyl<-fam.5sum1[,1:5]
pylo<-pyl
######################################################
###make locations for heatmap
pyl$count<-c(1:length(pyl[,1]))
pyl$ord_count<-sequence(rle(as.character(pyl$Order))$lengths)
pyl$cla_count<-sequence(rle(as.character(pyl$Class))$lengths)
pyl$ph_count<-sequence(rle(as.character(pyl$Phylum))$lengths)

# for order
pyl_ord<-pyl[pyl$ord_count==1,]
pyl_ord$Cont_count_len<-1
pyl_ord$Cont_count_len[c(1:length(pyl_ord$Cont_count_len)-1)]<-tail(pyl_ord$count, -1) - head(pyl_ord$count, -1)
pyl_ord$cc_end<-(pyl_ord$Cont_count_len-1)+pyl_ord$count
mm<-cbind(pyl_ord$count,pyl_ord$cc_end)
pyl_ord$yloc<-rowMeans(mm)
pyl_ord$yloc2<-(length(pyl[,1])+1)-pyl_ord$yloc
pyl_ord$top<-(length(pyl[,1])+1)-pyl_ord$count
pyl_ord$bot<-(length(pyl[,1])+1)-pyl_ord$cc_end
pyl_ord$top<-pyl_ord$top+0.1
pyl_ord$bot<-pyl_ord$bot-0.1
# for class
pyl_cla<-pyl[pyl$cla_count==1,]
pyl_cla$Cont_count_len<-1
pyl_cla$Cont_count_len[c(1:length(pyl_cla$Cont_count_len)-1)]<-tail(pyl_cla$count, -1) - head(pyl_cla$count, -1)
pyl_cla$cc_end<-(pyl_cla$Cont_count_len-1)+pyl_cla$count
mm<-cbind(pyl_cla$count,pyl_cla$cc_end)
pyl_cla$yloc<-rowMeans(mm)
pyl_cla$yloc2<-(length(pyl[,1])+1)-pyl_cla$yloc
pyl_cla$top<-(length(pyl[,1])+1)-pyl_cla$count
pyl_cla$bot<-(length(pyl[,1])+1)-pyl_cla$cc_end
pyl_cla$top<-pyl_cla$top+0.1
pyl_cla$bot<-pyl_cla$bot-0.1
#####################################################
###############################################
######     combine all SV by treatment   get mean per sample and SD and length(N)
tsam<-t(sam)  # rownames()     x["1",]
t.fam<-data.frame(t(fam.5sum1[,-c(1:5)]))
t.s<-merge(sam[,c(10:13,15:17)], t.fam, by="row.names") ##  colnames(t.s)   class(t.s$X1)   class(t.s$for_mer)
t.s<-data.frame(t.s,stringsAsFactors=T)
#fam.5len1 for nubmer of SV
tl.fam<-data.frame(t(fam.5len1[,-c(1:5)]))
tl.s<-merge(sam[,c(10:13,15:17)], tl.fam, by="row.names") ##  
tl.s<-data.frame(tl.s,stringsAsFactors=T)

el<-c(2:8)
sl<-c(9:length(colnames(t.s)))
t.s.mean<-aggregate(t.s[,sl],t.s[,el], FUN=mean)    ## sum reads for al Seq Var    colnames(t.s)
t.s.sd<-aggregate(t.s[,sl],t.s[,el], FUN=sd)    ## sd for al Seq Var
#   check length need use   ---   FUN=sum  of   fam.5len1!!!!!!!!
t.s.length<-aggregate(t.s[,sl],t.s[,el], FUN=length)    ##   now number of reps
t.s.sv<-aggregate(tl.s[,sl],tl.s[,el], FUN=mean)      ## number of sv al Seq Var  
asl<-c(8:length(colnames(t.s.sd)))   ##  colnames(t.s.sd)
t.s.se<-t.s.sd[,asl]/sqrt(t.s.length[,asl])
##
## keep only target treatments   ###############################
## mean   ################
x<-t.s.mean
x1<-x[!x$sample_type=="trap",]#remove trap
x1<-x1[x1$subsamples=="1" | x1$subsamples=="5" | x1$subsamples=="blank",]# 1 5 and blank

x1$for_mer<-factor(x1$for_mer, levels=c("coral spatial extracellular", "coral spatial total","coral temporal extracellular",
                                        "seagrass spatial extracellular","seagrass spatial total","seagrass temporal extracellular",
                                        "blank positive control blank","blank extraction blank blank","blank pcr blank blank"))
x1<-x1[order(x1$for_mer),] #### re order to factor created above
ttx<-t(x1[,8:length(names(x1))])  ### transpose dv rows with column treatments    colnames(x1)
ttx2<-cbind(pylo,ttx) ### add phylogeny
ttmean<-ttx2
## creat matrix
zz<-ttmean
zzz<-c(6:14) #row.names(zz) <- c(1:16
colnames(zz)[zzz]<- levels(x1$for_mer)
rn<-as.character(pyl$Family)
rn[16:18]<-""
cn<-levels(x1$for_mer)
cc<-data.matrix(zz[,zzz])   #   dd<-describe(cc)
cc<-log(cc+1) # describe(cc)
tcc<-t(cc[nrow(cc):1,])  # to make image plot correct need to turn matrix
mmcc<-tcc
## sd  ################
x<-t.s.sd
x1<-x[!x$sample_type=="trap",]#remove trap
x1<-x1[x1$subsamples=="1" | x1$subsamples=="5" | x1$subsamples=="blank",]# 1 5 and blank
x1$for_mer<-factor(x1$for_mer, levels=c("coral spatial extracellular", "coral spatial total","coral temporal extracellular",
                                        "seagrass spatial extracellular","seagrass spatial total","seagrass temporal extracellular",
                                        "blank positive control blank","blank extraction blank blank","blank pcr blank blank"))
x1<-x1[order(x1$for_mer),] #### re order to factor created above
ttx<-t(x1[,8:length(names(x1))])  ### transpose dv rows with column treatments    colnames(x1)
ttx2<-cbind(pylo,ttx) ### add phylogeny
ttsd<-ttx2
## creat matrix
zz<-ttsd
zzz<-c(6:14) #row.names(zz) <- c(1:16
colnames(zz)[zzz]<- levels(x1$for_mer)
rn<-as.character(pyl$Family)
rn[16:18]<-""
cn<-levels(x1$for_mer)
cc<-data.matrix(zz[,zzz])   #   summary(cc)
cc<-log(cc+1) # describe(cc)
tcc<-t(cc[nrow(cc):1,])  # to make image plot correct need to turn matrix
smcc<-tcc
## len #########
x<-t.s.sv
x1<-x[!x$sample_type=="trap",]#remove trap
x1<-x1[x1$subsamples=="1" | x1$subsamples=="5" | x1$subsamples=="blank",]# 1 5 and blank
x1$for_mer<-factor(x1$for_mer, levels=c("coral spatial extracellular", "coral spatial total","coral temporal extracellular",
                                        "seagrass spatial extracellular","seagrass spatial total","seagrass temporal extracellular",
                                        "blank positive control blank","blank extraction blank blank","blank pcr blank blank"))
x1<-x1[order(x1$for_mer),] #### re order to factor created above
ttx<-t(x1[,8:length(names(x1))])  ### transpose dv rows with column treatments    colnames(x1)
ttx2<-cbind(pylo,ttx) ### add phylogeny
ttsv<-ttx2
## creat matrix
zz<-ttsv
zzz<-c(6:14) #row.names(zz) <- c(1:16
colnames(zz)[zzz]<- levels(x1$for_mer)
rn<-as.character(pyl$Family)
rn[9:10]<-""   # remove last few unknowns
cn<-levels(x1$for_mer)
cc<-data.matrix(zz[,zzz])   #   dd<-describe(cc)
tcc<-t(cc[nrow(cc):1,])  # to make image plot correct need to turn matrix
svmcc<-tcc
################################################################################################################
##############################################################################################################
##############################################################################################################
####    heat map
## creat list for points
#  get lsit of means for plot############
mcc<-melt(mmcc)
mcc$value[mcc$value == 0] <- NA  ###better for plotting ######
mcc$x=rep(1:nrow(mmcc), len=nrow(mcc))
mcc$y=rep(1:ncol(mmcc), each=nrow(mmcc))
zlim=c(min(mcc$value,na.rm=T),max(mcc$value,na.rm=T))
nlevels<-50
levels <- seq(zlim[1],zlim[2],length.out = nlevels)
my_palette_p <- colorRampPalette(c("light blue","green", "yellow","red"))(nlevels)  
mcc$col <- my_palette_p[cut(mcc$value,nlevels)]  
my_palette <- colorRampPalette(c("black","white"))(nlevels)
#  get number of sv   ############
lcc<-melt(svmcc)
lcc$value[lcc$value == 0] <- NA  ###better for plotting
lcc$x=rep(1:nrow(svmcc),len=nrow(lcc))
lcc$y=rep(1:ncol(svmcc),each=nrow(svmcc))
zlim=c(min(lcc$value,na.rm=T),max(lcc$value,na.rm=T))
####################################################################################################
##################################################
## image.plot        use image if no legend (or if have mulipte plots), but use image.plot is want legend
tit<-c("Coral reef, extracellular DNA","Coral reef, intra- and extracellular DNA","Seagrass meadow, extracellular DNA","Seagrass meadow, intra- and extracellular DNA", "Blanks and control")
tit2<-c("Spatial replicates","Temporal replicates","Positive control","Extraction blank","PCR blank") #
ticks<-c(0,10,100,1000,10000,100000)
ticks2<-c(0,10,100,1000,10000)
yy = c(1:length(rownames(cc)))# extract latitudes
xx = c(1:length(colnames(cc))) # extract longitudes
ccc<-c("red","dark red","green","dark green","black") ##
tit3<-c("Spatial replicates","Spatial replicates","Temporal replicates","Spatial replicates","Spatial replicates","Temporal replicates","Positive control","Extraction blank","PCR blank") #
######### x y matrix of x   PLOT   ####################################################################
####################################################################
####################################################################
  par(mfrow=c(1,1),mar=c(10,20,1,10),oma=c(0,0,1,2)) 
  image(xx,yy,smcc, axes=F ,col=rev(my_palette), xlab ="", ylab ="")
  #points(0,0) # hack so axis will display properly
  axis(side=2, at= 1:length(rn), labels=rev(rn), las=1, cex.axis=0.8)
  box()
 ## x labels    ########
  tt<-c("red","dark red","red","green","dark green","green","black","black","black")
  axis(side=1, at=1:9, labels=F, las=2, cex.axis = 0.8)
  mtext(tit3,side=1,line=1,at=c(1:9),las=2, adj=1,cex=.8,col=tt)

########     plot circles  ##########################  
points(mcc$x,mcc$y,col=mcc$col,bg=mcc$col, pch=21, cex=2.5)  

########    add text to circles  ##########################  
 text(lcc$x,lcc$y,labels=round(lcc$value, digits=1), cex=.55) 
  
#######  add order and class y axis labels ########################
## add order 
  mm<-pyl_ord[-nrow(pyl_ord),]## remove last row
  yt<-(mm$yloc2)
  lt<-(mm$Order)
mtext(lt,side=2,line=8.5,at=yt,las=1, adj=.5,cex=.8)
## add class
mm<-pyl_cla
yt<-c(pyl_cla$yloc2)
lt<-as.character(pyl_cla$Class)
lt[2]<-"unknown Mammalia"
lt[3]<-"unknown Cordata"
#lt<-c(lt,"unknown Metazoa","unknown")
mtext(lt,side=2,line=13,at=yt,las=1, adj=.5,cex=.8)
##  label phylogeny  #######
lt<-c("Class","Order","Family")
xt<-c(13.1,8.5,3)
mtext(lt,side=2,line=xt,at=10.5,las=1, adj=.5,cex=0.8,font=2)
###add lines to empasise groupings
## order
mm<-pyl_ord
mm2<-mm[mm$Cont_count_len>1,]
for (i in 1:length(mm2$top)){
axis(2,at=c(mm2$bot[i],mm2$top[i]),col="black",line=6,tick=T,labels=rep("",2),lwd=1,lwd.ticks=0)  }
## class
mm<-pyl_cla
mm2<-mm[mm$Cont_count_len>1,]
for (i in 1:length(mm2$top)){
  axis(2,at=c(mm2$bot[i],mm2$top[i]),col="black",line=11.3,tick=T,labels=rep("",2),lwd=1,lwd.ticks=0)  }

##x axis label legend 
par(xpd=TRUE)
legend(-12.5,-0.5, col=ccc, legend=tit, pch = 15, cex = 0.85, bty='n', y.intersp=1.5)  
mtext("X axis labels",side=2,line=7.7,at=-.4,las=2, cex=.85, font=2)
###
text(11.5,10,labels="Num of reads",srt=45,cex=0.8,font=1, pos=1)
text(12.4,9.7,labels="St. dev.",srt=45,cex=0.8,font=1, pos=1)
###square legend
image.plot(zlim =c(0,log(max(ticks2))),  nlevel = nlevels,legend.only = TRUE, 
            horizontal=FALSE, col=rev(my_palette),  axis.args=list( at=log(ticks2+1),labels=ticks2, cex.axis=.8),
                smallplot= c(.81,.83,.35,.85), graphics.reset=T) 
## point legend
image.plot(zlim =c(0,log(max(ticks2))), nlevel = nlevels ,legend.only = TRUE,
           horizontal=FALSE, col=my_palette_p ,axis.args=list(at=log(ticks2+1), labels=F, cex.axis=.8),
           legend.args = list(cex = .8, side = 4, line=0, title="Num of reads"), smallplot= c(.775,.795,.35,.85), 
           graphics.reset=T)


################################################################################################################
  dev.off()


  
  
  
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
### other heatmap examples
################################################################################################################
######     corplot
par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,0))
corrplot(c, is.corr = FALSE, method = "circle", mar=c(1,0,0,0)) ## bg=background color


#####################################################
library(ggplot2)
##normal
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
## e cirles
ggplot(mat, aes(Name, variable)) +
  geom_point(aes(size = value, colour=value))
###squares and circles
ggplot(dataf, aes(y = factor(rowv),
                  x = factor(columnv))) +        ## global aes
  geom_tile(aes(fill = rectheat)) +         ## to get the rect filled
  geom_point(aes(colour = circlefill, 
                 size =circlesize))  +    ## geom_point for circle illusion
  scale_color_gradient(low = "yellow",  
                       high = "red")+       ## color of the corresponding aes
  scale_size(range = c(1, 20))+             ## to tune the size of circles
  theme_bw()



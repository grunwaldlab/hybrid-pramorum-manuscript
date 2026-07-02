library(stringr)
library(ggpubr)
library(IRanges)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggpmisc)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(viridis)
library(hrbrthemes)
library(svglite)
library(tidyr)

blocksize <- 20000
sliding_increment <- 20000

allbed <- read.delim(file="shared_phased_regions.bed",
                     sep='\t', header=F)

chrsizes <- read.delim(file="pr102v4_chrom_sizes.tsv", sep='\t', header=F)

#recoding the big vcf
hybvcf <- read.delim(file='freebayes_wh.vcf', sep='\t', header = T, row.names = NULL, skip=88)

#colnames(allbed) <- gsub("EUxNA", "CR1", colnames(allbed))
hybvcf$QUAL <- NULL
hybvcf$FILTER <- NULL
hybvcf$FORMAT <- NULL
hybvcf$INFO <- NULL

colnames(hybvcf) <- gsub("NAxEU", "EUxNA", colnames(hybvcf))
colnames(hybvcf) <- gsub("X.CHROM", "CHROM", colnames(hybvcf))
hybvcf <- hybvcf[hybvcf$REF%in%c("A", "C", "G", "T"),]
hybvcf <- hybvcf[hybvcf$ALT%in%c("A", "C", "G", "T"),]

names_list <- colnames(hybvcf)[6:ncol(hybvcf)]
hybvcf1 <-hybvcf

#getting rid of phase info. That info has been stored elsewhere - the bed multiiner file
#the bed multiiner file is used to select haplotigs
for(i in 1:ncol(hybvcf)){
  hybvcf[,i] <- str_split_fixed(hybvcf[,i], ":", 2)[,1]
}


#Separate VCF by haplotype and remove some characters
hybvcf[] <- lapply(hybvcf, function(x) gsub("/", "|", x))

hybvcf2 <- hybvcf
colnames(hybvcf) <- paste(colnames(hybvcf), "_1", sep="")
colnames(hybvcf2) <- paste(colnames(hybvcf2), "_2", sep="")

hybvcf[] <- lapply(hybvcf, function(x) gsub("\\Q|\\E0", "", x)) 
hybvcf[] <- lapply(hybvcf, function(x) gsub("\\Q|\\E1", "", x))
hybvcf[] <- lapply(hybvcf, function(x) gsub("\\Q|\\E.", "", x))

hybvcf2[] <- lapply(hybvcf2, function(x) gsub("0\\Q|\\E", "", x))
hybvcf2[] <- lapply(hybvcf2, function(x) gsub("1\\Q|\\E", "", x))
hybvcf2[] <- lapply(hybvcf2, function(x) gsub(".\\Q|\\E", "", x))

hybvcf <- cbind(hybvcf, hybvcf2[,c(5:ncol(hybvcf2))])

hybvcf[] <- lapply(hybvcf, function(x) gsub("\\Q.\\E", NA, x)) 

combinations <- combn(names_list, 2)


#subsetting 

combinations_paste=paste0(combinations[1,],combinations[2,])

for(i in 1:(ncol(hybvcf)-3)){
  k<-i+3
  hybvcf[,k] <- as.numeric(hybvcf[,k])
}
#


hybvcf$ALT_1 <- NULL

colnames(hybvcf) <- gsub("EUxNA", "CR1", colnames(hybvcf))

hybvcf$REF_1 <- NULL
hybvcf$ALT_1 <- NULL
hybvcf$ID_1 <- NULL

hybvcf$POS_1 <- as.numeric(hybvcf$POS_1)
hybvcf$ALT_2 <- NULL


hapnames <- colnames(hybvcf)[3:ncol(hybvcf)]



#subset to regions present in both lineages
allbed <- allbed[, 1:3]
allbed[,2]
bed_gr <- GRanges(seqnames = allbed[,1],
              ranges = IRanges(start = allbed[,2], end = allbed[,3]))
bed_gr <- reduce(bed_gr)
allbed1 <- as.data.frame(bed_gr)
allbed1 <- allbed1[!allbed1$width<blocksize,]
allbed1$windows<-floor(x = allbed1$width/blocksize) 

comb_df <- combn(c(paste0(names_list, '_1'), paste0(names_list, '_2')), 2)

hyb_slide_df <- data.frame(hap1 = character(), hap2 = character(), chromosome = character(), start = numeric(),
                           end=numeric(), identity=numeric(),
                           notphased_start = numeric(), notphased_end=numeric(), pairname=character(), block_snps=numeric())

#three loops:
#first loop: get every phaseable region from rows in bed file
for(j in 1:nrow(allbed1)){
  region_j <- allbed1[j,]
  
  #second loop: 
  #break phaseable region into blocks
  
  
  for(k in 1:region_j$windows[1]){
    vcf_k <- hybvcf[hybvcf$CHROM_1 == region_j$seqnames[1]&
                    hybvcf$POS_1>=(region_j$start[1]+((k-1)*blocksize))&
                    hybvcf$POS_1<(region_j$start[1]+((k)*blocksize)),] 
    nsnps_k <- nrow(vcf_k)
    hybvcf$CHROM_1
    #third loop
    #compute pairwise SNP differences for every combination
    for(i in 1:ncol(comb_df)){
      
      hap1 <- comb_df[1,i]
      hap2 <- comb_df[2,i]
      
      chr_i <- as.character(region_j$seqnames[1])
      
      vcf_i <- vcf_k[,c(hap1,hap2)]
      
      vcf_i
      dif_i <- length(which(as.numeric(vcf_i[,1])+as.numeric(vcf_i[,2])==1))
      x <- c(hap1, hap2, chr_i, 
             (region_j$start[1]+((k-1)*blocksize)),
             (region_j$start[1]+((k)*blocksize)),dif_i,NA,NA,paste0(hap1,hap2), nsnps_k)
      
      hyb_slide_df[nrow(hyb_slide_df)+1,] <- x
     
    }
    
  }
  print(j/nrow(allbed1))
}

hyb_slide_df$identity  <- as.numeric(hyb_slide_df$identity)
write.table(hyb_slide_df, file = 'hyb_slide_df', sep='\t', quote=F)


###########################
hyb_slide_df <- read.delim(file = 'hyb_slide_df', sep='\t')
hyb_slide_df$start <- as.numeric(hyb_slide_df$start )

ggplot(hyb_slide_df, aes(x = identity))+
  geom_histogram(bins=351)+
  scale_x_continuous(limits=c(-1,350))

hyb_slide_df$block_snps <- as.numeric(hyb_slide_df$block_snps)
ggplot(hyb_slide_df, aes(x = block_snps))+
  geom_histogram()

hyb_slide_df$start <- as.numeric(hyb_slide_df$start)
hyb_slide_df$end <- as.numeric(hyb_slide_df$end)


#I need to make a df that shows xmin and xmax positions only when the target haplotype has an identity under the threshhold
#all others need to have na values


#first step: order haplotypes
#type1 X8959, X9389, X9202, X9203,
#type2 X9170, X9390

hyb_list <- c('X8959', 'X9389', 'X9202', 'X9203', 'X9170', 'X9390')

#To make plot look better, I added some dummy rows
hyb_order <- c('space0','space00','space000','space000','space0000','space00000','space000000',
  'X9390_1','X9390_2','space1','space11','X9389_1','X9389_2','space2','space22','X9170_1','X9170_2',
               'space3','space33', 'X8959_1','X8959_2',
               'space4','space44',
               'X9203_1','X9203_2','space5','space55','X9202_1','X9202_2')

parent_haps <- c('X9381_1','X9381_2','X9391_1','X9391_2')

new_hyb_df <- hyb_slide_df[0,]
new_hyb_df$parent_hap <- NA

hyb_col <- ncol(hyb_slide_df)

for(i in 1:length(parent_haps)){
  #only take parental haplotypes
  df1 <- hyb_slide_df[hyb_slide_df$hap1 %in% parent_haps[i] | hyb_slide_df$hap2 %in% parent_haps[i] , ]
  
  df1$parent_hap <- parent_haps[i]
  
  new_hyb_df <- rbind(new_hyb_df,df1)
  
}

new_hyb_df$index <- NA
new_hyb_df$off_hap <- NA


for(i in 1:length(hyb_order)){
  #
 
  new_hyb_df[new_hyb_df$hap1 %in% hyb_order[i] | new_hyb_df$hap2 %in% hyb_order[i],'off_hap'] <- hyb_order[i]
  new_hyb_df[new_hyb_df$hap1 %in% hyb_order[i] | new_hyb_df$hap2 %in% hyb_order[i],'index'] <- i
}

new_hyb_df <- new_hyb_df[!is.na(new_hyb_df$index),]
new_hyb_df$pairname <- paste0(new_hyb_df$parent_hap, new_hyb_df$off_hap)

new_hyb_df$loca <- paste(new_hyb_df$chromosome, new_hyb_df$start, sep=":")

pos_n <- length(unique(new_hyb_df$loca))

fourcomp_df <- data.frame(chr=character(), start=numeric(), end=numeric(), 
                             off_hap = character(), off_index=numeric(),
                             EU1_H1 =numeric(), EU1_H2 =numeric(), NA2_H1=numeric(), NA2_H2=numeric())


for(i in 1:length(unique(new_hyb_df$loca))){
  pos_i <- unique(new_hyb_df$loca)[i]
  df_i <- new_hyb_df[new_hyb_df$loca==unique(new_hyb_df$loca)[i],]
  
  for(k in(1:length(hyb_order))){
    hyb_order
    
    df_k <- df_i[df_i$off_hap == hyb_order[k],]
    new_df_k <- data.frame(chr=df_k$chromosome[1], start=df_k$start[1], end=df_k$end[1],
                   off_hap = df_k$off_hap[1], off_index=df_k$index[1],
                   EU1_H1 =NA, EU1_H2 =NA, NA2_H1=NA, NA2_H2=NA)
    new_df_k$EU1_H1 <- df_k$identity[1]
    new_df_k$EU1_H2 <- df_k$identity[2]
    new_df_k$NA2_H1 <- df_k$identity[3]
    new_df_k$NA2_H2 <- df_k$identity[4]
    
    fourcomp_df[nrow(fourcomp_df)+1,] <- new_df_k
  }
  print(i/pos_n)
}

unique(new_hyb_df$loca)

cutoff=5

fourcomp_df_1 <- fourcomp_df

fourcomp_df_1$EU1_H1[fourcomp_df_1$EU1_H1>cutoff] <- 11
fourcomp_df_1$EU1_H1[fourcomp_df_1$EU1_H1<=cutoff] <- 0

fourcomp_df_1$EU1_H2[fourcomp_df_1$EU1_H2>cutoff] <- 11
fourcomp_df_1$EU1_H2[fourcomp_df_1$EU1_H2<=cutoff] <- 0

fourcomp_df_1$NA2_H1[fourcomp_df_1$NA2_H1>cutoff] <- 11
fourcomp_df_1$NA2_H1[fourcomp_df_1$NA2_H1<=cutoff] <- 0

fourcomp_df_1$NA2_H2[fourcomp_df_1$NA2_H2>cutoff] <- 11
fourcomp_df_1$NA2_H2[fourcomp_df_1$NA2_H2<=cutoff] <- 0

fourcomp_df_1$EU1_H1 <- gsub(11, 'Y', fourcomp_df_1$EU1_H1)
fourcomp_df_1$EU1_H1 <- gsub(0, 'N', fourcomp_df_1$EU1_H1)

fourcomp_df_1$EU1_H2 <- gsub(11, 'Y', fourcomp_df_1$EU1_H2)
fourcomp_df_1$EU1_H2 <- gsub(0, 'N', fourcomp_df_1$EU1_H2)

fourcomp_df_1$NA2_H1 <- gsub(11, 'Y', fourcomp_df_1$NA2_H1)
fourcomp_df_1$NA2_H1 <- gsub(0, 'N', fourcomp_df_1$NA2_H1)

fourcomp_df_1$NA2_H2 <- gsub(11, 'Y', fourcomp_df_1$NA2_H2)
fourcomp_df_1$NA2_H2 <- gsub(0, 'N', fourcomp_df_1$NA2_H2)


options(scipen=999)

fourcomp_df_1 <- fourcomp_df_1[!is.na(fourcomp_df_1),]

library(ggh4x)


############
#dynamic scaling of x axis so that the margins are less skewed
chr_i<-'pchr1'
chrsize_i <- chrsizes[chrsizes$V1==chr_i, 'V2']
fourcomp_df_subset <- fourcomp_df_1[fourcomp_df_1$chr==chr_i,]

xr_kb <- chrsize_i

in_per_mb <- 0.000001
overhead_in <- 1.70   # axis labels, margins, etc.

plot_width <- overhead_in + xr_kb * in_per_mb

#overlap boxes

b <- sort(unique(fourcomp_df_subset$start))

start <- b
end   <- b + sliding_increment

# start a new group when the next start is > previous end
grp <- cumsum(c(TRUE, start[-1] > end[-length(end)]))

xmin <- tapply(start, grp, min)
xmax <- tapply(end,   grp, max)

boxes <- data.frame(
  xmin = as.numeric(xmin),
  xmax = as.numeric(xmax)
)


fourcomp_df_10 <- fourcomp_df_1[fourcomp_df_1$chr==chr_i,]

figure1 <- ggplot(fourcomp_df_10, aes())+
  geom_rect(aes(xmin = start+5000, xmax = start+20000,
                ymin = off_index, ymax = off_index+1,
                alpha=NA2_H1), fill="blue",
            color='black', linewidth =0.00025)+
  geom_rect(aes(xmin = start+5000, xmax = start+20000,
                ymin = off_index, ymax = off_index+1,
                alpha=NA2_H2), fill="green",
            color='black', linewidth =0.00025)+
  geom_rect(aes(xmin = start+5000, xmax = start+20000,
                ymin = off_index, ymax = off_index+1,
                alpha=EU1_H1), fill="orange",
            color='black', linewidth =0.00025)+
  geom_rect(aes(xmin = start+5000, xmax = start+20000,
                ymin = off_index, ymax = off_index+1,
                alpha=EU1_H2), fill="red",
            color='black', linewidth =0.00025)+
  scale_alpha_manual(values = c('Y'=0, 'N'=0.75))+
  theme_classic()+
  scale_y_continuous(
    limits=c(8,40),
    breaks = c(9, 13, 17, 21, 25, 29, 38, 34),
    labels = c("9390","9389",'9170','8959','9203','9202',"NA2","EU1")
  )+
  geom_rect(
    data = boxes,
    aes(xmin = 0, xmax = chrsize_i, ymin = 33, ymax = 34),
    inherit.aes = FALSE,
    fill='#ffbc57',color='black',  alpha=1 )+
  geom_rect(
    data = boxes,
    aes(xmin = 0, xmax = chrsize_i, ymin = 34, ymax = 35),
    inherit.aes = FALSE,
    fill='#ff5a56',color='black',  alpha=1 )+  
  geom_rect(
    data = boxes,
    aes(xmin = 0, xmax = chrsize_i, ymin = 37, ymax = 38),
    inherit.aes = FALSE,
    fill='#66f956',color='black',  alpha=1 )+
  geom_rect(
    data = boxes,
    aes(xmin = 0, xmax = chrsize_i, ymin = 38, ymax = 39),
    inherit.aes = FALSE,  alpha=1 ,
    fill='#535fff',color='black')+  
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = 32.5, ymax = 39.5),
    inherit.aes = FALSE, fill='grey',
    alpha=0, color='black', linetype = 'dotted')+
  xlab('mb')+
    scale_x_continuous(breaks=seq(from=0, to=10000000, by=1000000), labels = function(x) x / 1000000)+
  ggtitle(chr_i)+
  theme(text=element_text(size=20),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.title.y = element_blank(),
       # axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
       axis.line.y = element_blank(),
       legend.position = 'none'
        )
 # add this section if looking at chromosome breakpoints from particular regions
 # scale_x_continuous(limits=c(8200000,8500000))+ 
 # geom_text(
  #  data = fourcomp_df_1[fourcomp_df_1$chr==chr_i,],
  #  aes(
  #    x = start + 12500,  # midpoint: (5000 + 20000)/2
  #    y = off_index + 0.5,
  #    label = paste0(start, "-", end)
  #  ),
  #  size = 3, angle = 45
  #)

###########
# figure 2: histogram summarizing shared haplotypes per block

detach('gdtools')

#in paper: 
#sibs = X9202, X9170
#clones = X9202, X9203
#parent/offspring = X9391, X9390
#parent/parent = X9381, X9391

identity_cutoff = 11

testhapA <- 'X9202'
testhapB <- 'X9170'

testhapA1 <- paste0(testhapA, '_1')
testhapA2 <- paste0(testhapA, '_2')
testhapB1 <- paste0(testhapB, '_1')
testhapB2 <- paste0(testhapB, '_2')


df_comp <- data.frame(pos = unique(hyb_slide_df$pos))
df_comp$count <- NA

#loop counting each comp, with 
for(i in 1:nrow(df_comp)){
  comp1 = as.numeric(hyb_slide_df[hyb_slide_df$hap1%in%c(testhapA1, testhapB1) & 
                 hyb_slide_df$hap2%in%c(testhapA1, testhapB1)&
                    hyb_slide_df$pos==df_comp$pos[i],'identity']<identity_cutoff)
  comp2 = as.numeric(hyb_slide_df[hyb_slide_df$hap1%in%c(testhapA2, testhapB1) & 
                                    hyb_slide_df$hap2%in%c(testhapA2, testhapB1)&
                                    hyb_slide_df$pos==df_comp$pos[i],'identity']<identity_cutoff)
  comp3 = as.numeric(hyb_slide_df[hyb_slide_df$hap1%in%c(testhapA1, testhapB2) & 
                                    hyb_slide_df$hap2%in%c(testhapA1, testhapB2)&
                                    hyb_slide_df$pos==df_comp$pos[i],'identity']<identity_cutoff)
  comp4 = as.numeric(hyb_slide_df[hyb_slide_df$hap1%in%c(testhapA2, testhapB2) & 
                                    hyb_slide_df$hap2%in%c(testhapA2, testhapB2)&
                                    hyb_slide_df$pos==df_comp$pos[i],'identity']<identity_cutoff)
  df_comp$count[i] <- comp1+comp2+comp3+comp4
  print(i/nrow(df_comp))
}

figure2 <- ggplot(df_comp, aes(x=count))+
  geom_bar(stat='count', fill='grey', color='black')+
  scale_y_continuous(limits=c(0,1200))+
  geom_hline(yintercept = nrow(df_comp), linetype='dotted')+
  geom_hline(yintercept = (nrow(df_comp)/2), linetype='dotted')+
  geom_hline(yintercept = (nrow(df_comp)/4), linetype='dotted')+
  xlab(' ')+
  theme_classic()+
  theme(text=element_text(size=20))
figure2

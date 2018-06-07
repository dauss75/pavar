library(HardyWeinberg)
library(optparse)

option_list = list(
  make_option(c("-v", "--vcf"),  type="character", help="gnomad vcf for genes"),
  make_option(c("-i", "--intervar"),  type="character", help="intervar output for pathogenic classification"),
  make_option(c("-o", "--output"),  type="character", help="Output file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(OptionParser(option_list=option_list))

if ( is.null(opt$v) & is.null(opt$i)) {
    print_help(opt_parser)
    stop("Missing input files path!\n")
} else {
    gnomad<-read.table(opt$v, fill=TRUE, sep="\t")
    paver<-read.table(opt$i, fill=TRUE, sep="\t")
}

combo<-merge(gnomad,paver,by="V2")
HWE.test<-as.data.frame(matrix(, nrow = nrow(combo), ncol = 13))
us.pop=326766748
colnames(HWE.test)<-c("Chr","Pos","GC_MM","GC_MN","GC_NN","p-val","HWE","AF","AF_male","AF_female","LCA_AF","LCA_AF_male","LCA_AF_female")
alpha=0.05
for (i in 1:nrow(combo)){
  string=combo$V8.x[i]
  t<-as.numeric(unlist(strsplit(gsub('.*GC=(.*);AF_raw=.*','\\1',string),",")))
  AF<-as.numeric(unlist(strsplit(gsub('.*;AF=(.*);AN=.*','\\1',string),",")))
  lca.AF<-floor(AF[1]*us.pop)
  AF_male<-as.numeric(unlist(strsplit(gsub('.*AF_Male=(.*);AF_Female=.*','\\1',string),",")))
  lca.AF_male<-floor(AF_male[1]*(us.pop/2))
  AF_female<-as.numeric(unlist(strsplit(gsub('.*;AF_Female=(.*);GC_Male=.*','\\1',string),",")))
  lca.AF_female<-floor(AF_female[1]*(us.pop/2))

  HWE.test[i,1:2]<-c(combo$V1.y[i],combo$V3.y[i])
  HWE.test[i,3:5]<-c(MM=t[1],MN=t[2],NN=t[3]); x<-c(MM=t[1],MN=t[2],NN=t[3])
  # chi-square test without Yates’ continuity correction, which is not recommended for low minor allele frequencies
  HWE.test[i,6] <- HWChisq(x, cc = 0, verbose = FALSE)$pval
  if (HWE.test[i,6] < alpha){
    HWE.test[i,7]<-"FALSE"
  } else {
    HWE.test[i,7]<-"TRUE"
  }
  HWE.test[i,8:10]<-c(AF[1],AF_male[1],AF_female[1])
  HWE.test[i,11:13]<-c(lca.AF,lca.AF_male,lca.AF_female)
}

write.table(HWE.test, opt$o, sep="\t", quote=FALSE, row.names=FALSE)

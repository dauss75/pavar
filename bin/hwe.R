packages <- c("data.table", "HardyWeinberg", "optparse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos = "http://cran.us.r-project.org")
}
suppressPackageStartupMessages(library(HardyWeinberg))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

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
  gnomad<-fread(opt$v); paver<-fread(opt$i)
  if(nrow(paver)==0){
    quit()
  }
}

## temporary variables
# gnomad<-fread("/Users/sjung/Documents/GitHub/pavar/output/gnomad/wgs/chr1/input_vcf/CRB1.vcf");
# paver<-fread("/Users/sjung/Documents/GitHub/pavar/output/gnomad/wgs/chr1/CRB1.pathogenic.intervar.txt")

##

colnames(gnomad)[1:5]<-c("Chr","Start","dbsnp","Ref","Alt")
combo<-merge(gnomad,paver,by="Start")

# HWE.test<-as.data.frame(matrix(, , ncol = 15))
HWE.test<-as.data.frame(matrix(, , ncol = 22))
us.pop=326766748
colnames(HWE.test)<-c("Chr","Pos","GC_MM","GC_MN","GC_NN","p-val","HWE","AF","AF_male","AF_female","LCA_AF",
                      "LCA_AF_male","LCA_AF_female","GC_AFR", "GC_AMR","GC_ASJ","GC_EAS","GC_FIN","GC_NFE","GC_OTH",
                      "q=Sum(AF)","Estimated LCA")
alpha=0.05
for (i in 1:nrow(combo)){
  if (length(unlist(strsplit(combo$Alt.x[i],","))) > 1){
    # find position of the allele
    string=combo$V8[i]
    index=which(unlist(strsplit(combo$Alt.x[i],",")) == combo$Alt.y[i])
    
    GC<-unlist(strsplit(gsub('.*;GC=(.*);.*','\\1',string),";"))[1];
    GC_content=as.numeric(unlist(strsplit(GC,","))[(3*(index-1)+1):(3*(index-1)+3)])
    AF<-as.numeric(unlist(strsplit(gsub('.*;AF=(.*);AN=.*','\\1',string),",")))[index]
    # print(AF)
    lca.AF<-floor(AF[1]*us.pop)
    AF_male<-as.numeric(unlist(strsplit(gsub('.*AF_Male=(.*);AF_Female=.*','\\1',string),",")))[index]
    lca.AF_male<-floor(AF_male[1]*(us.pop/2))
    AF_female<-as.numeric(unlist(strsplit(gsub('.*;AF_Female=(.*);GC_Male=.*','\\1',string),",")))[index]
    lca.AF_female<-floor(AF_female[1]*(us.pop/2))
    
    HWE.test[i,1:2]<-c(combo$Chr[i],combo$Start[i])
    HWE.test[i,3:5]<-c(MM=GC_content[1],MN=GC_content[2],NN=GC_content[3]); x<-c(MM=GC_content[1],MN=GC_content[2],NN=GC_content[3])
    # print(x)
    
    # chi-square test without Yates’ continuity correction, which is not recommended for low minor allele frequencies
    if (sum(x) > 0 & x[1]>5000){  # temporary value
      HWE.test[i,6] <- HWChisq(x, cc = 0, verbose = FALSE)$pval
      if (HWE.test[i,6] < alpha){
        HWE.test[i,7]<-"FALSE"
      } else {
        HWE.test[i,7]<-"TRUE"
      }
      HWE.test[i,8:10]<-c(AF[1],AF_male[1],AF_female[1])
      HWE.test[i,11:13]<-c(lca.AF,lca.AF_male,lca.AF_female)
      HWE.test[i,14]<-unlist(strsplit(gsub('.*;GC_AFR=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,15]<-unlist(strsplit(gsub('.*;GC_AMR=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,16]<-unlist(strsplit(gsub('.*;GC_ASJ=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,17]<-unlist(strsplit(gsub('.*;GC_EAS=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,18]<-unlist(strsplit(gsub('.*;GC_FIN=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,19]<-unlist(strsplit(gsub('.*;GC_NFE=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,20]<-unlist(strsplit(gsub('.*;GC_OTH=(.*);.*','\\1',string),";"))[1];

    } else {
      HWE.test<-HWE.test[-i,]
      next
    }
    
  } else {
    # find position of the allele
    string=combo$V8[i]
    GC<-unlist(strsplit(gsub('.*;GC=(.*);.*','\\1',string),";"))[1];
    GC_content=as.numeric(unlist(strsplit(GC,",")))
    AF<-as.numeric(unlist(strsplit(gsub('.*;AF=(.*);AN=.*','\\1',string),";")))
    # print(AF)
    lca.AF<-floor(AF[1]*us.pop)
    AF_male<-as.numeric(unlist(strsplit(gsub('.*AF_Male=(.*);AF_Female=.*','\\1',string),";"))[1])
    lca.AF_male<-floor(AF_male[1]*(us.pop/2))
    AF_female<-as.numeric(unlist(strsplit(gsub('.*;AF_Female=(.*);GC_Male=.*','\\1',string),";"))[1])
    GC_population<-gsub('.*;(GC_AFR=.*);GC_Male=.*','\\1',string)
    lca.AF_female<-floor(AF_female[1]*(us.pop/2))
    
    HWE.test[i,1:2]<-c(combo$Chr[i],combo$Start[i])
    HWE.test[i,3:5]<-c(MM=GC_content[1],MN=GC_content[2],NN=GC_content[3]); x<-c(MM=GC_content[1],MN=GC_content[2],NN=GC_content[3])
    
    # chi-square test without Yates’ continuity correction, which is not recommended for low minor allele frequencies
    if (sum(x) > 0 & x[1]>1000){
      HWE.test[i,6] <- HWChisq(x, cc = 0, verbose = FALSE)$pval
      if (HWE.test[i,6] < alpha){
        HWE.test[i,7]<-"FALSE"
      } else {
        HWE.test[i,7]<-"TRUE"
      }
      HWE.test[i,8:10]<-c(AF[1],AF_male[1],AF_female[1])
      HWE.test[i,11:13]<-c(lca.AF,lca.AF_male,lca.AF_female)
      HWE.test[i,14]<-unlist(strsplit(gsub('.*;GC_AFR=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,15]<-unlist(strsplit(gsub('.*;GC_AMR=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,16]<-unlist(strsplit(gsub('.*;GC_ASJ=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,17]<-unlist(strsplit(gsub('.*;GC_EAS=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,18]<-unlist(strsplit(gsub('.*;GC_FIN=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,19]<-unlist(strsplit(gsub('.*;GC_NFE=(.*);.*','\\1',string),";"))[1];
      HWE.test[i,20]<-unlist(strsplit(gsub('.*;GC_OTH=(.*);.*','\\1',string),";"))[1];
    } else {
      HWE.test<-HWE.test[-i,]
      next
    }
    
  }
}
#print(na.omit(HWE.test[,8]))
sig_q <- sum(na.omit(HWE.test[,8]))
# print(sig_q)
LCA <- floor(sig_q*sig_q*us.pop)
# print(carrier)
HWE.test[1,21:22]<-c(sig_q,LCA)
fwrite(HWE.test, opt$o, sep="\t")

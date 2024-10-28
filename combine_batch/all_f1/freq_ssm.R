setwd("~/Desktop/normal/revisit_freq_filtering_May/all_f1/")
# dir <- list.dirs()[-c(1:3,18:24)]
# dir <- gsub("\\.\\/","",dir)
#vaf_cutoff = 7
require(ggplot2)
library(data.table)
library(plyr)

files <- list.files(pattern =  "_*.ssm")
files <- files[-c(1,36)]
gene <- c()
loci <- c()

for (i in 1:length(files)) {
  

temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)

locus <- temp_file$locus
names(locus) <- temp_file$name


gene <- c(gene,unique(temp_file$name))
loci <- c(loci,locus)

  
}

a <- data.frame(loci =loci,gene=names(loci))
a_gene <- data.frame(table(gene))

#b <- data.frame(table(gene))
b <- ddply(a,.(loci,gene),nrow)




#mannually checked all the mutations in all 34 patients, remove these
remove <- c("chr7 92244630","chr2 48033272","chr13 32936644","chr1 156785616","chr17 37879761","chr19 15281458",
            "chr22 24167512","chr14 105239893","chr7 55233088","chr8 90982555","chr2 47693958",
            "chr2 48032882", "chr19 15273381","chr11 108163562","chr17 7579471",
            "chr1 12047799","chr19 15273380","chr6 152332983","chr2 48032882","chr20 36030937",
            "chr17 29508774","chr19 1218493","chr19 15281458","chr2 48010557")



for (i in 1:length(files)) {
  fn <- gsub("_freq.ssm","",files[i])
  temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)
  temp_file <- temp_file[!temp_file$locus %in% remove,]
  temp_file <- temp_file[!temp_file$name == "POLE",]
  temp_file <- temp_file[!temp_file$name == "ATR",]
  temp_file$id <- paste("s",seq(0,nrow(temp_file)-1,1),sep = "")
  write.table(temp_file,file = paste(fn,".ssm",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
  
  json <- readLines(paste("~/Desktop/normal/jan_final/all_params/",fn,".name.json",sep = ""))
  names <- gsub(".*samples","",json)
  names  <- gsub("clusters.*","",names)
  names <- unlist(strsplit(names,split = ","))
  names  <- gsub("[^[:alnum:]]","",names)
  names <- names[-length(names)]
  sink(paste(fn,".name.json",sep = ""))
  
  cat('{"samples": [')
  for (j in 1:length(names)) {
    #  name <- paste(gsub("\\..*","",files[i]))
    cat('"')
    cat(names[j])
    if(j==length(names)){
      cat('"')
      break
    }else{
      cat('",')
    }
  }
  
  
  #updating for mutation tree
  cat('], "clusters": [',paste("[",'"',temp_file$id,'"',"]",sep = "",collapse = ","),'], "garbage": []}')
  sink() 
  
}


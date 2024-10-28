setwd("~/Desktop/normal/full")
dir <- list.dirs()[-c(1:3,18:24)]
dir <- gsub("\\.\\/","",dir)
#vaf_cutoff = 7
library(ggplot2)
library(data.table)



for (d in 1:length(dir)) {
  
  ll=dir[d]
  
  setwd(paste("~/Desktop/normal/full/",ll,sep = ""))
  #library(data.table)
  files <- list.files(pattern =  "*.tsv")
cnv <- c()
  
  
  for (i in 1:length(files)) {
    #show file(sample) name
    name <- gsub("\\..*","",files[i])
    
    cat("sample: ")
    cat(paste(gsub("\\..*","",files[i]),"\n"))

    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)
    

    if(length(grep("full",name))>0){
      
      next
    }
    
    
    temp_file$locus <- paste(temp_file$chr,temp_file$start)
    temp_file$VAF <- temp_file$VAF/100
    
    

    temp_cna <- temp_file[temp_file$ISCN>=4 | temp_file$ISCN <=1,]
    if(nrow(temp_cna)>0){
      cnv <- c(cnv,unique(temp_cna$gene))
      
    }
    
    
    
  }
  cnv <- unique(cnv)
  if(length(cnv)==0){
    print(paste(name,"no cnv!"))
    next
  }else{
    print(paste(name,"yes cnv!"))
  }
  
  
  
  l <- length(files)
  names <- c()
  no <- c()
  count <- 0
  
  
  for (i in 1:length(files)) {
    #show file(sample) name
    name <- gsub("\\..*","",files[i])
    

    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)

    if(length(grep("full",name))>0){
      temp_file$id <- name
      temp_file <- temp_file[temp_file$gene %in% cnv,c(14,1,2,5,6,7,8,9)]
      germ <- temp_file
      colnames(germ)[6] <- "ISCN"

      next
    }
    
    
    #print(i)
    name <- gsub("\\..*","",files[i])
    names <- c(names,name)
    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"))
    #temp_file <- temp_file[temp_file$VAF>vaf_cutoff,]
    temp_file$locus <- paste(temp_file$chr,temp_file$start)
    temp_file$id <- name
    temp_file <- temp_file[temp_file$gene %in% cnv,c(21,20,15,6,12,13,7,9)]
    if(i==2){
      all <- rbind(germ,temp_file)
    }else{
      all <- rbind(all,temp_file)
    }
    


    
    
    
    
  }


  write.table(all,file = paste(ll,"_cnv.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)

  
  
  
  
  
  
  
  
  
}








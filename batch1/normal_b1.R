setwd("~/Desktop/normal/full")
dir <- list.dirs()[-c(1:4,19:25)]
dir <- gsub("\\.\\/","",dir)
vaf_cutoff = 7
library(ggplot2)
library(data.table)





for (d in 1:length(dir)) {
  
  ll=dir[d]
  
  setwd(paste("~/Desktop/normal/full/",ll,sep = ""))
  library(data.table)
  files <- list.files(pattern =  "*.tsv")
  snps <- list()
  vaf <- list()
  cov <- list()

  for (i in 1:length(files)) {
    name <- gsub("\\..*","",files[i])
    
    cat("sample: ")
    cat(paste(gsub("\\..*","",files[i]),"\n"))
    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)
    

    if(length(grep("full",name))>0){
      
      temp_file <- temp_file[,c(1,2,5,6,8,9,11,13)] #remove unwanted columns to reduce time
      temp_file$parse <- temp_file$locus
      temp_file[temp_file$dbsnp != "",]$dbsnp<-"YES"
      if(length(unique(temp_file$dbsnp == ""))==2){temp_file[temp_file$dbsnp == "",]$dbsnp<-"NO"}
      cat("sample: ")
      cat(paste(gsub("\\..*","",files[i]),"\n"))
      cat("\n")
      cat("\n")
      if(length(grep(",",temp_file$VAF))>0){
        for (k in grep(",",temp_file$VAF)) {
          k_temp <- sum(as.numeric(gsub(".*\\=","",unlist(strsplit(temp_file$VAF[k],split = ",")))))
          temp_file$VAF[k] <- k_temp
        }
        
        temp_file$VAF <- as.numeric(temp_file$VAF)
      }
      
      #quality filter
      temp_file$VAF <- as.numeric(temp_file$VAF)
      temp_file <- temp_file[temp_file$pvalue<0.05,]
       temp_file <- temp_file[temp_file$coverage>100,]
       temp_file <- temp_file[temp_file$VAF>7,]
      temp_file$VAF <- temp_file$VAF/100
      temp_germ <- temp_file
      
      germ_snp <- temp_file$parse
      names(germ_snp) <- temp_file$gene
      snps[[i]] <- germ_snp
      names(snps)[i] <- "germ"
      vaf[[i]] <- temp_file$VAF
      names(vaf)[i] <- "germ"
      cov[[i]] <- temp_file$coverage
      names(cov)[i] <- "germ"
      
      next
    }
    
    
    
    cat("_unfiltered total: ")
    cat(nrow(temp_file),"\n")
    temp_file$locus <- paste(temp_file$chr,temp_file$start)
    temp_file <- temp_file[temp_file$type != "CNV",]
    temp_file$parse <- temp_file$locus
    temp_file$VAF <- temp_file$VAF/100
    
    
    #check vaf diff
    a <- merge(temp_file,temp_germ,by="locus")
    a$diff <-a$VAF.x-a$VAF.y
   if(nrow(a[a$diff>0.3 && a$VAF.y<0.2,])>0){
     add_diff <- a[a$diff>0.3 && a$VAF.y<0.2,]$locus
     
     names(add_diff) <- a[a$diff>0.3,]$gene.x
     cat("diff>0.5:","\n")
     cat(add_diff,"\n\n")
       }else{
     add_diff <- c()
   }
  
    
    temp_snp <- temp_file$parse[!temp_file$parse %in% unname(unlist(snps["germ"]))]
    
    names(temp_snp) <- temp_file$gene[!temp_file$parse %in% unname(unlist(snps["germ"]))]
    
    temp_snp <- c(temp_snp,add_diff)
    snps[[i]] <- temp_snp
    
    names(snps)[i]<-name
    
    
    
    vaf[[i]] <- temp_file$VAF
    names(vaf)[i] <- name
    cov[[i]] <- temp_file$coverage
    names(cov)[i] <- name
    
    
    
  }
  
  #parse germline
  germ <- snps[["germ"]]
  
  snps["germ"] <- NULL
  germ_vaf <- vaf[["germ"]]
  vaf["germ"] <- NULL
  
  pdf(file=paste(ll,"pdf",sep = "."))
  
  #plot overall info
  for (v in 1:length(vaf)) {
    if(v==1){
      vaf_d <- data.frame(VAF=vaf[[v]],label=names(vaf)[v])
    }else{
      vaf_temp <- data.frame(VAF=vaf[[v]],label=names(vaf)[v])
      vaf_d<-rbind(vaf_d,vaf_temp)
    }
  }
  
  print(ggplot(vaf_d, aes(x=VAF, color=label)) +
          geom_density()+ggtitle(paste("VAF density plot for patient",ll)))
  
  print(ggplot(vaf_d, aes(x=VAF,color=label,fill=label)) +
          geom_histogram(binwidth = 0.02, position="dodge",alpha=0.5)+ggtitle(paste("VAF hist plot for patient",ll)))
  
  for (c in 1:length(cov)) {
    if(c==1){
      cov_d <- data.frame(COV=cov[[c]],label=names(cov)[c])
    }else{
      cov_temp <- data.frame(COV=cov[[c]],label=names(cov)[c])
      cov_d<-rbind(cov_d,cov_temp)
    }
  }
  
  print(ggplot(cov_d, aes(x=COV, color=label)) +
          geom_density()+ggtitle(paste("coverage density plot for patient",gsub("B.*","",name))))
  
  
  dev.off()
  
  germ_vaf <- vaf[["germ"]]
  vaf["germ"] <- NULL
  

  
  temp_allsnp <- data.frame(loci=gsub("\\#.*","",unlist(snps)),gene=gsub(".*\\.","",names(unlist(snps))),stringsAsFactors = F)
  temp_allsnp <- temp_allsnp[!duplicated(temp_allsnp$loci),]

  allsnp <- as.character(temp_allsnp$loci)
  names(allsnp) <- temp_allsnp$gene

  
  
  l <- length(allsnp)-1
  final <- data.frame(id=paste("s",seq(0,l,1),sep = ""),
                      locus = allsnp,
                      name=names(allsnp),
                      var_reads="",
                      total_reads="",
                      var_read_prob=paste(rep("0.5",length(files)-1),collapse = ","),
                      vaf="",
                      pval="")
  final <- final[order(final$locus),]
  names <- c()
  no <- c()
  count <- 0
  
  
  for (i in 1:length(files)) {
    #show file(sample) name
    name <- gsub("\\..*","",files[i])
    

    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)

    if(length(grep("full",name))>0){
      temp_file <- temp_file[,c(1,2,5,6,8,9,11,13)]
      temp_file$parse <- temp_file$locus
      temp_file[temp_file$dbsnp != "",]$dbsnp<-"YES"
      if(length(unique(temp_file$dbsnp == ""))==2){temp_file[temp_file$dbsnp == "",]$dbsnp<-"NO"}
      cat("sample: ")
      cat(paste(gsub("\\..*","",files[i]),"\n"))
      cat("\n")
      cat("\n")
      
      if(length(grep(",",temp_file$VAF))>0){
        for (k in grep(",",temp_file$VAF)) {
          k_temp <- sum(as.numeric(gsub(".*\\=","",unlist(strsplit(temp_file$VAF[k],split = ",")))))
          temp_file$VAF[k] <- k_temp
        }
        
        temp_file$VAF <- as.numeric(temp_file$VAF)
      }
      
      
      temp_file$VAF <- as.numeric(temp_file$VAF)
      temp_file <- temp_file[temp_file$pvalue<0.05,]
      temp_file <- temp_file[temp_file$coverage>100,]
      temp_file <- temp_file[temp_file$VAF>7,]
      temp_file$VAF <- temp_file$VAF/100
      temp_germ <- temp_file
      temp_germ <- temp_germ[,c(1,2,3,5,6,7)]
      next
    }
    
    
    print(i)
    name <- gsub("\\..*","",files[i])
    names <- c(names,name)
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"))
    temp_file$locus <- paste(temp_file$chr,temp_file$start)
    n <- length(temp_file$locus[temp_file$locus %in% allsnp])
    temp_file$`dbSNP151_common_20180423=1` <- gsub("\\_.*","",temp_file$`dbSNP151_common_20180423=1`)
    temp_file$`GNOMAD_EXOME_r2.1.1=1` <- gsub("\\_.*","",temp_file$`GNOMAD_EXOME_r2.1.1=1`)
    
    extras <- temp_file[temp_file$`dbSNP151_common_20180423=1`=="NO" & temp_file$`GNOMAD_EXOME_r2.1.1=1`=="NO",]
    extras <- extras[extras$pvalue<0.05,]
    extras <- extras[extras$coverage>100,]
    extras <- extras[,c("locus","gene","type","coverage","VAF","ISCN" )]
    extras$VAF <- extras$VAF/100
    
    
    extras <- merge(extras,temp_germ,by="locus")
    extras_vaf <- extras[extras$VAF.x-extras$VAF.y>0.3,]
    if(nrow(extras_vaf)>0){
      ex<-extras_vaf$locus
      names(ex)<- extras_vaf$gene.x
      cat("---not in db and diff>30: ")
      cat(ex,"\n\n")
      
    }
    
    if(length(temp_file$locus[temp_file$locus %in% allsnp]) != 0){
      
      parse_file <- temp_file[temp_file$locus %in% allsnp,]
      parse_file$altread <- round(parse_file$VAF*parse_file$coverage*0.01,0)
      parse_file <- parse_file[,c("locus","gene","coverage","altread","VAF","ISCN","pvalue" )]
      parse_file$VAF <- parse_file$VAF/100
      
    }else{
      no <- c(no,i)
      next
    }
    
    
    
    # assign(name,parse_file)
    if(length(unique(allsnp %in% temp_file$locus)) == 2){
      temp <- data.frame(locus=allsnp[!allsnp %in% parse_file$locus],
                         gene= names(allsnp[!allsnp %in% parse_file$locus]),
                         coverage=0,
                         altread=1,
                         VAF=0,
                         ISCN=2,
                         pvalue=0.04
                         #silent="NA",
                         #funct="NA"
      )
      all <- rbind(parse_file,temp)
      all <- all[order(all$locus),]
    }else{
      all <- parse_file
      all <- all[order(all$locus),]
    }
    
    
    all <- all[!duplicated(all$locus),]
    
    count <- count +1
    
    if(count==1){
      #print(final$var_reads)
      final$var_reads <- all$altread
      final$total_reads <- all$coverage
      final$vaf <- all$VAF
      final$pval <- all$pvalue
    }else{
      
      
      final$var_reads <- paste(final$var_reads,all$altread,sep = ",")
      final$total_reads <- paste(final$total_reads,all$coverage,sep = ",")
      final$vaf <- paste(final$vaf,all$VAF,sep = ",")
      final$pval <- paste(final$pval,all$pvalue,sep = ",")
    }

    
  }
  final$var_read_prob <- paste(rep("0.5",length(files)-length(no)-1),collapse = ",")
  cat('filter by normal:',nrow(final),"\n")
  
  #filtering by vaf
  vaf_del<-c()
  for (p in 1:nrow(final)) {
    if(length(grep("\\,",final[1,4]))>0){test_vaf <- unique(unlist(strsplit(final$vaf[p],split = ","))<0.07)
    }else{test_vaf <- final$vaf[p] <0.07}
    
    
    
    if(isTRUE(test_vaf)){
      vaf_del <- c(vaf_del,as.character(final$locus[p]))
      
    }
    
  }
  
  
  final <- final[!final$locus %in% vaf_del,]
  cat('filter by VAF:',nrow(final),"\n")
  
  #filtering by coverage
  cov_del<-c()
  for (p in 1:nrow(final)) {
    if(length(grep("\\,",final[1,4]))>0){test_cov <- unique(as.numeric(unlist(strsplit(final$total_reads[p],split = ",")))<100)
    }else{test_cov <-final$total_reads[p]<100}
    if(isTRUE(test_cov)){
      cov_del <- c(cov_del,as.character(final$locus[p]))
      
    }
    
  }
  final <- final[!final$locus %in% cov_del,]
  cat('filter by coverage normal:',nrow(final),"\n")
  
  #filtering by p val
  if(length(grep("\\,",final[1,4]))>0){
    cov_del<-c()
    for (p in 1:nrow(final)) {
      test_cov <- unique(as.numeric(unlist(strsplit(final$pval[p],split = ",")))>0.05)
      if(length(test_cov) == 2){
        cov_del <- c(cov_del,as.character(final$locus[p]))
        
      }
      
    }
    final <- final[!final$locus %in% cov_del,]
  }else{final <- final[final$pval<0.05,]}
  cat('filter by p-value:',nrow(final),"\n")
  final <- final[,-c(7,8)]
  
  
  #add coverage 0
  if(length(grep("\\,",final[1,4]))>0){
    for (p in 1:nrow(final)) {
      covs <- as.numeric(unlist(strsplit(final$total_reads[p],split = ",")))
      if(length(unique(covs==0))==2){
        covs[covs==0] <- round(mean(covs[covs!=0]),0)
        final$total_reads[p] <-  paste(covs,collapse = ",") 
      }
    }
  }
  
  final$id <- paste("s",seq(0,nrow(final)-1,1),sep = "")
  write.table(final,file = paste(ll,".ssm",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
  sink(paste(ll,".name.json",sep = ""))
  if(is.null(no)==FALSE){
    new_files <- files[-no]
  }else{
    new_files <- files
  }
  
  new_files <- new_files[-grep("full",new_files)]
  
  cat('{"samples": [')
  for (i in 1:length(new_files)) {
    #  name <- paste(gsub("\\..*","",files[i]))
    cat('"')
    cat(paste(gsub(".tsv","",new_files[i])))
    if(i==length(new_files)){
      cat('"')
      break
    }else{
      cat('",')
    }
  }
  #updating for mutation tree
  cat('], "clusters": [',paste("[",'"',final$id,'"',"]",sep = "",collapse = ","),'], "garbage": []}')
  sink() 
  
  
  
  
  
  
  
  
  
}








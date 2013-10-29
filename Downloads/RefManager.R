# Mathew McLean
# Functions defining a citation manager for R
# January 30, 2013

add.refDB <- function(refdf=NULL,bibEntry,location,filename='',topic){
  
  # create new database of references if one not provided
  if(is.null(refdf)){
    refdf <- data.frame(matrix(nr=1000,nc=9))
    colnames(refdf) <- c('title','author','journal','year','bibtype','topic',
                           'location','filename','bibtexEntry')
    attr(refdf,'size') <- 0
    class(refdf) <- c('refDB','data.frame')
  }
  
  initproc <- strsplit(bibEntry,'(\n)[[:space:]]*')[[1]]
  bibtype <- strsplit(substring(initproc[1],2),'[{]')[[1]][1]
  numentries <- length(initproc)-1
  fields <- numeric(numentries-1)
  vals <- numeric(numentries-1)
  for(i in 2:numentries){
    temp <- strsplit(initproc[i],'[=]')[[1]]
    fields[i-1] <- toupper(temp[1])
    vals[i-1] <- strsplit(substring(temp[2],2),'[}]')[[1]][1]
  }
  ind <- attr(refdf,'size') + 1
  
  tind <- which(fields%in%'AUTHOR')[1]
  if(!is.na(tind)) refdf$author[ind] <- vals[tind]
  tind <- which(fields%in%'TITLE')[1]
  if(!is.na(tind)) refdf$title[ind] <- vals[tind]
  tind <- which(fields%in%'JOURNAL')[1]
  if(!is.na(tind)) refdf$journal[ind] <- vals[tind]
  tind <- which(fields%in%'YEAR')[1]
  if(!is.na(tind)) refdf$year[ind] <- vals[tind]
  
  refdf$bibtype[ind] <- bibtype
  if(!missing(location)) refdf$location[ind] <- location
  refdf$filename[ind] <- filename
  if(!missing(topic)) refdf$topic[ind] <- topic
  refdf$bibtexEntry[ind] <- bibEntry
  attr(refdf,'size') <- ind
  
  return(refdf)
}

search.refDB <- function(dtb,searchterms,searchfields){
  size <- attr(dtb,'size')
  if(!all(searchfields%in%colnames(dtb))){
    stop(cat('Some of the supplied search fields are not valid.  Must be one of',colnames(dtb)))
  }
  if(length(searchterms)!=length(searchfields)){
    stop('searchterms and searchfields must be the same length')
  }
  nst <- length(searchterms)
  ind <- 1:size
  for(i in 1:nst){
    temp <- grep(toupper(searchterms[i]),toupper(dtb[[searchfields[i]]][ind]))
    if(!length(temp)){
      print('No match found')
      ind <- temp
      break
    }
    ind <- ind[temp]
  }
  if(length(ind)){
    print.refDB(dtb,ind)
  }
  return(ind)
}

open.refDB <- function(dtb,entrynum){
  if(is.na(dtb$location[entrynum]) | is.na(dtb$filename[entrynum])){
    stop(paste('No location or filename information for database entry',entrynum,sep=''))
  }
  location <- gsub('[[:space:]]','%20',dtb$location[entrynum])
  if(dtb$filename[entrynum]==''){
    shell.exec(location)
  }else{
    system(paste('open file:///',location,dtb$filename[entrynum],sep=''))
  }
}

print.refDB <- function(dtb,entrynums){
  cnames <- toupper(colnames(dtb))
  size <- attr(dtb,'size')
  if(missing(entrynums)) entrynums <- 1:size
  
  for(i in 1:length(entrynums)){
    if(entrynums[i]>size){
      stop(paste('Invalid entrynum=',entrynums[i,],'. entrynum must be less than',size,sep=''))
    }
    cat('ENTRY #   :    ',entrynums[i],'\n')
    for(j in 1:(ncol(dtb)-1)){    
      cat(cnames[j],rep('',9-nchar(cnames[j])),':    ',dtb[entrynums[i],j],'\n')
    }
    cat('////////////////////////////////////////////////////////////////////////////////\n')
  }
}

getbibtex.refDB <- function(dtb,entrynums){
  size <- attr(dtb,'size')
  if(any(entrynums>size)){
    stop('Invalid entrynum provided')
  }
  for(i in 1:length(entrynums)){

    cat(dtb$bibtexEntry[entrynums[i]],'\n')
  }
}

remove.refDB <- function(dtb,entrynum){
  size <- attr(dtb,'size')
  if(entrynum>size){
    stop('Invalid entrynum provided')
  }
  if(entrynum==size){
    dtb[entrynum:size,] <- rep(NA,9)
  }else{
    dtb[entrynum:size,] <- rbind(dtb[(entrynum+1):size,],rep(NA,9))
  }
  attr(dtb,'size') <- size-1
  
  return(dtb)
}
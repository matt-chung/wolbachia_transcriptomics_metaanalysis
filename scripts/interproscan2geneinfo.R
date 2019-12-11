args <- commandArgs(TRUE)
interproscan.path <- args[1]
counts.path <- args[2]

# interproscan.path <- "Z:/EBMAL/mchung_dir/EHPYL/NC_000915.1.cds.fa.interproscan.tsv"
# counts.path <- "Z:/EBMAL/mchung_dir/EHPYL/transcriptome_analysis/hpylori_counts.tsv"
# 
gomap.path <- "/home/mattchung/maps/goid2godescription.map"
iprmap.path <- "/home/mattchung/maps/interproid2description.map"

interproscan <- read.delim(interproscan.path,
                           header = F,
                           col.names = seq(1,14),
                           sep = "\t")
counts <- read.delim(counts.path,
                     sep = "\t")
gomap <- read.delim(gomap.path,
                    header = F,
                    sep = "\t")
gomap[,1] <- as.character(gomap[,1])
gomap[,2] <- as.character(gomap[,2])
gomap[,3] <- as.character(gomap[,3])

iprmap <- read.delim(iprmap.path,
                     header = F,
                     sep = "\t")
iprmap[,1] <- as.character(iprmap[,1])
iprmap[,2] <- as.character(iprmap[,2])

interproscan[,1] <- gsub(".[^_]+$","",interproscan[,1])

geneinfo <- as.data.frame(matrix(nrow=nrow(counts),
                                 ncol=5))
colnames(geneinfo) <- c("gene","interpro_description","go_biologicalprocess","go_cellularcomponent","go_molecularfunction")

goterms <- sort(unique(unlist(strsplit(paste(interproscan[,14],collapse="|"),split="[|]"))))
iprterms <- sort(unique(unlist(strsplit(paste(interproscan[,12],collapse="|"),split="[|]"))))

missinggoterms <- goterms[!(goterms %in% gomap[,1]) & goterms != "" ]
missingiprterms <- iprterms[!(iprterms %in% iprmap[,1]) & iprterms != ""]

while(length(missinggoterms) > 0){
  go <- readLines(paste0("https://gowiki.tamu.edu/wiki/index.php/Category:",missinggoterms[1]),warn=F)
  gocategory <- grep("namespace:",go,value=T)
  gocategory <- gsub(".* ! ","",gocategory)
  gocategory <- gsub("\".*","",gocategory)
  
  godescription <- grep("wgTitle",go,value=T)
  godescription <- gsub(".* ! ","",godescription)
  godescription <- gsub("\",","",godescription)
  
  gomap <- as.data.frame(rbind(gomap,
                               c(missinggoterms[1],gocategory,godescription)))
  
  missinggoterms <- missinggoterms[-1]
}

while(length(missingiprterms) > 0){
  interprodescription <- readLines(paste0("https://www.ebi.ac.uk/interpro/entry/",missingiprterms[1]))
  interprodescription <- grep("h2 class", interprodescription,value = T)
  interprodescription <- gsub(".*<h2 class=\"strapline\">","",interprodescription)
  interprodescription <- gsub(" <span>.*","",interprodescription)
  iprmap <- as.data.frame(rbind(iprmap,
                                c(missingiprterms[1],interprodescription)))
  
  missingiprterms <- missingiprterms[-1]
}

write.table(gomap,
            gomap.path,
            quote = F,
            col.names = F,
            row.names = F,
            sep= "\t")
write.table(iprmap,
            iprmap.path,
            quote = F,
            col.names = F,
            row.names = F,
            sep= "\t")

missing_from_map <- c()
geneinfo$gene <- rownames(counts)
for(i in 1:nrow(counts)){
  interproscan.subset <- interproscan[interproscan[,1] == rownames(counts)[i],]
  goterms <- sort(unique(unlist(strsplit(paste(interproscan.subset[,14],collapse="|"),split="[|]"))))
  iprterms <- sort(unique(unlist(strsplit(paste(interproscan.subset[,12],collapse="|"),split="[|]"))))
  
  goterms <- goterms[goterms != ""]
  iprterms <- iprterms[iprterms != ""]
  while(length(iprterms) > 0){
    geneinfo$interpro_description[i] <- paste0(geneinfo$interpro_description[i],"|",iprterms[1],":",iprmap[iprmap[,1] == iprterms[1],2])
    if(length(iprmap[iprmap[,1] == iprterms[1],2]) == 0){
      missing_from_map <- c(missing_from_map,iprterms[1])
    }
    iprterms <- iprterms[-1]
  }
  while(length(goterms) > 0){
    if(gomap[gomap[,1] == goterms[1],2] == "biological_process"){
      geneinfo$go_biologicalprocess[i] <- paste0(geneinfo$go_biologicalprocess[i],"|",goterms[1],":",gomap[gomap[,1] == goterms[1],3])
    }else if(gomap[gomap[,1] == goterms[1],2] == "cellular_component"){
      geneinfo$go_cellularcomponent[i] <- paste0(geneinfo$go_cellularcomponent[i],"|",goterms[1],":",gomap[gomap[,1] == goterms[1],3])
    }else if(gomap[gomap[,1] == goterms[1],2] == "molecular_function"){
      geneinfo$go_molecularfunction[i] <- paste0(geneinfo$go_molecularfunction[i],"|",goterms[1],":",gomap[gomap[,1] == goterms[1],3])
    }else{
      missing_from_map <- c(missing_from_map,goterms[1])
    }
    goterms <- goterms[-1]
  }
}
missing_from_map <- unique(missing_from_map)
print(missing_from_map)

for(i in 1:ncol(geneinfo)){
  geneinfo[,i] <- gsub("NA[|]","",geneinfo[,i])
}

write.table(geneinfo,
            gsub("[.]tsv$",".geneinfo.tsv",interproscan.path),
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")

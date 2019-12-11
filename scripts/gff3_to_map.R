args <- commandArgs(TRUE)
gff.path <- args[1]
#gff.path <- "X:/mattchung/EBMAL/wBm/references/GCF_000008385.1_ASM838v1_genomic.gff"

gff <- read.delim(gff.path,
                  comment.char = "#",
                  header = F,
                  sep = "\t")
map <- unique(gff[,c(1,4,5,7)])

id <- unique(substr(unlist(strsplit(paste0(gff[,9],collapse=";"),split=";")),
                    1,
                    regexpr("=",unlist(strsplit(paste0(gff[,9],collapse=";"),split=";"))) - 1
                    ))

map <- as.data.frame(cbind(map,
                           as.data.frame(matrix(nrow=nrow(map),
                                                ncol=length(id)))))
colnames(map) <- c("contig","start","stop","strand",id)


for(i in 1:nrow(map)){
  gff.subset <- gff[gff[,1] == map[i,1] &
                    gff[,4] == map[i,2] &
                    gff[,5] == map[i,3] &
                    gff[,7] == map[i,4],]
  functional_terms <- unique(unlist(strsplit(paste0(gff.subset[,9],collapse=";"),split=";")))
  id.subset <- substr(functional_terms,
                             1,
                             regexpr("=",functional_terms) - 1)
  term.subset <- substr(functional_terms,
                               regexpr("=",functional_terms) + 1,
                               nchar(functional_terms))
  for(j in 1:length(id.subset)){
    map[i,colnames(map) == id.subset[j]] <- paste0(map[i,colnames(map) == id.subset[j]],"|",term.subset[j])
  }
}

map <- as.data.frame(apply(map,2,function(x){gsub("^NA[|]","",x)}))

write.table(map,
            paste0(gff.path,".map"),
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")

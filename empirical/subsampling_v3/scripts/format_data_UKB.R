library(tidyverse)

args<-commandArgs(trailingOnly = TRUE)

data<-read.table(args[1],header=TRUE) %>% select(SNP,CLST,MAC,NCHROBS)
data<-data %>% rename("AC"="MAC","AN"="NCHROBS")
data<-data %>% pivot_wider(names_from=CLST,values_from=c(AC,AN))

meta<-read.table(args[3]) %>% select(V1,V3) %>% rename("Annot"="V3","SNP"="V1")

data_annot<-left_join(data,meta,by="SNP")%>% separate(col=SNP,into=c("CHROM","POS","REF","ALT"),sep=":") %>% mutate("Effect"=NA)

write_delim(data_annot,args[2],delim="\t")

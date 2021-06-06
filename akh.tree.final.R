### set up environment ----
#setwd("../AKH tree/")

#install.packages("phylotools")
#install.packages("myTAI")
#install.packages("taxize")
library(tidyverse)
library(Biostrings)
library(ape)
library(DECIPHER)
library(taxize)
library(myTAI)
library(ggtree)

### load and clean data ----
df <- phylotools::read.fasta("20210312.ncbi-protein-'Gonadotropin-releasing hormone receptor'.fasta")
CrzR <- phylotools::read.fasta("NP_648571.1-dmCrzR.fasta")

# parse fasta header -----

# accession numbers
df$accession <- sub(pattern = "\\s.*", "", df$seq.name)
CrzR$accession <- sub(pattern = "\\s.*", "", CrzR$seq.name)

# gene name
df$name <- gsub(
                pattern = "\\sisoform\\s\\w?.?\\s|LOW QUALITY PROTEIN:.|PREDICTED:.|.(pseudogene).|,.partial.|\\s$",
                replacement = "",
                gsub(
                    pattern = "\\[\\w+\\s\\w+]|\\[\\w+\\s\\w+\\s\\w+]",
                    replacement = "",
                    gsub(pattern = "^.._\\d+.\\d\\s", "", df$seq.name,perl = T)
                )
)
CrzR$name <- gsub(
    pattern = ",\\sisoform\\s\\w?.?\\s|LOW QUALITY PROTEIN:.|PREDICTED:.|.(pseudogene).|,.partial.|\\s$",
    replacement = "",
    gsub(
        pattern = "\\[\\w+\\s\\w+]|\\[\\w+\\s\\w+\\s\\w+]",
        replacement = "",
        gsub(pattern = "^.._\\d+.\\d\\s", "", CrzR$seq.name,perl = T)
    )
)

# species
df$species <- gsub( pattern = "^.*.\\s\\[|]", replacement = "", df$seq.name, perl = T)
CrzR$species <- gsub( pattern = "^.*.\\s\\[|]", replacement = "", CrzR$seq.name, perl = T)    

# aa length 
df$aa.length <- nchar(df$seq.text)
CrzR$aa.length <- nchar(CrzR$seq.text)

# isoform
df$isoform <- ifelse(
    grepl("isoform",df$seq.name),
    gsub( pattern = ".+\\sisoform\\s|\\s.*", replacement = "", df$seq.name, perl = T),
    NA
)
CrzR$isoform <- ifelse(
    grepl("isoform",CrzR$seq.name),
    gsub( pattern = ".+\\sisoform\\s|\\s.*", replacement = "", df$seq.name, perl = T),
    NA
)

# pseudogene, partial, PREDICTED:, LOW QUALITY PROTEIN:
df$pseudogene <- ifelse( grepl("pseudogene",df$seq.name), "pseudogene", NA )
df$partial <- ifelse( grepl("partial",df$seq.name), "partial", NA )
df$predicted <- ifelse( grepl("PREDICTED",df$seq.name), "predicted", NA )
df$low.quality <- ifelse( grepl("LOW QUALITY PROTEIN",df$seq.name), "Low Quality Protein", NA )

CrzR$pseudogene <- ifelse( grepl("pseudogene",CrzR$seq.name), "pseudogene", NA )
CrzR$partial <- ifelse( grepl("partial",CrzR$seq.name), "partial", NA )
CrzR$predicted <- ifelse( grepl("PREDICTED",CrzR$seq.name), "predicted", NA )
CrzR$low.quality <- ifelse( grepl("LOW QUALITY PROTEIN",CrzR$seq.name), "Low Quality Protein", NA )

# make new short name
df$name <- sub(pattern = "gonadotropin.releasing hormone", 
               replacement = "GnRH", 
               x = df$name,ignore.case = T )
df$name <- sub(pattern = "adipokinetic hormone", 
               replacement = "AKH", 
               x = df$name,ignore.case = T )
CrzR$name <- sub(pattern = "corazonin receptor", 
                   replacement = "Crz receptor", 
                   x = CrzR$name,ignore.case = T )

df <- rbind.data.frame(df,CrzR)

### load manually generated index determining which sequences to drop ----
drop.df <- read.csv("GnRHreceptors.droplist.csv")
drop.df$include <- as.logical(drop.df$include)

df <- left_join(df,drop.df,"accession")
df <- df[df$include,]
df <- df[ !is.na(df$seq.name), ]

df$species.short <- stringr::word(df$species, 1, 2)

rm(CrzR,drop.df)

# get phylum ----

#### Change this entry to your unique key  ####
key <- readLines("key_entrez.txt")
rentrez::set_entrez_key(key)

species.info <- lapply(
    unique(df$species),
    function(x) {taxonomy(organism = x, db = "ncbi", output = "classification" ) %>%
        filter( rank == "species" | rank == "family"| rank == "order"|
                rank == "class "| rank == "phylum") %>%
        select(-id) %>% pivot_wider( names_from = rank,values_from = name)}
)

species.info[[52]] <- taxonomy(organism = df$species[52], db = "ncbi", output = "classification" ) %>% 
    filter( rank == "species" | rank == "family"| rank == "order"| rank == "class "| rank == "phylum") %>% 
    select(-id) %>% pivot_wider( names_from = rank,values_from = name)
species.info <- bind_rows(species.info)

colnames(species.info)[4] <- "species.short"

species.info[ species.info$species.short == "Larimichthys crocea",]$order <- "Acanthuriformes"
species.info[ species.info$species.short == "Chrysochloris asiatica",]$order <- "Afrosoricida"
species.info[ species.info$species.short == "Echinops telfairi",]$order <- "Afrosoricida"

species.info <- rbind.data.frame(
                species.info,
                data.frame(phylum ="Chordata", order ="Perciformes", 
                    family ="Moronidae", species.short ="Morone saxatilis"))

#save(species.info,file = "species.info.RData")
#load("species.info.RData")

df$species.short <- stringr::word(df$species, 1, 2)
df <- left_join(df,species.info,by = "species.short")

### prepare for tree ----
ls <- df$seq.text
#names(ls) <- paste(df$accession,df$name, ifelse(is.na(df$isoform),"  ",df$isoform), df$species)
names(ls) <- paste(df$accession,df$species)

df$tip.name <- paste(df$accession,df$species)

set.seed(0319)
seqs <- AAStringSet(x = ls,use.names = TRUE)
seqs.align <- AlignSeqs(seqs)
seqs.dm <- DistanceMatrix(seqs.align)
seqs.njs <- njs(seqs.dm)

# Make trees ----
tip_metadata <- df %>% select( tip.name,species.short:family)
colnames(tip_metadata)[3] <- "Phylum"

p <- ggtree(seqs.njs, branch.length='none', layout='circular') +
    geom_tiplab2(color="black", offset = 1, align=TRUE,size=1.5)

mp <- p %<+% tip_metadata + 
    geom_tippoint(aes(shape=Phylum, color=Phylum, fill=Phylum),size=2,position = position_nudge(x=13)) +
    geom_hilight(node=844, fill="dodgerblue", alpha=.4) +
    geom_hilight(node=859, fill="darkgreen", alpha=.4) +
    scale_shape_manual(values = c("Arthropoda"=19,"Chordata"=19,"Brachiopoda"=19,
                    "Cnidaria"=23,"Echinodermata"=24,"Mollusca"=25,"Nematoda"=22),) +
    scale_color_manual(values = c("Arthropoda"="#003f5c","Chordata"="#ff0000","Brachiopoda"="#cd007c",
                    "Cnidaria"="#d45087","Echinodermata"="#f95d6a","Mollusca"="#d45087","Nematoda"="#ff7c43"))+
    scale_fill_manual(values = c("Arthropoda"="#003f5c","Chordata"="#ff0000","Brachiopoda"="#cd007c",
                    "Cnidaria"="#d45087","Echinodermata"="#f95d6a","Mollusca"="#d45087","Nematoda"="#ff7c43"))+
    guides(shape = guide_legend(override.aes = list(size = 6))) +
    theme(legend.title=element_text(size=28), 
          legend.text=element_text(size=24))

mp <- rotate_tree(mp,-150)

ggsave("GnRHReceptor.tree.pdf",mp,height = 24,width = 24)

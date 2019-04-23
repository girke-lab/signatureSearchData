## For making "lincs", "lincs_expr", "cmap", "cmap_expr", "cmap_rank" datasets
## from scratch, please refer to the vignette of this package by running
## `browseVignettes("signatureSearchData")` in R. 

################################
## make dtlink_db_clue_sti.db ##
################################

# This SQLite database contains drug-target links (dtlink) in DrugBank, CLUE, 
# and STITCH databases

## Get dtlink in DrugBank database (version 5.1.2)

### Download the `drugbank_all_target_uniprot_links.csv.zip` at 
### <https://www.drugbank.ca/releases/latest#external-links> under 
### `Target Drug-UniProt Links` section where `DRUG GROUP` is `All`.
### Extract the downloaded zip file, rename the extracted csv file as 
### `drug_target_uniprot_links_5.1.2.csv`, save it under `data` directory 

dt_uni_link <- read.delim("./data/drug_target_uniprot_links_5.1.2.csv",
                        sep=",", stringsAsFactors=FALSE)
#### Add gn_sym column for dt_uni_link
library(org.Hs.eg.db)
k <- keys(org.Hs.eg.db, keytype = "UNIPROT")
uni_gnsym <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys=k, columns=c("UNIPROT", "SYMBOL"), keytype="UNIPROT"))
dt_uni_link <- dplyr::as_tibble(dt_uni_link)
dt_uni_sym <- dplyr::left_join(dt_uni_link[,c("Name","UniProt.ID")], 
                               uni_gnsym, by = c("UniProt.ID"="UNIPROT"))
dtlink_db <- na.omit(unique(data.frame(drug_name=tolower(dt_uni_sym$Name), 
                        t_gn_sym=dt_uni_sym$SYMBOL, stringsAsFactors=FALSE)))
#### exclude drugs that have more than 100 targets in dtlink_db
dtlist_db <- split(dtlink_db$t_gn_sym, dtlink_db$drug_name)
dtnum_db <- sapply(dtlist_db, length)
sum(dtnum_db>100)
dtlist_db <- dtlist_db[dtnum_db<=100]
dtlink_db <- list2link(dtlist_db, type="dt")
#### exclude proteins that are targeted by more than 100 drugs in dtlink_db
tdlist_db <- split(dtlink_db$drug_name, dtlink_db$t_gn_sym)
tdnum_db <- sapply(tdlist_db, length)
sum(tdnum_db>100)
tdlist_db <- tdlist_db[tdnum_db<=100]
tdlink_db <- list2link(tdlist_db, type="td")
dtlink_db <- tdlink_db[,c("drug_name","t_gn_sym")]
length(unique(dtlink_db$drug_name)) # 5250
length(unique(dtlink_db$t_gn_sym)) # 2443
write.table(dtlink_db, "./data/dtlink_db_5.1.2.xls", 
            quote=FALSE, row.names=FALSE, sep="\t")

## Get dtlink in CLUE touchstone database at <https://clue.io/touchstone>
### Get perturbation information in CLUE touchstone (version 1.1.1.2) from API

#### Query the CLUE website at <https://clue.io/query> to get drugs in 
#### Touchstone database:
#### At Query window, select `Gene expression (L1000)` as query parameters. 
#### Select `individual query` as query mode. Then click the `example` link and
#### choose `PI3K/MTOR inhibitor` as an example query, the name of your query,
#### UP-regulated genes and DOWN-regulated genes would be automatically filled.
#### Then click `SUBMIT` button.
#### At the Query History window, select the query result, then click the 
#### `DETAILED LIST` button
#### At the result window, select `Compound x 2429` under `Perturbagen Type`.
#### Click `Export` button to save the result as `clue_example_mtor_res.txt` 
#### to your data directory under current R session. 

clue_res <- read.delim("./data/clue_example_mtor_res.txt",sep="\t") # 2429 X 6
clue_drugs <- as.character(unique(clue_res$Name)) # 2397
clue_pert_info_api <- getPertInfo(lincs_drugs)
saveRDS(clue_pert_info_api, "./data/clue_pert_info_api.rds")
clue_dt <- clue_pert_info_api[,c("pert_iname","target")]
clue_dtlist <- sapply(split(clue_dt$target, clue_dt$pert_iname), 
                       function(x) unique(unlist(x))) # 2397
dtlink_clue <- list2link(clue_dtlist, type="dt")
dtlink_clue$drug_name <- tolower(dtlink_clue$drug_name)
#### exclude drugs that have more than 100 targets, 
#### then exclude genes that have more than 100 targeted drugs
dtlist_clue <- split(dtlink_clue$t_gn_sym, dtlink_clue$drug_name)
dtnum_clue <- sapply(dtlist_clue, length)
sum(dtnum_clue > 100) # 0

tdlist_clue <- split(dtlink_clue$drug_name, dtlink_clue$t_gn_sym)
tdnum_clue <- sapply(tdlist_clue, length)
sum(tdnum_clue > 100)
tdlist_clue <- tdlist_clue[tdnum_clue <= 100]
tdlink_clue <- list2link(tdlist_clue, type="td")
dtlink_clue <- tdlink_clue[,c("drug_name","t_gn_sym")]
length(unique(dtlink_clue$drug_name)) # 1949
length(unique(dtlink_clue$t_gn_sym)) # 1325
write.table(dtlink_clue, "data/dtlink_clue_1.1.1.2.xls", 
            quote=FALSE, sep="\t", row.names=FALSE)

## Get dtlink in STITCH database (v5.0) 
### At STITCH Download page 
### <http://stitch.embl.de/cgi/download.pl?UserId=Z159fH9RYlUj&sessionId=2yGVo2kYC0ri>,
### Choose `Homo sapiens` organism, download and decompress the 
### `9606.protein_chemical.links.detailed.v5.0.tsv.gz` and
### `chemicals.v5.0.tsv.gz` as 
### `9606_protein_chemical_links_detailed_v5.0.tsv` and
### `chemicals.v5.0.tsv`, respectively, to the data directory 
### under your current working directory of R session. 

### Need to allocate large number of memory to read the table, 20Gb works
sti_dtlink_detail <- readr::read_tsv(
  "data/9606_protein_chemical_links_detailed_v5.0.tsv") # 15,473,939 X 7
sti_dtlink <- dplyr::filter(sti_dtlink_detail, 
                            experimental > 700 | database > 700) # 346,909 X 7
### Eliminate drugs that have more than 100 targets in `sti_dtlink`
sti_dtlist <- split(sti_dtlink$protein, sti_dtlink$chemical)
sti_dtnum <- sapply(sti_dtlist, length)
invalid_drugs <- names(sti_dtlist[sti_dtnum>100]) # 356
sti_dtlink <- dplyr::filter(sti_dtlink, ! chemical %in% invalid_drugs)
sti_dtlink$protein <- sub("9606.", "", sti_dtlink$protein) # 227,948 X 7
### get PubChem CID to drug name mappings from `chemicals.v5.0.tsv` file
sti_chem_info <- readr::read_tsv("data/chemicals.v5.0.tsv")
### Map Ensembl protein ids to gene symbol 
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
gnames <- keys(edb, keytype = "GENENAME")
ensb_gn_protid <- na.omit(select(edb, keys = gnames, keytype = "GENENAME", 
                                 columns = c("GENENAME", "PROTEINID")))
sti_dtlink %<>% left_join(ensb_gn_protid, by = c("protein" = "PROTEINID")) %>% 
    left_join(sti_chem_info[,c("chemical","name")]) %>% 
    dplyr::select(drug_name=name, t_gn_sym=GENENAME) # 227,948 X 2
sti_dtlink <- unique(na.omit(sti_dtlink)) # 142291 X 2
sti_dtlink$drug_name <- sapply(sti_dtlink$drug_name, 
                               function(x) unlist(strsplit(x,", "))[1])
sti_dtlink$drug_name <- tolower(sti_dtlink$drug_name)
### Exclude proteins that are targeted by more than 100 drugs in sti_dtlink
#### check whether there are drugs that have more than 100 targets 
#### after ENSP id to gene SYMBOL mappings
dtlist_sti <- split(sti_dtlink$t_gn_sym, sti_dtlink$drug_name)
dtnum_sti <- sapply(dtlist_sti, length)
sum(dtnum_sti>100) # 0
tdlist_sti <- split(sti_dtlink$drug_name, sti_dtlink$t_gn_sym)
tdnum_sti <- sapply(tdlist_sti, length)
sum(tdnum_sti>100) # 255
tdlist_sti <- tdlist_sti[tdnum_sti<=100]
tdlink_sti <- list2link(tdlist_sti, type="td")
dtlink_sti <- tdlink_sti[,c("drug_name","t_gn_sym")] # 50,012 X 2
length(unique(dtlink_sti$drug_name)) # 18,847
length(unique(dtlink_sti$t_gn_sym)) # 5,119
write.table(dtlink_sti, "./data/dtlink_sti_v5.0.xls", 
            quote=FALSE, sep="\t", row.names=FALSE)

## Get dtlink tables in DrugBank/CLUE/STITCH databases and store in 
## SQLite database
### dtlink: drug_name-t_gn_sym and drug_name-ENTREZ 
### (Note: drug_name in dtlink are all lowercase)
#### Combine dtlink in DrugBank/CLUE/STITCH databases
dtlink_db <- read.delim("data/dtlink_db_5.1.2.xls", stringsAsFactors=FALSE)
dtlink_clue<-read.delim("data/dtlink_clue_1.1.1.2.xls", stringsAsFactors=FALSE)
dtlink_sti <- read.delim("data/dtlink_sti_v5.0.xls", stringsAsFactors=FALSE)
dtlink <- rbind(dtlink_db, dtlink_clue, dtlink_sti) %>% unique
#### Exclude proteins or drugs that have more than 100 targets
dtlist <- split(dtlink$t_gn_sym, dtlink$drug_name)
tdlist <- split(dtlink$drug_name, dtlink$t_gn_sym)
dt_num <- sapply(dtlist, length)
sum(dt_num>100) # 0
td_num <- sapply(tdlist, length)
sum(td_num>100) # 29
tdlist <- tdlist[td_num<=100]
tdlink <- list2link(tdlist, type="td")
dtlink <- tdlink[,c("drug_name","t_gn_sym")] # 62,228 X 2
write.table(dtlink, "data/dtlink3.xls", quote=FALSE, row.names=FALSE, sep="\t")
length(unique(dtlink$drug_name)) # 22647
length(unique(dtlink$t_gn_sym)) # 5729
#### Map gene SYMBOL to entrez id to get drug_name-entrez links
library(org.Hs.eg.db)
sym <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
sym_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = sym, 
                                    keytype = "SYMBOL", columns = "ENTREZID")
dtlink_entrez <- left_join(as_tibble(dtlink), sym_entrez, 
                           by = c("t_gn_sym"="SYMBOL")) %>% 
                 dplyr::select(drug_name, ENTREZID) %>% 
                 distinct(.keep_all = TRUE)
write.table(dtlink_entrez, "data/dtlink_entrez3.xls", quote=FALSE, 
            row.names=FALSE, sep="\t")
### Store dtlinks in SQLite database
library(RSQLite)
mydb <- dbConnect(SQLite(), "data/dtlink_db_clue_sti.db")
dbWriteTable(mydb, "dtlink_db", dtlink_db, overwrite = TRUE)
dbWriteTable(mydb, "dtlink_clue", dtlink_clue, overwrite = TRUE)
dbWriteTable(mydb, "dtlink_sti", dtlink_sti, overwrite = TRUE)
dbWriteTable(mydb, "dtlink", dtlink, overwrite = TRUE)
dbWriteTable(mydb, "dtlink_entrez", dtlink_entrez, overwrite = TRUE)
dbListTables(mydb)
dbDisconnect(mydb)

#####################
## make goAnno.rds ##
#####################

# It is an intermediate file storing the GO to GENE SYMBOLE mapping information
# for `dtnetplot` function in `signatureSearch` package

OrgDb <- GOSemSim::load_OrgDb(OrgDb = "org.Hs.eg.db")
kk <- keys(OrgDb, keytype="SYMBOL")
goAnno <- AnnotationDbi::select(OrgDb, keys=kk, keytype="SYMBOL",
                                columns=c("GOALL", "ONTOLOGYALL"))
goAnno <- unique(na.omit(goAnno)) # 3,430,396 X 4
saveRDS(goAnno, "./data/goAnno.rds")

##########################
## make goAnno_drug.rds ##
##########################

# It is an intermediate file storing the GO to drug name mapping information
# for `dsea` functions in `signatureSearch` package

## convert GOALL-SYMBOL in goAnno to GOALL-drug
### Get drug-target mapping information in DrugBank, CLUE and STITCH databases 
### from `dtlink_db_clue_sti.db` created above
library(RSQLite); library(dplyr)
conn <- dbConnect(SQLite(), "./data/dtlink_db_clue_sti.db")
dtlink <- dbGetQuery(conn, 'SELECT * FROM dtlink') # 62,228 X 2
dbDisconnect(conn)
### Convert GO-GENE mappings to GO-drug mappings
goAnno <- readRDS("./data/goAnno.rds")
goAnno_drug <- left_join(as_tibble(goAnno[,c(2,1,4)]), dtlink, 
                         by = c("SYMBOL"="t_gn_sym"))
goAnno_drug <- as.data.frame(goAnno_drug)[,c("GOALL","ONTOLOGYALL","drug_name")]
goAnno_drug <- na.omit(goAnno_drug)
goAnno_drug <- goAnno_drug[!duplicated(goAnno_drug[,c("GOALL","drug_name")]),]
saveRDS(goAnno_drug, "./data/goAnno_drug.rds") # 7,517,403 X 3

######################
## make GO_DATA.rds ##
######################

# It is an intermediate file storing the GO annotation environment (GO term to 
# gene sets, gene to GO terms, GO term ID to description, GO ID to ontology) 
# used for `tsea` methods to accelerate the speed by avoiding building this 
# environment every time running the tsea functions 

GO_DATA <- clusterProfiler:::get_GO_data(
  OrgDb="org.Hs.eg.db", ont="ALL", keytype="SYMBOL")
saveRDS(GO_DATA, "./data/GO_DATA.rds")

###########################
## make GO_DATA_drug.rds ##
###########################

# It is an intermediate file storing the GO annotation environment (GO term to 
# drug sets, drug to GO terms, GO term ID to description, GO ID to ontology) 
# used for `dsea` methods to accelerate the speed by avoiding building this 
# environment every time running the dsea functions 

## Read in `goAnno_drug.rds` generated above
goAnno_drug <- readRDS("data/goAnno_drug.rds")
GO2GENE <- goAnno_drug[,c("GOALL","drug_name")]
GO_DATA_drug <- clusterProfiler:::build_Anno(GO2GENE, 
                                  clusterProfiler:::get_GO2TERM_table())
saveRDS(GO_DATA_drug, "./data/GO_DATA_drug.rds")

#########################
## make taurefList.rds ##
#########################

# It uses all signatures in the reference database (such as LINCS) to query 
# against itself as Qref to compute tau score of `gess_lincs` 
# method in `signatureSearch` package. Tau score compares observed enrichment 
# score to all others in Qref. It represents the percentage of reference queries
# with a lower |NCS| than |NCSq,r|, adjusted to retain the sign of NCSq,r. 
# NCSq,r is the normalized connectivity score for signature r relative to 
# query q. A tau of 90 indicates that only 10 percent of reference perturbations
# showed stronger connectivity to the query. For more details, please refer to 
# Subramanian et al., 2017, Cell, A Next Generation Connectivity Map: L1000 
# Platform and the First 1,000,000 Profiles

## Create Query Reference DB for Tau Score Computation of `gess_lincs` method
### Load `lincs` database created above
library(HDF5Array); library(SummarizedExperiment)
se = loadHDF5SummarizedExperiment("./data/lincs")
score_mat = assay(se)
### Create query list for all signatures in se
query_list <- lapply(colnames(score_mat), function(x) {
  vec <- as.matrix(score_mat[,x])
  names(vec) <- rownames(score_mat)
  sigvec = sort(vec, decreasing = TRUE)
  list(upset=utils::head(names(sigvec), 150), 
       downset=utils::tail(names(sigvec), 150))
})
names(query_list) = colData(se)$pert_cell_factor
###  Query signatures in the query_list against lincs database
###  To save time, the processing is paralleled with BiocParallel to run on 
###  CPU cores of a computer cluster with a scheduler (e.g. Slurm). 
#### Define submission function
f <- function(x, se, query_list, dest_dir) {
  require(signatureSearch)
  chunkno <- x 
  sz <- 10 # small enough to use short queue 
  qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz))
  myMA <- matrix(NA, length(query_list), sz, 
                 dimnames=list(names(query_list), seq_len(sz)))
  qlistone <- qlist[[chunkno]] 
  for(i in seq_along(qlistone)){
    qsig <- suppressMessages(qSig(qsig=list(upset=qlistone[[i]]$up, 
                                            downset=qlistone[[i]]$down), 
                 gess_method = "LINCS", refdb = se))
    lincs <- signatureSearch::gess_lincs(qsig, sortby=NA)
    resultDF <- result(lincs)
    ncs <- resultDF$NCS
    mynames <- paste(resultDF$pert, resultDF$cell, resultDF$type, sep="__")
    names(ncs) <- mynames
    myMA[,i] <- ncs[rownames(myMA)]
    colnames(myMA)[i] <- names(qlistone[i])
  }
  myMA <- myMA[, seq_along(qlistone)] 
  ## Only relevant for last entry that may not have as many columns as sz
  if(! dir.exists(dest_dir)) dir.create(dest_dir)
  write.table(as.data.frame(myMA), 
              file=paste0(dest_dir, "/result_", sprintf("%03d", chunkno)), 
              col.names=NA, quote=FALSE, sep="\t")
}
#### Split query into small chunks, each chunk will be processed on computer 
#### cluster as one process. Use `batchtools` to manage job submission.
#### Please make sure that `data` directory exists under your current working
#### directory of R session
library(batchtools)
sz <- 10 # small enough to finish each chunk within 2 hours
qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz))
if(! file.exists("slurm.tmpl")) 
  download.file("https://goo.gl/tLMddb", "slurm.tmpl", quiet=TRUE)
if(! file.exists(".batchtools.conf.R")) 
  download.file("https://goo.gl/5HrYkE", ".batchtools.conf.R", quiet=TRUE)
if(dir.exists("data/tau_queries_reg")) 
  unlink("data/tau_queries_reg", recursive=TRUE)
reg = makeRegistry(file.dir="data/tau_queries_reg", 
                   conf.file=".batchtools.conf.R")
dest_dir="data/tau_queries" 
## directory to store the intermediate results generated by function f
ids = batchMap(f, x = seq(along=qlist), 
            more.args = list(se=se, query_list=query_list, dest_dir=dest_dir))
#testJob(id = 1)
done <- submitJobs(ids, resources=list(walltime=36000, ncpus=1, memory=10240), 
                   partition="short", reg=reg)
#waitForJobs()
getStatus()

#### Inspect result and resubmit jobs for missing and empty files
fileDF <- file.info(list.files(dest_dir, pattern="result_*", full.names=TRUE))
index_empty <- as.numeric(gsub(".*_", "", row.names(fileDF[fileDF$size==0, ])))
qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz))
index_all_files <- seq_along(qlist)
index_exist <- as.numeric(gsub(".*_", "", row.names(fileDF)))
index_missing <- index_all_files[!index_all_files %in% index_exist]
index_repeat <- unique(sort(c(index_empty, index_missing)))
if(length(index_repeat)!=0){
  if(dir.exists("data/tau_queries_reg")) 
    unlink("data/tau_queries_reg", recursive=TRUE)
  reg = makeRegistry(file.dir="data/tau_queries_reg", 
                     conf.file=".batchtools.conf.R")
  ids = batchMap(f, x = index_repeat, more.args = 
                   list(se=se, query_list=query_list, dest_dir=dest_dir))
  done <- submitJobs(ids, resources=list(walltime=36000, ncpus=1, 
                                    memory=10240, partition="short"), reg=reg)
}
## Repeat the above steps until all jobs are successfully completed 
## (index_repeat is empty)

#### Organize results in list where each component contains data.frame with 
#### query results from one cell type
pathDF <- data.frame(query=names(query_list), 
                     path=rep(paste0("data/tau_queries/result_", 
                                     sprintf("%03d", seq_along(qlist))), 
                              sapply(qlist, length)))
pathDF <- data.frame(pathDF, target=gsub("^.*?__", "", pathDF$query))
pathList <- split(as.character(pathDF$path), factor(pathDF$target))
pathList <- sapply(pathList, unique) # eliminates duplicated imports of files
taurefList <- rep(NA, length(pathList)); names(taurefList) <- names(pathList)
taurefList <- as.list(taurefList) 
for(celltype in names(pathList)) {
  for(i in seq_along(pathList[[celltype]])) {
    tmpDF <- read.delim(pathList[[celltype]][i], row.names=1, check.names=FALSE)
    tmpDF <- round(tmpDF, 2) # Reduces data size
    colindex <- gsub("^.*?__", "", colnames(tmpDF)) %in% celltype
    tmpDF <- tmpDF[, colindex, drop=FALSE] 
    if(i==1) { 
      rowindex <- gsub("^.*?__", "", rownames(tmpDF)) %in% celltype
      containerDF <- tmpDF[rowindex, , drop=FALSE]
    } else {
      containerDF <- cbind(containerDF, tmpDF[rownames(containerDF),])
    }
    print(paste("Finished", i, "of", length(pathList[[celltype]]), 
                "results collected in", ncol(containerDF), "columns."))
  }
  taurefList[[celltype]] <- containerDF
  rm(containerDF); gc()
}
saveRDS(taurefList, file="data/tau_queries/taurefList.rds")

######################
## make ES_NULL.txt ##
######################

#  ES null distribution is generated with random queries for computing nominal 
#  P-values for Enrichment Score in gess_lincs result. 
#  To save time, the codes below parallel the process
#  with BiocParallel to run on CPU cores of a computer cluster with Slurm 
#  scheduler. It uses 1000 random queries searching against the lincs database
#  to generate the null distribution of ES.

## Create list of random queries
se = loadHDF5SummarizedExperiment("./data/lincs")
idnames <- rownames(se)
query_list <- signatureSearch:::randQuerySets(
  id_names=idnames, N_queries=1000, set_length=150)
## Define submission function
f <- function(x, query_list, se, dest_dir) {
  require(signatureSearch)
  chunkno <- x 
  sz <- 10 # small enough to run within 2 hours 
  qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz)) 
  qlistone <- qlist[[chunkno]] 
  for(i in seq_along(qlistone)){
    esout <- signatureSearch:::lincsEnrich(
      se, upset=qlistone[[i]]$up, downset=qlistone[[i]]$down, 
      sortby=NA, output="esonly", type=1)
    names(esout) <- colnames(se)
    wtcs = esout
    if(i==1) 
      myMA <- matrix(NA, length(esout), length(qlistone), 
                     dimnames=list(names(wtcs), seq_along(qlistone)))
    myMA[,i] <- wtcs[rownames(myMA)]
    colnames(myMA)[i] <- names(qlistone[i])
  }
  write.table(as.data.frame(myMA), 
              file=paste0(dest_dir, "/result_", sprintf("%03d", chunkno)), 
              col.names = NA, quote=FALSE, sep="\t")
}
## Split query into chunks, each chunk will be processed on cluster as 
## one process
sz <- 10 # small enough to finish each chunk within 2 hours
qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz)) 
dest_dir="data/es_queries" 
if(! file.exists("slurm.tmpl")) 
  download.file("https://goo.gl/tLMddb", "slurm.tmpl", quiet=TRUE)
if(! file.exists(".batchtools.conf.R")) 
  download.file("https://goo.gl/5HrYkE", ".batchtools.conf.R", quiet=TRUE)
if(dir.exists("data/es_queries_reg")) 
  unlink("data/es_queries_reg", recursive=TRUE)
reg = makeRegistry(file.dir="data/es_queries_reg", 
                   conf.file=".batchtools.conf.R")
ids = batchMap(f, x = seq(along=qlist), more.args = 
                 list(query_list=query_list, se=se, dest_dir=dest_dir))
#testJob(id = 1)
done <- submitJobs(ids, resources=list(walltime=36000, ncpus=1, memory=10240, 
                                       partition="short"), reg=reg)
getStatus()

## Inspect result and resubmit jobs for missing and empty files
fileDF <- file.info(list.files(dest_dir, pattern="result_*", full.names=TRUE))
index_empty <- as.numeric(gsub(".*_", "", row.names(fileDF[fileDF$size==0, ])))
qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz))
index_all_files <- seq_along(qlist)
index_exist <- as.numeric(gsub(".*_", "", row.names(fileDF)))
index_missing <- index_all_files[!index_all_files %in% index_exist]
index_repeat <- unique(sort(c(index_empty, index_missing)))
if(length(index_repeat)!=0){
  if(dir.exists("data/es_queries_reg")) 
    unlink("data/es_queries_reg", recursive=TRUE)
  reg = makeRegistry(file.dir="data/es_queries_reg", 
                     conf.file=".batchtools.conf.R")
  ids = batchMap(f, x = index_repeat, more.args = 
                   list(se=se, query_list=query_list, dest_dir=dest_dir))
  done <- submitJobs(ids, resources=list(walltime=36000, ncpus=1, memory=10240,
                                         partition="short"), reg=reg)
}
## Repeat the above steps until all jobs are successfully completed 
## (index_repeat is empty)
## 
## If all jobs are successfully completed, collect results in frequency table 
## with 3 digit accuracy
esMA <- data.frame(WTCS=as.character(round(rev(seq(-1, 1, by=0.001)), 3)), 
                   Freq=0, stringsAsFactors=FALSE) 
myfiles <- list.files(dest_dir, "result_\\d{3,3}", full.names=TRUE)
for(i in myfiles) {
  df <- read.delim(i, row.names=1)
  freq <- table(round(as.numeric(as.matrix(df),3),3)) 
  # processes entire data.frame
  freq <- freq[as.character(esMA[,1])]
  freq[is.na(freq)] <- 0
  esMA[,"Freq"] <- as.numeric(esMA[,"Freq"]) + as.numeric(freq)
  print(paste("Processed", i))
}
write.table(esMA, file=paste0(dest_dir,"/ES_NULL.txt"), quote=FALSE, row.names=FALSE, sep="\t")

normalizeCel <- function(chiptype_list, rerun = FALSE){
  mydir <- getwd()
  if(rerun==TRUE){
    for(i in names(chiptype_list)){
      dir.create(paste0("./data/",i))
      setwd(paste0("./data/",i))
      mydata <- affy::ReadAffy(filenames=chiptype_list[[i]], 
                               celfile.path="../CEL")
      eset <- affy::rma(mydata)
      write.table(affy::exprs(eset), file="./rma_exprs.xls", 
                  quote=FALSE, sep="\t",col.names=NA)
      setwd("../../")
    }
  } else {
    print("To execute function, set 'rerun=TRUE'")
  }
  setwd(mydir)
}

list2link <- function(list, type){
  link <- data.frame()
  for(i in 1:length(list)){
    if(length(list[[i]])!=0){
      if(type=="dt"){
        tmp <- data.frame(drug_name=names(list[i]), 
                          t_gn_sym=as.character(list[[i]]), 
                          stringsAsFactors = FALSE)
      }
      if(type=="td"){
        tmp <- data.frame(t_gn_sym=names(list[i]), 
                          drug_name=as.character(list[[i]]), 
                          stringsAsFactors = FALSE)
      }
      link <- rbind(link, tmp)
    }
  }
  return(link)
}

getPertInfo <- function(pert_iname){
  pert_info = NULL
  col = c("canonical_smiles", "description", "inchi_key", "inchi_string", 
          "moa", "molecular_formula", "pert_id", "pert_iname", "pert_type", 
          "pubchem_cid", "target")
  for(iname in pert_iname){
    json <- system(paste0(
                 'curl -X GET --globoff --header "Accept: application/json" --header "user_key: 6e09440bcbfa67ebfac72b928d4b9366" "https://api.clue.io/api/perts?filter=%7B%22where%22%3A%7B%22pert_iname%22%3A%22',
                  iname, '%22%7D%7D"'),
                  intern = TRUE, ignore.stderr = TRUE)
    
    if(length(json)==0){
      message("Curl error happens on ", iname)
      next()
    }
    if(grepl("error",json) | json=="[]"){
      message("API error happens on ", iname)
      next()
    }
    require(jsonlite)
    jsonr <- fromJSON(json)
    # print(i)
    # i=i+1
    # rownames(jsonr) <- pert_id
    for(item in col){
      if(! item %in% colnames(jsonr)){
        new <- data.frame(NA)
        colnames(new) <- item
        jsonr = cbind(jsonr, new)
      }
    }
    pert_info <- rbind(pert_info, jsonr[,col])
  }
  return(pert_info)
}

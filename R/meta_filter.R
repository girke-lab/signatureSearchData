#' LINCS signature meta data filter
#' 
#' Filter LINCS level5 signatures at specific concentation and treatment time
#' 
#' @param meta tibble, read in from LINCS signature info meta data 
#' @param dose dose/concentration of compound treatment, need to match elements
#' in 'pert_idose' column of 'meta'
#' @param time compound treatment time, need to match elements in 'pert_itime' 
#' column of 'meta'
#' @return tibble
#' @importFrom magrittr %<>%
#' @importFrom dplyr filter
#' @importFrom dplyr bind_cols
#' @importFrom dplyr distinct
#' @export
sig_filter <- function(meta, dose, time="24 h"){
    meta %<>% filter(pert_type=="trt_cp" & pert_idose==dose & pert_itime==time) 
    meta %<>% 
        bind_cols(alt_id=paste(meta$pert_iname, meta$cell_id, sep="_")) %>%
        bind_cols(pert_cell_factor=paste(meta$pert_iname, meta$cell_id, 
                                         meta$pert_type, sep="__")) %>% 
        distinct(alt_id, .keep_all=TRUE) # eliminate technical duplicates
    return(meta)
}

#' LINCS instance meta data filter
#' 
#' Filter LINCS level3 instances at specific concentation and treatment time
#' 
#' @param meta tibble, read in from LINCS instance info meta data
#' @param dose dose/concentration of compound treatment, need to match elements
#' in 'pert_dose' column of 'meta'
#' @param dose_unit unit of dose of compound treatment, need to match elements
#' in 'pert_dose_unit' column of 'meta'
#' @param time compound treatment time, need to match elements in 'pert_time' 
#' column of 'meta'
#' @param time_unit unit of time, need to match elements in 'pert_time_unit' 
#' column of 'meta'
#' @return tibble
#' @export
inst_filter <- function(meta, dose=10, dose_unit="um", time=24, time_unit="h"){
    inst42 %<>% filter(pert_type=="trt_cp" & pert_dose==dose & pert_dose_unit== dose_unit
                       & pert_time==time & pert_time_unit==time_unit)
    inst42 %<>% bind_cols(alt_id=paste(inst42$pert_iname, inst42$cell_id, sep="_")) %>% 
        bind_cols(pert_cell_factor=paste(inst42$pert_iname, inst42$cell_id, inst42$pert_type, sep="__"))
    return(inst42)
}
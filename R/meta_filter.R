#' Filter LINCS Level 5 by Condition
#' 
#' Function to filter Level 5 data from LINCS by specific conditions, such as compound 
#' concentrations and treatment times.
#' @param meta tibble containing experimental conditions of LINCS data
#' @param pert_type perturbation type, 'trt_cp' refers to treatment ('trt') with compound ('cp'). 
#' Description of other perturbation types ('pert_type') can be found in the GEO CMap LINCS 
#' User Guide v2.1 URL: https://docs.google.com/document/d/1rbHBy3DKekFm9lZouRG-ZcfLmCsfkUKzGPxxjqxPlYw/edit#
#' @param dose concentration of compound used for treatment, needs to match elements
#' in 'pert_idose' column of 'meta'
#' @param time compound treatment time, needs to match elements in 'pert_itime' 
#' column of 'meta'
#' @return tibble
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr bind_cols
#' @importFrom dplyr distinct
#' @examples 
#' meta <- data.frame(pert_idose=c("2 um","4 um", "10 um"), 
#'                    pert_itime="24 h",
#'                    pert_type="trt_cp", pert_iname=c("p1","p2","p3"),
#'                    cell_id="MCF7")
#' sig_filter(meta, dose="2 um")
#' @export
sig_filter <- function(meta, pert_type="trt_cp", dose, time="24 h"){
    meta %<>% dplyr::filter(pert_type==pert_type & pert_idose==dose & pert_itime==time) 
    meta %<>% 
        bind_cols(alt_id=paste(meta$pert_iname, meta$cell_id, sep="_")) %>%
        bind_cols(pert_cell_factor=paste(meta$pert_iname, meta$cell_id, 
                                         meta$pert_type, sep="__")) %>% 
        distinct(alt_id, .keep_all=TRUE) # eliminate technical duplicates
    return(meta)
}

#' Filter LINCS Level 3 by Condition
#' 
#' Function to filter Level 3 data from LINCS by specific conditions, such as compound 
#' concentrations and treatment times.
#' @param meta tibble containing experimental conditions of LINCS data
#' @param pert_type perturbation type, 'trt_cp' refers to treatment ('trt') with 
#' compound ('cp'). Description of other perturbation types ('pert_type') can 
#' be found in the GEO CMap LINCS User Guide v2.1 
#' URL: https://docs.google.com/document/d/1rbHBy3DKekFm9lZouRG-ZcfLmCsfkUKzGPxxjqxPlYw/edit#
#' @param dose concentration of compound used for treatment, needs to match elements
#' in 'pert_dose' column of 'meta'
#' @param dose_unit unit of dose of compound treatment, needs to match elements
#' in 'pert_dose_unit' column of 'meta'
#' @param time compound treatment time, needs to match elements in 'pert_time' 
#' column of 'meta'
#' @param time_unit unit of time, needs to match elements in 'pert_time_unit' 
#' column of 'meta'
#' @return tibble
#' @examples 
#' meta <- data.frame(pert_dose=c(2,4,10), pert_dose_unit="um", 
#'                    pert_time=24, pert_time_unit="h",
#'                    pert_type="trt_cp", pert_iname=c("p1","p2","p3"),
#'                    cell_id="MCF7")
#' inst_filter(meta)
#' @export
inst_filter <- function(meta, pert_type="trt_cp", dose=10, dose_unit="um", time=24, time_unit="h"){
    meta %<>% dplyr::filter(pert_type==pert_type & pert_dose==dose & pert_dose_unit== dose_unit
                       & pert_time==time & pert_time_unit==time_unit)
    meta %<>% bind_cols(alt_id=paste(meta$pert_iname, meta$cell_id, sep="_")) %>% 
        bind_cols(pert_cell_factor=paste(meta$pert_iname, meta$cell_id, meta$pert_type, sep="__"))
    return(meta)
}

alt_id = pert_dose = pert_dose_unit = pert_idose = pert_itime =
pert_time = pert_time_unit = pert_type = NULL

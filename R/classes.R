################################################################################
# Copyright (c) 2018, 2019 Genome Research Ltd.
#
# Author: Howard Lightfoot <cancerrxgene@sanger.ac.uk>
#
# This file is part of gdscmatrixanalyser.
#
# gdscmatrixanalyser. is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
################################################################################
#' @import dplyr
#' @import purrr
#' @import tidyr
NULL

# ---- Basic classes -----

#' drugset class
#'
#' @param drugset_id numeric drugset identifier
#' @param layout dataframe with columns: position, lib1, lib1_dose, treatment1, lib2, lib2_dose, treatment2, cmatrix
#'
#' @return drugset
#' @export
drugset <- function(drugset_id, layout = NA){
  my_drugset <- list(drugset_id = drugset_id, layout = layout)
  attr(my_drugset, "class") <- "drugset"
  return(my_drugset)
}

#' cell_line class
#'
#' @param cell_id
#' @param cell_line_name
#' @param master_cell_id
#' @param cosmic_id
#' @param sidm
#' @param tissue
#' @param cancer_type
#'
#' @return cell_line
#' @export
cell_line <- function(cell_id,
                      cell_line_name,
                      master_cell_id,
                      cosmic_id,
                      sidm = NULL,
                      tissue = NULL,
                      cancer_type = NULL){
  my_cell_line <- list(
    cell_id = cell_id,
    cell_line_name = cell_line_name,
    master_cell_id = master_cell_id,
    cosmic_id = cosmic_id,
    sidm = NULL,
    tissue = NULL,
    cancer_type = NULL
  )
  attr(my_cell_line, "class") <- "cell_line"
  return(my_cell_line)
}

#' plate class
#'
#' @param barcode
#' @param plate_date
#' @param scan_id
#' @param drugset_id
#' @param cell_line
#' @param intensities
#' @param mean_B
#' @param mean_nc1
#' @param mean_nc0
#' @param cmatrix_list
#'
#' @return
#' @export
#'
#' @examples
plate <- function(barcode,
                  plate_date = NULL,
                  scan_id = NULL,
                  drugset_id = NULL,
                  # A cell line object
                  cell_line = list(),
                  intensities = NULL,
                  mean_B = NULL,
                  mean_nc1 = NULL,
                  mean_nc0 = NULL,
                  cmatrix_list = list()
                  ){
  # intensities is a data frame of [position, intensity]
  my_plate <- list(barcode = barcode,
                   plate_date = plate_date,
                   scan_id = scan_id,
                   drugset_id = drugset_id,
                   cell_line = cell_line,
                   intensities =  intensities,
                   mean_B =mean_B,
                   mean_nc1 = mean_nc1,
                   mean_nc0 = mean_nc0,
                   cmatrix_list = cmatrix_list)
  attr(my_plate, "class") <- "plate"
  return(my_plate)
}

#' cmatrix class
#'
#' @param cmatrix
#' @param lib1
#' @param lib2
#'
#' @return
#' @export
cmatrix <- function(cmatrix, lib1, lib2) {
  my_matrix <- list(cmatrix = cmatrix,
                    lib1 = dose_response(lib = lib1),
                    lib2 = dose_response(lib = lib2)
                    )
  attr(my_matrix, "class") <- "cmatrix"
  return(my_matrix)
}

#' dose_repsonse object
#'
#' an attribute of cmatrix ???
#'
#' @param lib
#' @param ID
#' @param name
#' @param nlme_stats
#' @param target
#' @param pathway
#' @param owner
#'
#' @return
#' @export
dose_response <- function(lib = lib,
                          ID = NULL,
                          name = NULL,
                          nlme_stats = NULL,
                          target = NULL,
                          pathway = NULL,
                          owner = NULL){
  my_dose_response <- list(lib = lib,
                           ID = ID,
                           name = name,
                           nlme_stats = nlme_stats,
                           target = target,
                           pathway = pathway,
                           owner = owner)
  attr(my_dose_response, "class") <- "dose_response"
  return(my_dose_response)
}

################################################################################

#' new_drugset
#' Instantiate a new drugset obect
#'
#' @param screen_data from screen_data <- gdscDbR::getScreenData
#'
#' @return drugset
#' @export
new_drugset <- function(screen_data){
  ds_id <- unique(screen_data$DRUGSET_ID)
  try(if(length(ds_id) != 1) stop("More than one drugset id in data. Try screen_data %>% split(.$DRUGSET_ID) %>% map(~ new_drugset)"))

  ds_layout <- screen_data %>% distinct(POSITION, TAG, DRUG_ID, CONC) %>%
    mutate(CONC = round(CONC, digits = 8)) %>%
    mutate(TAG = ifelse(TAG == "NC-1", "NC1", TAG)) %>%
    mutate(TAG = ifelse(TAG == "NC-0", "NC0", TAG)) %>%
    mutate(TAG = ifelse(TAG == "UN-USED", "UNUSED", TAG)) %>%
    filter(TAG != "DMSO") %>%
    arrange(POSITION, TAG) %>%
    separate(TAG, into = c("lib", "dose", "treatment"), sep = "-", extra = "merge", fill = "right") %>%
    mutate(treatment = ifelse(dose %in% c("C", "S"), dose, treatment)) %>%
    mutate(dose = ifelse(dose %in% c("C", "S"), NA, dose)) %>%
    unite(col = "treatment", lib, dose, treatment, DRUG_ID, CONC, sep = "_") %>%
    group_by(POSITION) %>%
    mutate(tment = 1:n()) %>%
    ungroup() %>%
    mutate(tment = paste("t", tment, sep = "")) %>%
    spread(key = tment, value = treatment) %>%
    separate(t1, into = c("lib1", "lib1_dose", "treatment1", "lib1_drug_id", "lib1_conc"), sep = "_") %>%
    separate(t2, into = c("lib2", "lib2_dose", "treatment2", "lib2_drug_id", "lib2_conc"), sep = "_") %>%
    mutate(lib1_conc = suppressWarnings(as.numeric(lib1_conc)),
           lib2_conc = suppressWarnings(as.numeric(lib2_conc))) %>%
    # TEST THAT treatment1 == treatment2 unless NA...tbc - should be redundant
    mutate(treatment = treatment1) %>% select(-treatment1, -treatment2)

   if (ds_layout %>% filter(grepl("L\\d+", lib1), grepl("L\\d+", lib2)) %>% nrow() > 0){
       ds_layout <- ds_layout %>%
         left_join({.} %>%
                filter(grepl("^L\\d+", lib1) &
                         grepl("^L\\d+", lib2)) %>%
                distinct(lib1, lib2) %>%
                mutate(cmatrix = 1:n()), by = c("lib1", "lib2"))
   } else {
     ds_layout <- ds_layout %>% mutate(cmatrix = NA)
   }

  my_drugset <- drugset(drugset_id = ds_id, layout = ds_layout)

  return(my_drugset)
}

# ---- Calculate mean controls for a plate ----


#' set_plate_controls
#'
#' @param my_plate
#'
#' @return my_plate plate object
#' @export
set_plate_controls <- function(my_plate){
  my_drugset <- get_drugset(my_plate)
  plate_map <- inner_join(my_plate$intensities, my_drugset$layout, by = "POSITION")
  stopifnot(nrow(plate_map) <= 1536)
  my_plate$mean_B <- plate_map %>% filter(lib1 == "B") %>%
    summarise(mean(INTENSITY)) %>% unlist()
  my_plate$mean_nc1 <- plate_map %>% filter(lib1 == "NC1") %>%
    summarise(mean(INTENSITY)) %>% unlist()
  my_plate$mean_nc0 <- plate_map %>% filter(lib1 == "NC0") %>%
    summarise(mean(INTENSITY)) %>% unlist()
  return(my_plate)
}

#' set_plate_matrix_list
#'
#' sets the cmatrix_list for a plate object based on the plates drugset attribute
#'
#' @param my_plate
#'
#' @return my_plate
#' @export
set_plate_matrix_list <- function(my_plate) {
  drugset <- get_drugset(my_plate)
  cmatrices <- drugset$layout %>%
    filter(!is.na(cmatrix)) %>%
    select(cmatrix, lib1, lib2) %>%
    distinct()
  my_plate$cmatrix_list <- cmatrices %>% pmap(cmatrix)
  names(my_plate$cmatrix_list) <- paste("m", cmatrices$cmatrix, sep = "")
  return(my_plate)
}

#' set_lib_nlme_stats
#'
#' Set the combination matrix dose response nlme stats for each library
#'
#' @param my_plate
#' @param nlme_stats
#'
#' @return my_plate
#'
#' @details The data frame is stored in a plate e.g. my_plate$cmatrix_list$m1$lib2$nlme_stats
#' and has columns: CL, POSITION, drug, xmid, scal, x, conc = x_micromol, IC50, auc, RMSE, maxc, y, yhat
#'
#' @export
set_lib_nlme_stats <- function(my_plate, nlme_stats){
  stopifnot(!purrr::is_empty(my_plate$cmatrix_list))
  for (lib in c("lib1", "lib2")){
    lib_responses <- my_plate$cmatrix_list %>%
      map( ~ nlme_stats %>%
             filter(BARCODE == my_plate$barcode, lib_drug == .x[[lib]]$lib) %>%
             select(CL, POSITION, drug, xmid, scal, x, conc = x_micromol, IC50, auc, RMSE, maxc, y, yhat) %>%
             rename_with(function(x) paste(lib, x, sep = "_"),
                         .cols = c(drug, xmid, scal, x, conc, IC50, auc, RMSE, maxc, y, yhat))
      )

    my_plate$cmatrix_list <- map2(my_plate$cmatrix_list, lib_responses, function(.x, .y){
      .x[[lib]]$nlme_stats <- .y
      return(.x)
    })
  }
  # Check neither is empty because of failed wells and remove if so
  for (cm in names(my_plate$cmatrix)){
    if(nrow(my_plate$cmatrix_list[[cm]]$lib1$nlme_stats) == 0 |
       nrow(my_plate$cmatrix_list[[cm]]$lib2$nlme_stats) == 0) {
      my_plate$cmatrix_list[[cm]] <- NULL
    }
  }

  return(my_plate)
}

# check_lib_nlme_stats <- function(cmatr){}

#' get_monotherapies
#'
#' Make gdscIC50 nlme input data frame for single agent fitting
#'
#' @param my_plate
#' @param project_name
#'
#' @details Curves are fitted for every unique combination of MASTER_CELL_ID (CL)
#' and BARCODE and lib_drug (drug treatment as specified by the drug set).
#'
#' @return data frame of monotherapies ready for nlme fitting
#' @export
get_monotherapies <- function(my_plate, project_name){
  my_drugset <- get_drugset(my_plate)
  plate_map <- inner_join(my_plate$intensities,
                          my_drugset$layout %>% filter(treatment == 'S', grepl("^L\\d+", lib1)) %>% select(-starts_with("lib2_")),
                          by = "POSITION")
  stopifnot(nrow(plate_map) <= 1536)
  monotherapies <- plate_map %>%
    mutate(normalized_intensity = calc_viability(INTENSITY, my_plate$mean_B, my_plate$mean_nc1, trim = TRUE),
           norm_neg_pos = "NC1+B",
           BARCODE = my_plate$barcode,
           SCAN_ID = my_plate$scan_id,
           DRUGSET_ID = my_plate$drugset_id,
           CELL_ID = my_plate$cell_line$cell_id,
           CELL_LINE_NAME = my_plate$cell_line$cell_line_name,
           MASTER_CELL_ID = my_plate$cell_line$master_cell_id,
           time_stamp = gdscIC50::getTimeStamp(),
           sw_version = "gdscMatrixanalyser") %>%
    select(BARCODE,
           SCAN_ID,
           DRUGSET_ID,
           CELL_ID,
           CELL_LINE_NAME,
           MASTER_CELL_ID,
           POSITION,
           lib_drug = lib1,
           lib1_conc,
           treatment,
           normalized_intensity,
           norm_neg_pos,
           time_stamp,
           sw_version) %>%
    mutate(RESEARCH_PROJECT = project_name) %>%
    gdscIC50::setConcsForNlme(conc_col = "lib1_conc") %>%
    gdscIC50::prepNlmeData(cl_id = "MASTER_CELL_ID", drug_specifiers = c("BARCODE", "lib_drug")) %>%
    left_join(plate_map %>%
                distinct(lib_drug = lib1, DRUG_ID_lib = lib1_drug_id),
              by = "lib_drug")
  return(monotherapies)
}


#' set_drug_annotation
#'
#' Set drug annotation on plate cmatrix libraries
#'
#' @param my_plate
#' @param drug_annotation_df
#'
#' @return my_plate
#' @export
set_drug_annotation <- function(my_plate, drug_annotation_df){
  stopifnot(!purrr::is_empty(my_plate$cmatrix_list))
  for (cm in names(my_plate$cmatrix_list)){
    for (libnum in c("lib1", "lib2")){
      da <- drug_annotation_df %>% filter(DRUGSET_ID == my_plate$drugset_id,
                                          lib == my_plate$cmatrix_list[[cm]][[libnum]]$lib)
      stopifnot(nrow(da) == 1)
      my_plate$cmatrix_list[[cm]][[libnum]]$ID <- da$DRUG_ID
      my_plate$cmatrix_list[[cm]][[libnum]]$name <- da$DRUG_NAME
      my_plate$cmatrix_list[[cm]][[libnum]]$target <- da$PUTATIVE_TARGET
      my_plate$cmatrix_list[[cm]][[libnum]]$pathway <- da$PATHWAY_NAME
      my_plate$cmatrix_list[[cm]][[libnum]]$owner <- da$OWNED_BY
    }
  }
  return(my_plate)
}

# ---- Add cell line annotation ----
#' set_cell_line_annotation
#'
#' Add cell line annotation to plate$cell_line
#'
#' @param my_plate
#' @param cell_line_annotation_df
#'
#' @return my_plate
#' @export
set_cell_line_annotation <- function(my_plate, cell_line_annotation_df){
  ca <- cell_line_annotation_df %>% filter(MASTER_CELL_ID == my_plate$cell_line$master_cell_id)
  stopifnot(nrow(ca) == 1)
  my_plate$cell_line$sidm <- ca$SIDM
  my_plate$cell_line$tissue <- ca$TISSUE
  my_plate$cell_line$cancer_type <- ca$CANCER_TYPE
  return(my_plate)
}


# ---- Helper functions - not exported ----

calc_viability <- function(intensity, mean_B, mean_NC1, trim = TRUE){
  viability <- (intensity - mean_B) / (mean_NC1 - mean_B)
  if(trim){
    viability = ifelse(viability < 0, 0, viability)
    viability = ifelse(viability > 1, 1, viability)
  }
  return(viability)
}


# ---- Get drugset for plate ----
get_drugset <- function(my_plate){
  if(!exists("drugset_list")) {
    stop("You must have a list with name 'drugset_list' to proceed - try make_drugset_list(screen_data).")
    }
  ds <- paste("ds", my_plate$drugset_id, sep = "")
  my_drugset <- drugset_list[[ds]]
  stopifnot(my_plate$drugset_id == my_drugset$drugset_id)
  return(my_drugset)
}

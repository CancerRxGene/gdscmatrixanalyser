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


# --- Functions here create or operate on a list of plate objects using purrr style map functions ----

#' make_drugset_list
#'
#' Make list of drugsets from screen_data
#'
#' Many other functions will look up drugset details from the list using my_plate$drugset_id
#'
#' @param screen_data
#'
#' @return A list of drugset objects
#' @export
make_drugset_list <- function(screen_data){
  drugset_list <- screen_data %>% split(.$DRUGSET_ID) %>% map(~ new_drugset(.x))
  names(drugset_list) <- paste("ds", names(drugset_list), sep ="")
  return(drugset_list)
}

#' make_plate_list
#'
#' Make a list of plate objects from screen_data
#'
#' @param screen_data
#'
#' @return a list of plate objects
#' @export
make_plate_list <- function(screen_data) {
  screen_data %>% split(.$BARCODE) %>%
    map(~ plate(barcode = unique(.x$BARCODE),
                plate_date = lubridate::ymd(unique(.x$DATE_CREATED)),
                scan_id = unique(.x$SCAN_ID),
                drugset_id = unique(.x$DRUGSET_ID),
                cell_line = cell_line(cell_id = unique(.x$CELL_ID),
                                    cell_line_name = unique(.x$CELL_LINE_NAME),
                                    master_cell_id = unique(.x$MASTER_CELL_ID),
                                    cosmic_id = unique(.x$COSMIC_ID)
                ),
                intensities = .x %>% distinct(POSITION, INTENSITY) %>% arrange(POSITION))
    )
}

# ---- Set all annotation of drugs and cell lines ----
#' set_all_annotation
#'
#' Set all drug and cell line annotation for a list of plate objects
#'
#' @param plate_list a list of plate object
#' @param drug_annotation_df data frame with columns: DRUGSET_ID, DRUG_ID,
#'  DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, OWNED_BY.
#' @param cell_line_annotation_df data frame with columns: MASTER_CELL_ID, SIDM, TISSUE, CANCER_TYPE
#'
#' @return a list of annotated plates
#' @export
set_all_annotation <- function(plate_list, drug_annotation_df, cell_line_annotation_df){
  plate_list <- plate_list %>%
    map(~ set_drug_annotation(.x, drug_annotation_df)) %>%
    map(~ set_cell_line_annotation(.x, cell_line_annotation_df))
  return(plate_list)
}

#' make_nlme_data
#'
#' Make a gdscIC50 nlme input data_frame from a list of plates
#'
#' @param plate_list
#' @param project_name
#'
#' @details This function maps \code_monotherapies to a list of plate objects.
#'
#' @return an nlme_data data frame suitable for gdscIC50 workflow.
#'
#' @seealso \code{\link{get_monotherapies}}
#' @export
make_nlme_data <- function(plate_list, project_name){
  nlme_data <- plate_list %>% map_dfr(~ get_monotherapies(.x, project_name))
  return(nlme_data)
}


#' make_df_well_data
#'
#' Make a single dataframe of all matrix well statistics from a list of plate
#' objects.
#'
#' @param plate_list
#' @param project_id
#' @param project_name
#'
#' @return a data frame of matrix well statistics
#' @export
make_df_well_data <- function(plate_list, project_id, project_name){
  all_matrix_well_stats <- plate_list %>%
    map_dfr(~ well_stats_to_df(my_plate = .,
                               project_id = project_id,
                               project_name = project_name))
  return(all_matrix_well_stats)
}

# ---- Make a df of all matrix summary data ----
#' make_df_matrix_summaries
#'
#' Make a single dataframe of all matrix level summary statistics from a list
#' of plate objects.
#'
#' @param plate_list
#' @param project_id
#' @param project_name
#'
#' @return a data frame of all matrix well statistics
#' @export
make_df_matrix_summaries <- function(plate_list, project_id, project_name){
  all_matrix_summaries <- plate_list %>%
    map_dfr(~ matrix_stats_to_df(my_plate = .,
                                 project_id = project_id,
                                 project_name = project_name)
            )
  return(all_matrix_summaries)
}




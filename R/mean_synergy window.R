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


#' get_window_positions
#'
#' Calculate cmatrix window positions for a drugset and window size.
#'
#' @param my_drugset
#' @param cmatrix_number
#' @param window_size
#'
#' @return data frame of window positions
get_window_positions <- function(my_drugset, cmatrix_number, window_size) {
  matrix_layout <- my_drugset$layout %>%
    filter(cmatrix == cmatrix_number) %>%
    select(POSITION, lib1_dose, lib2_dose) %>%
    mutate(dose1 = as.integer(stringr::str_extract(lib1_dose, "\\d+")),
              dose2 = as.integer(stringr::str_extract(lib2_dose, "\\d+")))

  matrix_size <- matrix_layout %>%
    summarise(maxd1 = max(dose1), maxd2 = max(dose2))

  # Window start positions
  window_doses <- tibble(dose1 = 1:(matrix_size$maxd1- window_size + 1),
                             dose2 = 1:(matrix_size$maxd2- window_size + 1)) %>%
    cross_df() %>%
    pmap(.f = function(dose1, dose2){
      tibble(dose1 = dose1:(dose1 + (window_size -1)),
                 dose2 = dose2:(dose2 + (window_size -1))
      ) %>% cross_df()
    })

  window_positions <- window_doses %>%
    map(~ inner_join(matrix_layout, ., by = c("dose1", "dose2")))

  names(window_positions) <- paste("w", 1:length(window_positions), sep = "")

  window_positions <- map2(window_positions,
                           names(window_positions),
                           ~ .x %>% mutate(window = .y)) %>%
    map_dfr(~.) %>%
    group_by(window) %>%
    mutate(min_dose1 = min(dose1), min_dose2 = min(dose2)) %>%
    mutate(min_dose1 = paste0("D", min_dose1), min_dose2 = paste0("D", min_dose2)) %>%
    select(POSITION, window, lib1_dose, lib2_dose, min_dose1, min_dose2) %>%
    mutate(cmatrix = paste("m", cmatrix_number, sep = ""))

  return(window_positions)
}


# ---- Set up the windows for each drugset ----
get_drugset_cmatrices <- function(my_drugset){
  ds_cmatrices <- my_drugset$layout %>%
    filter(!is.na(cmatrix)) %>%
    distinct(cmatrix) %>%
    rename(m = cmatrix) %>%
    unlist()
  return(ds_cmatrices)
}

#' get_drugset_windows
#'
#' Get all overlapping window positions for all combination matrices in a drugset
#'  for a given window size (default of 3).
#'
#' @param my_drugset
#' @param window_size numeric , e.g., window_size  = 3 results in 3 by 3 windows
#'
#' @return data frame of window positions with columns: POSITION, window,
#'  lib1_dose, lib2_dose, min_dose1m min_dose2, cmatrix.
#'
#'  @seealso \code{\link{drugset}}
#'
#' @export
get_drugset_windows <- function(my_drugset, window_size = 3) {
  get_drugset_cmatrices(my_drugset) %>%
    map_dfr(
      ~ get_window_positions(my_drugset = my_drugset,
                             window_size = window_size,
                             cmatrix_number = .x)
    )
}

get_all_plate_effects <- function(my_plate){
  all_plate_effects_df <- my_plate$cmatrix_list %>%
    map_dfr(~ .x$effects %>% mutate(cmatrix = paste("m", .x$cmatrix, sep = "")))
  return(all_plate_effects_df)
}

#' synergy_window
#'
#' Find highest scoring synergy windows in each cmatrix in a plate
#'
#' @param my_plate
#' @param synergy_metric e.g. HSA, Bliss
#' @param synergy_only
#' @param window_size
#' @param drugset_windows
#'
#' @return data frame with columns: barcode, cmatrix, output_window_size,
#'  output_var, output_min_d1, output_min_d2
#'
#' @export
synergy_window <- function(my_plate,
                           synergy_metric = HSA,
                           synergy_only = FALSE,
                           window_size = 3,
                           drugset_windows){
  synergy_metric <- enquo(synergy_metric)

  if (missing(drugset_windows)){
    my_drugset <- get_drugset(my_plate)
    drugset_windows <- get_drugset_windows(my_drugset, window_size = window_size)
  }

  if (synergy_only){
    var_suffix <- "window_SO"
  } else {
    var_suffix <- "window"
  }

  input_var <- paste(quo_name(synergy_metric), "excess", sep = "_")
  input_var <- sym(input_var)
  output_var <- paste(quo_name(synergy_metric), var_suffix, sep = "_")
  output_window_size <- paste(quo_name(synergy_metric), var_suffix, "size", sep = "_")
  output_min_d1 <-  paste(quo_name(synergy_metric), var_suffix, "dose1", sep = "_")
  output_min_d2 <-  paste(quo_name(synergy_metric), var_suffix, "dose2", sep = "_")

  window_data <- inner_join(get_all_plate_effects(my_plate),
                            drugset_windows,
                            by = c("POSITION", "cmatrix"))

  if (synergy_only) {
    window_data <- window_data %>%
      mutate(!!input_var := ifelse(!!input_var < 0, NA, !!input_var))
  }

  mean_synergy_windows <- window_data %>%
    group_by(cmatrix, window) %>%
    summarise(mean_synergy_metric = mean(!!input_var, na.rm = TRUE))

  highest_synergy_windows <- mean_synergy_windows %>%
    group_by(cmatrix) %>%
    summarise(!!output_var := ifelse(
      all(is.na(mean_synergy_metric)),
      max(mean_synergy_metric, na.rm = FALSE),
      max(mean_synergy_metric, na.rm = TRUE))
    ) %>%
    left_join(mean_synergy_windows, by = "cmatrix") %>%
    filter(!!sym(output_var) == mean_synergy_metric) %>%
    select(-mean_synergy_metric) %>%
    inner_join(drugset_windows %>%
                 distinct(cmatrix, window, min_dose1, min_dose2),
               by = c("cmatrix", "window")
               ) %>%
    mutate(!!output_window_size := window_size,
           barcode = my_plate$barcode) %>%
    select(barcode,
           cmatrix,
           !!output_window_size,
           !!output_var,
           !!output_min_d1 := min_dose1,
           !!output_min_d2 := min_dose2
           ) %>%
    split(.$cmatrix)

  highest_synergy_windows <- highest_synergy_windows %>%
    map(function(.x){
      if(nrow(.x) > 1 ){
        .x <- .x %>% slice(1)
      }
      return(.x)
    })

  # In case of null results i.e. all windows are NA - synergy only
  for (cm in names(my_plate$cmatrix_list)) {
    if (is.null(highest_synergy_windows[[cm]])){
      highest_synergy_windows[[cm]] <- tibble(
        barcode = my_plate$barcode,
        cmatrix = cm,
        !!output_window_size := window_size,
        !!output_var := NA,
        !!output_min_d1 := NA,
        !!output_min_d2 := NA
      )
    }
  }

  return(highest_synergy_windows)
}



#' set_plate_list_window_metrics
#'
#' For a list of plate objects set the window values for a given synergy metric
#'
#' @param plate_list
#' @param synergy_windows_list
#' @param synergy_metric
#' @param synergy_only
#'
#' @return a list of plates
#' @export
set_plate_list_window_metrics <- function(plate_list,
                                          synergy_windows_list,
                                          synergy_metric,
                                          synergy_only = FALSE){
  synergy_metric <- enquo(synergy_metric)
  output_var <- paste(quo_name(synergy_metric), "window", sep = "_")
  if (synergy_only) {
    output_var <- paste(output_var, "SO", sep = "_")
  }
  plate_list <- plate_list %>% map(function(.x) {
    # print(.x$barcode)
    for (cm in names(.x$cmatrix_list)){
      .x$cmatrix_list[[cm]][[output_var]] <-
        synergy_windows_list[[.x$barcode]][[cm]] %>% select(-cmatrix, -barcode)
    }
    return(.x)
  })
  return(plate_list)
}








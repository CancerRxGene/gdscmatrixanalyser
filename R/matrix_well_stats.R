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

#' get_matrix_data
#'
#' Assemble data for matrix synergy stats from drugset details and plate intensities
#'
#' @param my_plate
#' @param cmatrix_name
#'
#' @return dataframe
#' @export
get_matrix_data <- function(my_plate, cmatrix_name = "m1"){
  my_drugset <- get_drugset(my_plate)
  my_plate$intensities %>%
    right_join(
      my_drugset$layout %>%
        filter(cmatrix == my_plate$cmatrix_list[[cmatrix_name]]$cmatrix),
      by = "POSITION")
}

#' get_matrix_and_lib_data
#'
#' @param my_plate
#' @param cmatrix_name
#'
#' @return
#' @export
get_matrix_and_lib_data <- function(my_plate, cmatrix_name = "m1"){
  matrix_data <- get_matrix_data(my_plate, cmatrix_name)
  lib1_doses <- sort(as.numeric(sub("D(\\d+)", "\\1", unique(matrix_data$lib1_dose))))
  lib2_doses <- sort(as.numeric(sub("D(\\d+)", "\\1", unique(matrix_data$lib2_dose))))

  matrix_data %>%
    left_join(my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats %>%
                arrange(desc(lib1_conc)) %>%
                mutate(dose = paste("D", lib1_doses, sep = "")) %>%
                select(-c(lib1_conc, POSITION, CL)),
              by = c("lib1_dose"="dose")) %>%
    left_join(my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats%>%
                arrange(desc(lib2_conc)) %>%
                mutate(dose = paste("D", lib2_doses, sep = "")) %>%
                select(-c(lib2_conc, POSITION, CL)),
              by = c("lib2_dose"="dose"))
}

#' get_matrix_and_synergy_data
#'
#' @param my_plate
#' @param cmatrix_name
#'
#' @return data frame
#' @export
get_matrix_and_synergy_data <- function(my_plate, cmatrix_name = "m1"){
  get_matrix_data(my_plate, cmatrix_name) %>%
    left_join(my_plate$cmatrix_list[[cmatrix_name]]$effects,
              by = "POSITION")
}

#' get_matrix_complete_data
#'
#' Create data frame of data for a cmatrix on a plate
#'
#' @param my_plate
#' @param cmatrix_name
#'
#' @return data frame
#' @export
get_matrix_complete_data <- function(my_plate, cmatrix_name = "m1"){
  matrix_data <- get_matrix_data(my_plate, cmatrix_name)
  lib1_doses <- sort(as.numeric(sub("D(\\d+)", "\\1", unique(matrix_data$lib1_dose))))
  lib2_doses <- sort(as.numeric(sub("D(\\d+)", "\\1", unique(matrix_data$lib2_dose))))

  matrix_data %>%
    left_join(my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats %>%
                arrange(desc(lib1_conc)) %>%
                mutate(dose = paste("D", lib1_doses, sep = "")) %>%
                select(-c(lib1_conc, POSITION, CL)),
              by = c("lib1_dose"="dose")) %>%
    left_join(my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats%>%
                arrange(desc(lib2_conc)) %>%
                mutate(dose = paste("D", lib2_doses, sep = "")) %>%
                select(-c(lib2_conc, POSITION, CL)),
              by = c("lib2_dose"="dose")) %>%
    left_join(my_plate$cmatrix_list[[cmatrix_name]]$effects,
              by = "POSITION")
}

#' set_inhibtion_effects
#'
#' Calculate cmatrix effects  for a plate as inhibition values
#'
#' @param my_plate a plate object
#'
#' @return plate object with a data frame of effects set per combination matrix
#'
#' @details the data frame is stored in, e.g., my_plate$cmatrix_list$m1$effects.
#' The dataframe has columns:
#' \itemize{
#'   \item{\code{POSITION} well position in the plate}
#'   \item{\code{inhibition} - measured from 0 (no inhibition of cell growth) to 1 (complete cell kill)
#' }
#'
#' @export
set_inhibition_effects <- function(my_plate){
  my_drugset <- get_drugset(my_plate)
  for (cm in names(my_plate$cmatrix_list)){
    my_plate$cmatrix_list[[cm]]$effects <- get_matrix_data(my_plate, cm) %>%
      mutate(inhibition = 1- calc_viability(INTENSITY, my_plate$mean_B, my_plate$mean_nc1)) %>%
      select(POSITION, inhibition)
  }
  return(my_plate)
}

#' set_combo_MaxE
#'
#' Set the combination maximum effect size for a each cmatrix on a plate
#' This is the 2nd highest inhibition value for a cmatrix
#'
#' @param my_plate a plate object
#'
#' @return a plate object with MaxE effects set e.g. my_plate$$cmatrix_list$m1$MaxE
#'
#' @export
set_combo_MaxE <- function(my_plate){
  my_drugset <- get_drugset(my_plate)
  for (cm in names(my_plate$cmatrix_list)){
      my_plate$cmatrix_list[[cm]]$MaxE <-
        sort(my_plate$cmatrix_list[[cm]]$effects$inhibition,
             decreasing = TRUE)[2]
  }
  return(my_plate)
}

# ... HSA, HSA_excess and HSA_mean_excess
#' set_HSA
#'
#' Highest Single Agent activity is calculated by comparing the combination
#' treatment effect to the most effective of the two fitted monotherapies.
#'
#' @param my_plate a plate object
#'
#' @return a plate object with added HSA effects information.
#'
#' @details The combination matrix effects data frame e.g., my_plate$cmatrix_list$m1$effects,
#'  which previously contains at least columns POSITION and inhibition is expanded to
#'  contain the columns:
#' \itemize{
#'   \item{\code{HSA} The maximum of the library 1 or library 2 monotherapy fitted dose
#'   response effects at the respective concentrations in the combinations matrix.}
#'   \item{\code{HSA_excess} The inhibition effect of the combination treatment minus the HSA effect.}
#' }
#'
#'  An additional data frame is also added, e.g., my_plate$$cmatrix_list$m1$HSA.
#'   The data frame has the columns:
#' \itemize{
#'   \item{\code{HSA_synergistic_wells} The number of wells displaying synergy, i.e., HSA_excess > 0}
#'   \item{\code{HSA_matrix} The mean HSA_excess for the combination matrix .
#'   \item{\code{HSA_matrix_SO} The mean HSA_excess for wells where HSA_excess > 0, i.e. synergy only wells.}
#' }
#'
#' @seealso \code{\link{set_inhibition_effects}}, \code{\link{set_Bliss}}
#'
#' @export
set_HSA <- function(my_plate){
  for (cm in names(my_plate$cmatrix_list)){
    hsa_df <-
      get_matrix_and_lib_data(my_plate, cm) %>%
      mutate(HSA = pmap_dbl(., .f = function(lib1_yhat, lib2_yhat, ...) max(lib1_yhat, lib2_yhat))) %>%
      select(POSITION, HSA)
    my_plate$cmatrix_list[[cm]]$effects <-
      my_plate$cmatrix_list[[cm]]$effects %>%
      left_join(hsa_df, by = "POSITION") %>%
      mutate(HSA_excess = inhibition - HSA)

    HSA_synergistic_wells <- my_plate$cmatrix_list[[cm]]$effects %>% filter(HSA_excess > 0) %>% nrow()
    HSA_mean_excess <- mean(my_plate$cmatrix_list[[cm]]$effects$HSA_excess)
    HSA_mean_excess_SO <- my_plate$cmatrix_list[[cm]]$effects %>%
      filter(HSA_excess > 0) %>%
      summarise(HSA_mean_SO = ifelse(length(HSA_excess) == 0,
                                     0,
                                     mean(HSA_excess, na.rm = TRUE))
                ) %>%
      unlist()

    my_plate$cmatrix_list[[cm]]$HSA <-
      tibble(HSA_synergistic_wells = HSA_synergistic_wells,
                 HSA_matrix = HSA_mean_excess,
                 HSA_matrix_SO = HSA_mean_excess_SO)
  }
  return(my_plate)
}

# ... Bliss_additivity
#' set_Bliss
#'
#' Bliss metrics are caclulated from the fitted monotherapies (yhat values) for
#' the two library drug treatments in the combinations matrix using the formula:
#'
#' \code{lib1_yhat + lib2_yhat - (lib1_yhat * lib2_yhat)}
#'
#' Set Bliss metrics for each cmatrix on a plate:
#'   the number of Bliss synergistic wells in a cmatrix,
#'   the mean excess over Bliss for a cmatrix,
#'   the mean excess over Bliss for only the synergistic wells.
#'
#' @param my_plate a plate object.
#'
#' @return a plate object with with added Bliss additivity metrics.
#'
#' @details The combination matrix effects data frame e.g., my_plate$cmatrix_list$m1$effects,
#'  which previously contained columns POSITION and inhibition is expanded to
#'  contain the columns:
#' \itemize{
#'   \item{\code{Bliss_additivity} Calculated using the monotherapy fitted dose
#'   response effects at the respective concentrations in the combinations matrix.}
#'   \item{\code{Bliss_excess} The inhibition effect of the combination treatment minus the Bliss additivity.}
#' }
#'
#'  An additional data frame is also added, e.g., my_plate$$cmatrix_list$m1$HSA.
#'   The data frame has the columns:
#' \itemize{
#'   \item{\code{Bliss_synergistic_wells} The number of wells displaying synergy, i.e., Bliss_excess > 0}
#'   \item{\code{Bliss_matrix} The mean Bliss_excess for the combination matrix .
#'   \item{\code{Bliss_matrix_SO} The mean Bliss_excess for wells where Bliss_excess > 0, i.e., synergy only wells.}
#' }
#'
#' @seealso \code{\link{set_inhibition_effects}}, \code{\link{set_HSA}}
#'
#' @return a plate object
#' @export
set_Bliss <- function(my_plate){
  for (cm in names(my_plate$cmatrix_list)){
    bliss_df <-
      get_matrix_and_lib_data(my_plate, cm) %>%
      mutate(Bliss_additivity = lib1_yhat + lib2_yhat - (lib1_yhat * lib2_yhat)) %>%
      select(POSITION, Bliss_additivity)

    my_plate$cmatrix_list[[cm]]$effects <-
      my_plate$cmatrix_list[[cm]]$effects %>%
      left_join(bliss_df, by = "POSITION") %>%
      mutate(Bliss_excess = inhibition - Bliss_additivity)

    Bliss_synergistic_wells <- my_plate$cmatrix_list[[cm]]$effects %>% filter(Bliss_excess > 0) %>% nrow()
    Bliss_mean_excess <- mean(my_plate$cmatrix_list[[cm]]$effects$Bliss_excess)
    Bliss_mean_excess_SO <- my_plate$cmatrix_list[[cm]]$effects %>%
      filter(Bliss_excess > 0) %>%
      summarise(mean_Bliss_SO = ifelse(length(Bliss_excess) == 0,
                                       0,
                                       mean(Bliss_excess, na.rm = TRUE))
      ) %>%
      unlist()

    my_plate$cmatrix_list[[cm]]$Bliss <-
      tibble(Bliss_synergistic_wells = Bliss_synergistic_wells,
                 Bliss_matrix = Bliss_mean_excess,
                 Bliss_matrix_SO = Bliss_mean_excess_SO)

  }
  return(my_plate)
}

#
# set_Loewe_index <- function(my_plate, my_drugset){
#   for (cm in names(my_plate$cmatrix_list)){
#     loewe_df <-
#       get_matrix_and_lib_data(my_plate, my_drugset, cm) %>%
#       mutate(lib1_equiv_dose = gdscIC50::getConcFromX(x = getX(y = combo_y,
#                                                      lib1_xmid,
#                                                      lib1_scal),
#                                             maxc = lib1_maxc)
#       ) %>%
#       mutate(lib2_equiv_dose = getConcFromX(x = getX(y = combo_y,
#                                                      lib2_xmid,
#                                                      lib2_scal),
#                                             maxc = lib2_maxc)
#       ) %>%
#       mutate(Loewe_index = (lib1_conc / lib1_equiv_dose) + (lib2_conc / lib2_equiv_dose)) %>%
#       select(POSITION, Loewe_index)
#
#     my_plate$cmatrix_list[[cm]]$effects <-
#       my_plate$cmatrix_list[[cm]]$effects %>%
#       left_join(loewe_df, by = "POSITION")
#   }
#   return(my_plate)
# }

#' mean_synergy_window
#'
#' Calculate mean synergy metrics per window for all possible windows in a
#' cmatrix
#'
#' @param my_plate
#' @param cmatrix_name
#' @param synergy_metric
#' @param window_size
#' @param synergy_only
#'
#' @return dataframe with columns: min_dose1, min_dose2, mean_synergy
#' @export
mean_synergy_window <- function(my_plate, cmatrix_name, synergy_metric = HSA_excess, window_size = 3, synergy_only = FALSE){
  # my_drugset <- get_drugset(my_plate)
  synergy_matrix_stats <- get_matrix_and_synergy_data(my_plate, cmatrix_name)
  synergy_metric <- enquo(synergy_metric)

  # stopifnot(!purrr::is_null(expr(`$`(synergy_matrix_stats, !!synergy_metric))))
  if (synergy_only){
    synergy_matrix_stats <- synergy_matrix_stats %>%
      mutate(!!synergy_metric := ifelse(!!synergy_metric < 0, NA, !!synergy_metric))
  }

  synergy_matrix_stats <- synergy_matrix_stats %>%
    mutate(dose1 = as.integer(stringr::str_extract(lib1_dose, "\\d+"))) %>%
    mutate(dose2 = as.integer(stringr::str_extract(lib2_dose, "\\d+")))
  # Return the max synergy for a window size e.g. 2x2 in the Block/matrix
  window_means <- map(
    .x = 1:(max(synergy_matrix_stats$dose1) - (window_size -1)),
    ~ .x %>% map2(.x = .x , .y = 1:(max(synergy_matrix_stats$dose2) - (window_size -1)),
                  .f = ~ synergy_matrix_stats %>%
                    filter(dose1 %in% (.x:(.x + window_size -1))) %>%
                    filter(dose2 %in% (.y:(.y + window_size -1))) %>%
                    select(dose1, dose2, HSA) %>%
                    # min dose as in lowest number but actually highest concentration
                    summarise(min_dose1 = min(dose1),
                              min_dose2 = min(dose2),
                              mean_synergy =  mean(HSA, na.rm = TRUE)) %>%
                    # If SO and all wells are non -synergistic i.e. NA - mean will return NaN
                    mutate(mean_synergy = ifelse(is.nan(mean_synergy), 0, mean_synergy)) %>%

                    mutate(min_dose1 = paste("D", min_dose1, sep = ""),
                           min_dose2 = paste("D", min_dose2, sep = ""))
    )) %>%
    flatten() %>%
    map_dfr(rbind)

  return(window_means)
}

#' mean_synergy_window_c
#'
#' A byte code compiled version of mean_synergy_window
#'
#' @return
#' @export
#'
#' @examples
mean_synergy_window_c <- function(){
  return(compiler::cmpfun(mean_synergy_window))
}


#' highest_mean_window
#'
#' The highest scoring synergy window
#' @param window_means
#'
#' @return the highest scoring window
highest_mean_window <- function(window_means){
  if (all(is.na(window_means$mean_synergy))){
    highest_scoring_window <-
      tibble(min_dose1 = NA,
                 min_dose2 = NA,
                 mean_synergy = NA,
                 window_size = unique(window_means$window_size)
      )
  }
  else{
    highest_scoring_window <- window_means %>%
      filter(mean_synergy == max({.}$mean_synergy, na.rm = TRUE))
    if (nrow(highest_scoring_window) > 1) {
      # Get the lowest doses if this ever happens
      highest_scoring_window <-
        highest_scoring_window %>%
        arrange(min_dose1, min_dose2) %>%
        slice(n())
    }
  }
  return(highest_scoring_window)
}


#' set_HSA_window_scores
#'
#' Set HSA window scores for a plate including synergy-only scores
#'
#' @param my_plate
#' @param window_size
#'
#' @return a plate object
#' @export
set_HSA_window_scores <- function(my_plate, window_size = 3){
  # my_drugset <- get_drugset(my_plate)
  for (cm in names(my_plate$cmatrix_list)){
      highest_scoring_window <-
        mean_synergy_window(my_plate,
                            cmatrix_name = cm,
                            synergy_metric = HSA_excess,
                            window_size = window_size) %>%
        highest_mean_window()

      my_plate$cmatrix_list[[cm]]$HSA_window <-
        highest_scoring_window %>%
        select(HSA_window_size = window_size,
               HSA_window = mean_synergy,
               HSA_window_dose1 = min_dose1,
               HSA_window_dose2 = min_dose2)

      #  set the synergy only windows at the same time
      highest_scoring_window <-
        mean_synergy_window(my_plate,
                            cmatrix_name = cm,
                            synergy_metric = HSA_excess,
                            window_size = window_size,
                            synergy_only = TRUE) %>%
        highest_mean_window()

      my_plate$cmatrix_list[[cm]]$HSA_window_SO <-
        highest_scoring_window %>%
        select(HSA_window_SO_size = window_size,
               HSA_window_SO = mean_synergy,
               HSA_window_SO_dose1 = min_dose1,
               HSA_window_SO_dose2 = min_dose2)

  }
  return(my_plate)
}

#' set_Bliss_window_scores
#'
#' Set Bliss window scores for a plate including synergy-only scores
#'
#' @param my_plate
#' @param window_size
#'
#' @return a plate project
#' @export
set_Bliss_window_scores <- function(my_plate, window_size = 3){
  # my_drugset <- get_drugset(my_plate)
  for (cm in names(my_plate$cmatrix_list)){
    highest_scoring_window <-
      mean_synergy_window(my_plate,
                          cmatrix_name = cm,
                          synergy_metric = Bliss_excess,
                          window_size = window_size) %>%
      highest_mean_window()

    my_plate$cmatrix_list[[cm]]$Bliss_window <-
      highest_scoring_window %>%
      select(Bliss_window_size = window_size,
             Bliss_window = mean_synergy,
             Bliss_window_dose1 = min_dose1,
             Bliss_window_dose2 = min_dose2)

    #  set the synergy only windows at the same time
    highest_scoring_window <-
      mean_synergy_window(my_plate,
                          cmatrix_name = cm,
                          synergy_metric = Bliss_excess,
                          window_size = window_size,
                          synergy_only = TRUE) %>%
      highest_mean_window()

    my_plate$cmatrix_list[[cm]]$Bliss_window_SO <-
      highest_scoring_window %>%
      select(Bliss_window_SO_size = window_size,
             Bliss_window_SO = mean_synergy,
             Bliss_window_SO_dose1 = min_dose1,
             Bliss_window_SO_dose2 = min_dose2)

  }
  return(my_plate)
}


#' set_plate_attributes_pre_nlme
#'
#' Set the plate attributes necessary to fit dose response curves using GDSC nlme
#' model.
#'
#' @param my_plate
#' @param day1_data
#'
#' @details The following attributes are set:
#' \itemize{
#'   \item{control well mean values: mean_nc0 (media only), mean_nc1 (DMSO) and mean_B (blank no cells)}
#'   \item{data summary for a matched day 1 plate \}
#'   \item{cmatrix_list - combination matrices with library ids from the corresponding drugset}
#' }
#'
#' @return a plate object
#' @seealso \code{\link{set_plate_controls}}, \code{\link{set_day1_data}}, \code{\link{set_plate_matrix}}
#'
#' @export
set_plate_attributes_pre_nlme <- function(my_plate, day1_data){
  # my_drugset <- get_drugset(my_plate)
  my_plate <- set_plate_controls(my_plate)
  my_plate <- set_day1_data(my_plate, day1_data)
  my_plate <- set_plate_matrix_list(my_plate)
  return(my_plate)
}

#' set_plate_attributes_post_nlme
#'
#' Set plate attributes: nlme stats, inhibtion effects, combination max effect,
#'   HSA metrics, Bliss metrics, HSAwindow metrics, Bliss window metrics
#'
#' @param my_plate
#' @param day1_data
#' @param nlme_stats
#'
#' @return a plate object
#' @export
set_plate_attributes_post_nlme <- function(my_plate, day1_data, nlme_stats){
  my_plate <- tryCatch(set_lib_nlme_stats(my_plate, nlme_stats),
    error = function(e) print(paste(my_plate$barcode, "set_lib_nlme_stats", e, sep = " "))
  )
  my_plate <- tryCatch(set_inhibition_effects(my_plate),
    error = function(e) print(paste(my_plate$barcode, "set_inhibition_effects", e, sep = " "))
  )
  my_plate <- tryCatch(set_combo_MaxE(my_plate),
    error = function(e) print(paste(my_plate$barcode, "set_combo_MaxE", e, sep = " "))
  )
  my_plate <- tryCatch(set_HSA(my_plate),
    error = function(e) print(paste(my_plate$barcode, "set_HSA", e, sep = " "))
  )
  my_plate <- tryCatch(set_Bliss(my_plate),
    error = function(e) print(paste(my_plate$barcode, "set_Bliss", e, sep = " "))
  )
  my_plate <- tryCatch(set_HSA_window_scores(my_plate, window_size = 3),
    error = function(e) print(paste(my_plate$barcode, "set_HSA_window_scores", e, sep = " "))
  )
  my_plate <- tryCatch(set_Bliss_window_scores(my_plate, window_size = 3),
    error = function(e) print(paste(my_plate$barcode, "set_Bliss_window_scores", e, sep = " "))
  )
  return(my_plate)
}

#' well_stats_to_df
#'
#' Create a dataframe of well stats to write to csv and then input to
#' matrixexplorer database.
#'
#' @param my_plate
#' @param project_id
#' @param project_name
#'
#' @return a data frame of all well statistics
#' @export
well_stats_to_df <- function(my_plate, project_id, project_name){
  my_drugset <- get_drugset(my_plate)
  well_stats_df <-
    tibble(PROJECT_ID = project_id,
               PROJECT_NAME = project_name,
               BARCODE = my_plate$barcode,
               DRUGSET_ID = my_plate$drugset_id,
               CELL_LINE_NAME = my_plate$cell_line$cell_line_name,
               MASTER_CELL_ID = my_plate$cell_line$master_cell_id,
               COSMIC_ID = my_plate$cell_lin$cosmic_id,
               CELL_ID = my_plate$cell_line$cell_id,
               TISSUE = my_plate$cell_line$tissue,
               CANCER_TYPE = my_plate$cell_line$cancer_type
    ) %>%
    cbind(
      names(my_plate$cmatrix_list) %>%
        map_dfr(~ get_matrix_complete_data(my_plate, .x))
    )
  return(well_stats_df)
}

#' matrix_stats_to_df
#'
#' Create a dataframe of matrix stats to write to csv and then input to
#' matrixexplorer database.
#'
#' @param my_plate
#' @param project_id
#' @param project_name
#'
#' @return data frame
#' @export
matrix_stats_to_df <- function(my_plate, project_id, project_name){
  my_drugset <- get_drugset(my_plate)
  matrix_stats_df <-
    tibble(PROJECT_ID = project_id,
               PROJECT_NAME = project_name,
               BARCODE = my_plate$barcode,
               DRUGSET_ID = my_plate$drugset_id,
               CELL_LINE_NAME = my_plate$cell_line$cell_line_name,
               MASTER_CELL_ID = my_plate$cell_line$master_cell_id,
               COSMIC_ID = my_plate$cell_lin$cosmic_id,
               CELL_ID = my_plate$cell_line$cell_id,
               TISSUE = my_plate$cell_line$tissue,
               CANCER_TYPE = my_plate$cell_line$cancer_type
    ) %>%
    cbind(
      names(my_plate$cmatrix_list) %>%
        map_dfr(~ get_matrix_summary(my_plate, .x))
    )
  return(matrix_stats_df)
}

#' get_matrix_summary
#'
#' Create a dataframe of matrix summary stats to write to csv and then input to
#' matrixexplorer database.
#'
#' @param my_plate
#' @param project_id
#' @param project_name
#'
#' @return data frame
#' @export
get_matrix_summary <- function(my_plate, cmatrix_name){
  matrix_summary <-
    tryCatch(
      tibble(
        cmatrix = my_plate$cmatrix_list[[cmatrix_name]]$cmatrix,
        well_treatments = paste(my_plate$cmatrix_list[[cmatrix_name]]$lib1$lib, my_plate$cmatrix_list[[cmatrix_name]]$lib2$lib, sep = "-"),
        lib1 = my_plate$cmatrix_list[[cmatrix_name]]$lib1$lib,
        lib1_ID =  my_plate$cmatrix_list[[cmatrix_name]]$lib1$ID,
        lib1_name = my_plate$cmatrix_list[[cmatrix_name]]$lib1$name,
        lib1_conc = my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats$lib1_conc %>% max(),
        lib1_RMSE = my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats$lib1_RMSE %>% unique(),
        lib1_MaxE = my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats$lib1_yhat %>% max(),
        lib1_IC50_ln = my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats$lib1_IC50 %>% unique(),
        lib1_IC50_uM = my_plate$cmatrix_list[[cmatrix_name]]$lib1$nlme_stats$lib1_IC50 %>% unique() %>% exp(),
        lib1_target =  my_plate$cmatrix_list[[cmatrix_name]]$lib1$target,
        lib1_pathway =  my_plate$cmatrix_list[[cmatrix_name]]$lib1$pathway,
        lib1_owner = my_plate$cmatrix_list[[cmatrix_name]]$lib1$owner,
        lib2 = my_plate$cmatrix_list[[cmatrix_name]]$lib2$lib,
        lib2_ID =  my_plate$cmatrix_list[[cmatrix_name]]$lib2$ID,
        lib2_name =  my_plate$cmatrix_list[[cmatrix_name]]$lib2$name,
        lib2_conc = my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats$lib2_conc %>% max(),
        lib2_RMSE = my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats$lib2_RMSE %>% unique(),
        lib2_MaxE = my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats$lib2_yhat %>% max(),
        lib2_IC50_ln = my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats$lib2_IC50 %>% unique(),
        lib2_IC50_uM = my_plate$cmatrix_list[[cmatrix_name]]$lib2$nlme_stats$lib2_IC50 %>% unique() %>% exp(),
        lib2_target = my_plate$cmatrix_list[[cmatrix_name]]$lib2$target,
        lib2_pathway = my_plate$cmatrix_list[[cmatrix_name]]$lib2$pathway,
        lib2_owner  = my_plate$cmatrix_list[[cmatrix_name]]$lib2$owner,
        matrix_size =  my_plate$cmatrix_list[[cmatrix_name]]$effects %>% nrow(),
        combo_MaxE =  my_plate$cmatrix_list[[cmatrix_name]]$MaxE
      ) %>%
        mutate(Delta_MaxE_lib1 = combo_MaxE - lib1_MaxE,
               Delta_MaxE_lib2 = combo_MaxE - lib2_MaxE)
      %>% mutate(Delta_combo_MaxE_day1 =
                          combo_MaxE - my_plate$day1$day1_inhibition_scale) %>%
        cbind(my_plate$cmatrix_list[[cmatrix_name]]$HSA) %>%
        cbind(my_plate$cmatrix_list[[cmatrix_name]]$HSA_window) %>%
        cbind(my_plate$cmatrix_list[[cmatrix_name]]$HSA_window_SO) %>%
        cbind(my_plate$cmatrix_list[[cmatrix_name]]$Bliss) %>%
        cbind(my_plate$cmatrix_list[[cmatrix_name]]$Bliss_window) %>%
        cbind(my_plate$cmatrix_list[[cmatrix_name]]$Bliss_window_SO) %>%
        cbind(my_plate$day1),
      error = function(e) print(my_plate$barcode)
    )
  return(matrix_summary)
}


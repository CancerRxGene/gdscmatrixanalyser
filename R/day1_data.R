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

#' Process day 1 data
#'
#' Summarise screen day 1 data with mean and standard deviation of NC-0 wells
#' per barcode
#'
#' @param day1_data
#'
#' @return day_1_data summarised data frame
#' @export
process_day1_data <- function(day1_data) {
  day1_data <- day1_data %>%
    mutate(plate_date = lubridate::ymd(DATE_CREATED)) %>%
    select(BARCODE, plate_date, MASTER_CELL_ID, CELL_ID, INTENSITY, TAG) %>%
    filter(TAG == "NC-0") %>%
    group_by(BARCODE, plate_date, MASTER_CELL_ID) %>%
    summarise(day1_intensity_mean = mean(INTENSITY), day1_intensity_sd = sd(INTENSITY))
  return(day1_data)
}

#' Set day 1 data
#'
#' Set the day 1 data for a plate object from a dataframe of all day 1 data for
#' the screen
#'
#' @param my_plate matrixanalyser plate object
#' @param day1_data Data frame of all screen day 1 data previously processed by
#'  process_day1_data
#'
#' @return my_plate
#' @export
set_day1_data <- function(my_plate, day1_data){
  d1 <- day1_data %>%
    filter(plate_date == my_plate$plate_date,
           MASTER_CELL_ID == my_plate$cell_line$master_cell_id)

  if(nrow(d1) > 1) {
    stop(paste("nrow(d1) > 1 is not TRUE. Barcode: ", my_plate$barcode, sep = ""))
  }
  if(nrow(d1) == 0){
    my_plate$day1 <- tibble(day1_intensity_mean = NA,
                                day1_viability_mean = NA,
                                day1_inhibition_scale = NA,
                                growth_rate = NA,
                                doubling_time = NA)
  }
  else{
    my_plate$day1 <-
      tibble(day1_intensity_mean = d1$day1_intensity_mean) %>%
      mutate(day1_viability_mean = (day1_intensity_mean - my_plate$mean_B) / (my_plate$mean_nc1 - my_plate$mean_B),
             day1_inhibition_scale = 1 - day1_viability_mean,
             growth_rate = log(my_plate$mean_nc1 / day1_intensity_mean, base = 2),
             doubling_time = 72 / growth_rate)
  }

  return(my_plate)
}


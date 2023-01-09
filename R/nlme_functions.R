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

# This function is in gdscIC50 but needs exporting
getX <- function(y, xmid, scal){
  x <- xmid - scal*(log((1 - y) / y))
  # x <- xmid - scal*(log((y) / (1 -y)))
  return(x)
}


#' calculate_Emax
#'
#' Calculate the effect at the maximum concentration
#'
#' @param nlme_stats data frame of nlme stats
#'
#' @return data frame of nlme stats with column \code{Emax} added
#'
#' @seealso \code{\link[gdscIC50]{calcNlmeStats}
#'
#' @export
calculate_Emax <- function(nlme_stats){
  emax_stats <- nlme_stats %>%
    filter(maxc == x_micromol) %>%
    mutate(Emax = yhat) %>%
    distinct(drug, CL, Emax)
  nlme_stats <- left_join(nlme_stats, emax_stats)
  return(nlme_stats)
}

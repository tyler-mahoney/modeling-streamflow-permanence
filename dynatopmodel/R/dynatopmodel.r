#' Implementation of the Dynamic TOPMODEL hydrological model.
#'
#' @description
#' A native R implementation and enhancement of Dynamic TOPMODEL, an extension to
#' the semi-distributed hydrological model TOPMODEL. It includes some digital
#' terrain analysis functions for discretisation of catchments by topographic
#' indexes and other geo-referenced layers containing relevant landscape data.
#'
#' TOPMODEL (Beven & Kirkby, 1979) is a well-established and widely used
#' hydrological model that implements a spatial aggregation strategy
#' ("discretisation") in order to reduce its computational demands. Hydrological
#' similar areas identified by the discretisation procedure are referred to as
#' hydrological response units (HRUs). Beven and Freer (2001) introduced a
#' "dynamic" variant that addressed some of the limitations of the original
#' TOPMODEL but which retained its computational and parametric efficiency. In
#' particular, the original assumption of a quasi-steady water table was replaced
#' by time-dependent kinematic routing between and within HRUs. This allows a
#' more flexible discretisation approach, whereby any type of landscape data can
#' be used to identify the HRUs.
#'
#' With the introduction of a single new parameter SDmax specifying the maximum
#' storage deficit before downslope flow out of a HRU ceases, variable upslope
#' drainage areas can be simulated. This allows application to arid catchments
#' subject to seasonal drying of upslope areas.
#'
#' The 2001 version was implemented in FORTRAN and it, and its source code,
#' have not been made generally available. It has been applied in a number of
#' studies (see Metcalfe et al, 2015). A modified version with chemical stores
#' attached to each HRU was implemented by Page et al. (2007) and applied to
#' modelling the Cl signal in the Hafren, Mid-Wales.

#' This new version, described in detail in Metcalfe et al. (2015), retains the
#' core dynamics of the FORTRAN implementation but makes use of data storage and
#' vectorisation features of the R language to allow efficient scaling of the
#' problem domain. This version was utilised by Metcalfe et al. (2017) to supply
#' hillslope runoff to a new hydraulic channel routing model for evaluation
#' of the effectiveness of arrays of in-channel barriers on mitigating flood risk.
#'
#' A new, semi-distributed surface routing algorithm has been introduced in this
#' version. This allows examination of surface storages as they move downslope
#' and specification of different effective velocities throughout each unit. It
#' also allows for a modified routing matrix that can reflect situations where
#' the surface flow pathways differ from the topography. These could be used to
#' simulate a unit associated with one or more surface features designed to
#' intercept storm runoff, such as excavated ponds or bunds (see Hankin et al.,
#' 2016, 2017; Metcalfe, 2017).
#'
#' The preprocessing routines supplied incorporate handling of geo-referenced
#' spatial data to allows it to integrate with modern GIS through
#' industry-standard file formats, such as GEOTiff and ESRI Shapefiles.
#'
#' @seealso \code{\link{discretise}}
#' @seealso \code{\link{run.dtm}}
#'
#' @references
#' Beven, K. J. and M. J. Kirkby (1979). A physically based variable contributing area model of basin hydrology. Hydrol. Sci. Bull 24(1): 43-69.
#'
#' Beven, K. J. and J. Freer (2001). A Dynamic TOPMODEL. Hydrological Processes 15(10): 1993-2011.
#'
#' Hankin, B., Craigen I., Chappell, N., Metcalfe, P., Page, T. (2016). The Rivers Trust Life-IP Natural Course Project: Strategic Investigation of Natural Flood Management in Cumbria.Technical Report. Available at http://naturalcourse.co.uk/uploads/2017/04/2016s4667-Rivers-Trust-Life-IP-NFM-Opportunities-Technical-Report-v8.0.pdf. Rivers Trust, Callington, Cornwall, UK.
#'
#' Hankin, B., Metcalfe, P., Johnson, D., Chappell, N., Page, T., Craigen, I., Lamb, R., Beven, K. (2017). Strategies for Testing the Impact of Natural Flood Risk Management Measures. In Hromadka, T. & Rao, P. (eds). Flood Risk Management. InTech, Czech Republic. ISBN 978-953-51-5526-3.
#'
#' Metcalfe, P. (2017). Development of a modelling framework for integrated flood risk management (Doctoral dissertation). Lancaster University, UK.
#'
#' Metcalfe, P., Beven, K., & Freer, J. (2015). Dynamic TOPMODEL: a new implementation in R and its sensitivity to time and space steps. Environmental Modelling & Software, 72, 155-172.
#'
#' Metcalfe, P., Beven, K., Hankin, B., & Lamb, R. (2017). A modelling framework for evaluation of the hydrological impacts of nature?based approaches to flood risk management, with application to in-channel interventions across a 29 km2 scale catchment in the United Kingdom. Hydrological Processes, 31(9), 1734-1748.
#'
#' Page, T., Beven, K. J., Freer, J., & Neal, C. (2007). Modelling the chloride signal at Plynlimon, Wales, using a modified dynamic TOPMODEL incorporating conservative chemical mixing (with uncertainty). Hydrological Processes, 21(3), 292-307.
#'
#' @docType package
#' @name dynatopmodel
#' @import deSolve
#' @import xts
#' @import zoo
#' @import rgdal
#' @import sp
#' @import raster
#' @import rgeos
#' @import methods
#' @importFrom grDevices colorRampPalette colours dev.cur dev.flush dev.hold dev.list dev.new dev.off dev.set jpeg rainbow
#' @importFrom graphics abline axis box grconvertY layout legend mtext par title
#' @importFrom stats end start time xtabs runif
#' @importFrom utils capture.output read.csv read.table setTxtProgressBar txtProgressBar write.table
#' @importFrom lubridate wday mday round_date hour yday
NULL

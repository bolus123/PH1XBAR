#' Thickness measurement of silicon wafer
#'
#' A dataset containing the thickness measurements in nm
#' at different positions on the silicon wafer
#'
#' @format A matrix with 30 rows and 5 variables:
#' \describe{
#'   \item{pos1}{Thickness measurement at Position 1 (outer circle)}
#'   \item{pos2}{Thickness measurement at Position 2 (outer circle)}
#'   \item{pos3}{Thickness measurement at Position 3 (middle circle)}
#'   \item{pos4}{Thickness measurement at Position 4 (middle circle)}
#'   \item{pos5}{Thickness measurement at Position 5 (inner circle)}
#' }
#' @references Roes, Kit CB, and Ronald JMM Does. "Shewhart-type charts in nonstandard situations." Technometrics 37.1 (1995): 15-24
"grinder_data"

#' Bore diameter in manufacturing automotive driver gears
#'
#' A dataset cotaining bore diameter measurements in mm
#'
#' @format A matrix with 20 rows and 5 variables:
#' \describe{
#'   \item{X1}{Diameter measurement at Position 1}
#'   \item{X2}{Diameter measurement at Position 2}
#'   \item{X3}{Diameter measurement at Position 3}
#'   \item{X4}{Diameter measurement at Position 4}
#'   \item{X5}{Diameter measurement at Position 5}
#' }
#' @references Wooluru, Yerriswamy, D. R. Swamy, and P. Nagesh. "THE PROCESS CAPABILITY ANALYSIS-A TOOL FOR PROCESS PERFORMANCE MEASURES AND METRICS-A CASE STUDY." International Journal for Quality Research 8.3 (2014).
"bore_diameter_data"

#' Prescription fentanyl consumption in Preston county, WV
#'
#' A dataset containing prescription fentanyl consumption in Preston county, WV, measured using MME percapita. This is a subset from Rich et al. <doi: 10.21105/joss.02450>
#'
#' @format A vector with 60 elements
#' @references Rich, S., Tran, A. B., Williams, A., Holt, J., Sauer, J., & Oshan, T. M. (2020). arcos and arcospy: R and Python packages for accessing the DEA ARCOS database from 2006-2014. Journal of Open Source Software, 5(53), 2450.
"preston_data"

#' Seasonal snowfall in inches in Minneapolis/St. Paul, MN
#'
#' A dataset containing snowfalls measured in inches in Minneapolis/St. Paul, MN.
#'
#' @format A data frame with 82 rows and 4 variables:
#' \describe{
#'   \item{Year}{year of the snowfalls}
#'   \item{jan}{snowfalls in January}
#'   \item{feb}{snowfalls in February}
#'   \item{mar}{snowfalls in March}
#'  }
#' @references Mukherjee, P. S. (2016). On phase II monitoring of the probability distributions of univariate continuous processes. Statistical Papers, 57(2), 539-562.
"snowfall_data"


#' Pistonring data
#'
#' A dataset containing piston ring data
#'
#' @format A data frame with 25 rows and 5 variables:
#' \describe{
#'   \item{X1}{Observation 1 in subgroups}
#'   \item{X2}{Observation 2 in subgroups}
#'   \item{X3}{Observation 3 in subgroups}
#'   \item{X4}{Observation 4 in subgroups}
#'   \item{X5}{Observation 5 in subgroups}
#'  }
#' @references Montgomery, Douglas C. 2005. Introduction to Statistical Quality Control. John Wiley & Sons.
"pistonring_data"

#' Semiconductor data
#'
#' A dataset cotaining the 151st feature in SECOM dataset 
#'
#' @format A vector with 50 observations:
#' \describe{
#'   \item{obs}{the 151st feature}
#'  }
#' @references McCann, Michael, and Adrian Johnston. 2008. “SECOM.” UCI Machine Learning Repository. DOI: https://doi.org/10.24432/C54305.
"semiconductor_data"
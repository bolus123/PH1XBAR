#' @keywords package
#' @aliases PH1XBAR-package
"_PACKAGE"

#' Phase I Shewhart X-bar Control Chart
#'
#' The utility of this package is in building a Shewhart-type control chart based on new methods for subgrouped and individual data. The Phase I chart is based on the multivariate normal/t or ARMA process.
#'
#' @name PH1XBAR-package
#' @references Champ, C.W., and Jones, L.A. (2004) Designing Phase I X-bar charts with small sample sizes. Quality and Reliability Engineering International. 20(5), 497-510
#' @references Yao, Y., Hilton, C.W., and Chakraborti, S. (2017) Designing Phase I Shewhart X-bar charts: Extended tables and software. Quality and Reliability Engineering International. 33(8), 2667-2672.
#' @references Yao, Y., and Chakraborti, S. (2021). Phase I monitoring of individual normal data: Design and implementation. Quality Engineering, 33(3), 443-456.
#' @references Yao, Y., and Chakraborti, S. (2021). Phase I process monitoring: The case of the balanced one-way random effects model. Quality and Reliability Engineering International, 37(3), 1244-1265.
#' @references Yao, Y., Chakraborti, S., Yang, X., Parton, J., Lewis Jr, D., and Hudnall, M. (2023). Phase I control chart for individual autocorrelated data: application to prescription opioid monitoring. Journal of Quality Technology, 55(3), 302-317.
#' @import forecast mvtnorm pracma
#' @importFrom stats qt
#' @importFrom stats dchisq
#' @importFrom stats pnorm
#' @importFrom stats qchisq
#' @importFrom stats qnorm
#' @importFrom stats runif
#' @importFrom stats arima
#' @importFrom stats arima.sim
#' @importFrom stats dnorm
#' @importFrom stats integrate
#' @importFrom stats rbeta
#' @importFrom stats rchisq
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats rt
#' @importFrom stats sd
#' @importFrom stats uniroot
#' @importFrom stats var
#' @importFrom graphics axis
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics text
#' @examples
#' #Build a Phase I basic Shewhart control chart
#' data(grinder_data)
#' PH1XBAR(grinder_data, nsim=10)
#' 
#' # Build a Phase I individual control chart with an ARMA model
#' data(preston_data)
#' PH1ARMA(preston_data, nsim.process=10, nsim.coefs=10)
#' 
NULL
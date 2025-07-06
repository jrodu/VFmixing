#' Example Data
#'
#' An example data set
#'
#' @format ## `ExData`
#' A matrix with 501 rows and 3 columns.
#' 
#' @source The example in the RungeKutta function with some error.
"ExData"

#' Example Path
#'
#' An example path
#'
#' @format ## `ExPath`
#' A list with 2 elements: 
#' \describe{
#'   \item{Path}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#' }
#' 
#' @source N/A
"ExPath"

#' Example Trajectory
#'
#' An example trajectory
#'
#' @format ## `ExTraj`
#' A list with 2 elements: 
#' \describe{
#'   \item{Traj}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#' }
#' 
#' @source N/A
"ExTraj"

#' Example Estimated Path
#'
#' A path estimated using the ByFours method.
#'
#' @format ## `ExEstPath_ByFours`
#' A list with 5 elements: 
#' \describe{
#'   \item{Path}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#'   \item{n_splits}{The number of regions in the cubic spline.}
#'   \item{Error}{A matrix with the squared residuals in each dimension for each data point.}
#'   \item{PosSplits}{A vector of boolean describing which regions in the spline can be split in the next iteration.}
#' }
#' 
#' @source The example in the getEstimatedPaths_ByFours function.
"ExEstPath_ByFours"

#' Example Estimated Path
#'
#' A path estimated using the FullRegression method.
#'
#' @format ## `ExEstPath_FullRegression`
#' A list with 11 elements: 
#' \describe{
#'   \item{Path}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#'   \item{n_splits}{The number of regions in the cubic spline.}
#'   \item{Error}{A matrix with the squared residuals in each dimension for each data point.}
#'   \item{PosSplits}{A vector of boolean describing which regions in the spline can be split in the next iteration.}
#'   \item{Cov}{A list of lists containing covariance matrices for each region and dimension of the cubic spline coefficients.}
#'   \item{Sigma2}{A list of vectors containing the sigma2 estimates for each region and dimension of the cubic spline.}
#'   \item{SplitLabels}{A factor vector describing which region of the spline each data point lies.}
#'   \item{FullCov}{A list a covariance matrices for the cubic spline coefficients for each dimension.}
#'   \item{Free_Parameters}{A matrix containing the estimates of the free coefficients for each dimension's spline.}
#'   \item{Free_Parameter_Cov}{A list of covariance matrices for the free coefficients for each dimension's spline.}
#' }
#' 
#' @source The example in the getEstimatedPaths_FullRegression function.
"ExEstPath_FullRegression"

#' Example Estimated Path
#'
#' A path estimated using the Projection method.
#'
#' @format ## `ExEstPath_Projection`
#' A list with 14 elements: 
#' \describe{
#'   \item{Path}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#'   \item{n_splits}{The number of regions in the cubic spline.}
#'   \item{Error}{A matrix with the squared residuals in each dimension for each data point.}
#'   \item{PosSplits}{A vector of boolean describing which regions in the spline can be split in the next iteration.}
#'   \item{Cov}{A list of lists containing covariance matrices for each region and dimension of the cubic spline coefficients.}
#'   \item{Sigma2}{A list of vectors containing the sigma2 estimates for each region and dimension of the cubic spline.}
#'   \item{SplitLabels}{A factor vector describing which region of the spline each data point lies.}
#'   \item{FullCov}{A list a covariance matrices for the cubic spline coefficients for each dimension.}
#'   \item{Free_Parameters}{A matrix containing the estimates of the free coefficients for each dimension's spline.}
#'   \item{Free_Parameter_Cov}{A list of covariance matrices for the free coefficients for each dimension's spline.}
#'   \item{UnrestrictedCoef}{A list of unrestricted cubic coefficents without any induced smoothness in the spline.}
#'   \item{UnrestrictedCov}{A list of lists containing covariance matrices for each region and dimension of the unrestricted cubic spline coefficients.}
#'   \item{IntegralList}{A list of matrices used to project the unrestricted coefficients into the restricted space.}
#' }
#' 
#' @source The example in the getEstimatedPaths_Projection function.
"ExEstPath_Projection"

#' Example Estimated Trajectory
#'
#' A trajectory estimated using the ByFours method.
#'
#' @format ## `ExEstTraj_ByFours`
#' A list with 5 elements: 
#' \describe{
#'   \item{Traj}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#'   \item{n_splits}{The number of regions in the cubic spline.}
#'   \item{Error}{A matrix with the squared residuals in each dimension for each data point.}
#'   \item{PosSplits}{A vector of boolean describing which regions in the spline can be split in the next iteration.}
#' }
#' 
#' @source The example in the getEstimatedTrajectories_ByFours function.
"ExEstTraj_ByFours"

#' Example Estimated Trajectory
#'
#' A trajectory estimated using the Optimization method.
#'
#' @format ## `ExEstTraj_Optimization`
#' A list with 9 elements: 
#' \describe{
#'   \item{Traj}{A list of cubic coefficient matrices.}
#'   \item{TimeSplits}{A vector of node, starting, and ending times.}
#'   \item{n_splits}{The number of regions in the cubic spline.}
#'   \item{Error}{A vector with the integral squared errors for each region of the spline.}
#'   \item{PosSplits}{A vector of boolean describing which regions in the spline can be split in the next iteration.}
#'   \item{Cov}{A list of lists containing covariance matrices for each region and model of the cubic spline coefficients.}
#'   \item{Sigma2}{A list of vectors containing the sigma2 estimates for each region of the cubic spline.}
#'   \item{IntegralList_Traj}{A list of matrices used to transform the path coefficients into the trajectory coefficients.}
#'   \item{IntegralList_PathTraj}{A list of matrices used to transform the path coefficients into the trajectory coefficients.}
#' }
#' 
#' @source The example in the getEstimatedTrajectories_Optimization function.
"ExEstTraj_Optimization"

#' Example PropagationList
#'
#' An example propagation list used in the drifterVisualizationTool.
#'
#' @format ## `ExPropagationList`
#' A list with 2 elements: 
#' \describe{
#'   \item{First}{A data frame with the propagation estimates.}
#'   \item{Second}{A list of matrices with the intermediate Runge Kutta steps for each propagation.}
#' }
#' 
#' @source The example in the getPropagationList function.
"ExPropagationList"

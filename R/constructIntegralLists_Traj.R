#' constructIntegralList_Traj
#'
#' @param EstPath an EstimatedPath object.
#' @param baseVectorFields a baseVectorField function.
#' @param cur_integral_list a list of matrices. Each corresponds with a block of the previous overall trajectory integral matrix.
#' @param time_splits a numeric vector. The nodes of the spline with the starting and ending times. 
#' @param worst_split a numeric scalar. The region to be split.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#'
#' @returns a list of matrices.
#' @export
#'
constructIntegralList_Traj = function(EstPath, baseVectorFields, cur_integral_list, time_splits, worst_split, n_models){
  
  EstPathTimeSplits = EstPath$TimeSplits
  
  integrateMatrixElement = function(lower, upper, index_sum, baseVel_Indices){
    
    integrand = function(t){
      
      cur_path_index = min(c(max(c(sum((EstPathTimeSplits - t)<0),1)), length(EstPathTimeSplits)-1))
      cur_EstPos = evaluateCubic(t, EstPath$Path[[cur_path_index]], EstPathTimeSplits[cur_path_index])
      
      cur_baseVel = baseVectorFields(t, cur_EstPos)
      
      (t-lower)^(8-index_sum)*(t(cur_baseVel[,baseVel_Indices[1]]) %*% t(t(cur_baseVel[,baseVel_Indices[2]])))
      
    }
    
    
    err = tryCatch(expr = stats::integrate(f = Vectorize(integrand), lower = lower, upper = upper, subdivisions = 10000)$value, 
                   error = function(e){NULL})
    
    if(is.null(err)){
      err = cubature::cubintegrate(integrand, lower = lower, upper = upper, method = 'pcubature')$integral
    }
    
    err
  }
  
  model_comb = gtools::combinations(n_models, 2, repeats.allowed = T)
  
  first_block = matrix(rep(0, (4*n_models)^2), nrow = 4*n_models)
  
  for(i in 1:nrow(model_comb)){
    
    sub_elements = unlist(apply(matrix(2:8), MARGIN = 1, FUN = integrateMatrixElement, lower = time_splits[worst_split], upper = time_splits[worst_split+1], baseVel_Indices = model_comb[i,]))
    
    sub_block = matrix(sub_elements[c(1,2,3,4,2,3,4,5,3,4,5,6,4,5,6,7)], nrow = 4, byrow = T)
    
    first_block[1:4 + 4*(model_comb[i,1]-1), 1:4 + 4*(model_comb[i,2]-1)] = first_block[1:4 + 4*(model_comb[i,2]-1), 1:4 + 4*(model_comb[i,1]-1)] = sub_block
    
  }
  
  
  second_block = matrix(rep(0, (4*n_models)^2), nrow = 4*n_models)
  
  for(i in 1:nrow(model_comb)){
    
    sub_elements = unlist(apply(matrix(2:8), MARGIN = 1, FUN = integrateMatrixElement, lower = time_splits[worst_split+1], upper = time_splits[worst_split+2], baseVel_Indices = model_comb[i,]))
    
    sub_block = matrix(sub_elements[c(1,2,3,4,2,3,4,5,3,4,5,6,4,5,6,7)], nrow = 4, byrow = T)
    
    second_block[1:4 + 4*(model_comb[i,1]-1), 1:4 + 4*(model_comb[i,2]-1)] = second_block[1:4 + 4*(model_comb[i,2]-1), 1:4 + 4*(model_comb[i,1]-1)] = sub_block
    
  }
  
  cur_integral_list[[worst_split]] = first_block
  cur_integral_list = append(cur_integral_list, list(second_block), after = worst_split)
  
  cur_integral_list
  
}


#' constructIntegralList_PathTraj
#'
#' @param EstPath an EstimatedPath object.
#' @param baseVectorFields a baseVectorFields function.
#' @param cur_integral_list a list of matrices. Each corresponds with a block of the previous overall PathTrajectory integral matrix.
#' @param worst_split a numeric scalar. The region to be split.
#' @param Path_TimeSplits a numeric vector. The nodes of the path spline with the starting and ending times.
#' @param Traj_TimeSplits a numeric vector. The nodes of the currrent trajectory spline with the starting and ending times.
#' @param n_dim a numeric scalar. The number of dimensions in the path.
#' @param n_models a numeric scalar. The number of vector fields being mixed.
#'
#' @returns a list of matrices.
#' @export
#'
constructIntegralList_PathTraj = function(EstPath, baseVectorFields, cur_integral_list, worst_split, Path_TimeSplits, Traj_TimeSplits, n_dim, n_models){
  
  n_PathRegions = length(Path_TimeSplits) - 1
  n_TrajRegions = length(Traj_TimeSplits) - 1
  
  Total_Col = 4*n_models*n_TrajRegions
  Total_Row = 4*n_dim*n_PathRegions
  
  integrateMatrixCol = function(col){
    
    Traj_Region = ceiling(col/(4*n_models))
    
    lower_path_region = min(c(max(c(sum((Path_TimeSplits - Traj_TimeSplits[Traj_Region])<0),1)), length(Path_TimeSplits)-1))
    upper_path_region = min(c(max(c(sum((Path_TimeSplits - Traj_TimeSplits[Traj_Region+1])<0),1)), length(Path_TimeSplits)-1))
    
    Overlaps = lower_path_region:upper_path_region
    
    integrateMatrixElement = function(row){
      
      Path_Region = ceiling(row/(4*n_dim))
      
      lower = sort(c(Path_TimeSplits[0:1 + Path_Region], Traj_TimeSplits[0:1 + Traj_Region]))[2]
      upper = sort(c(Path_TimeSplits[0:1 + Path_Region], Traj_TimeSplits[0:1 + Traj_Region]))[3]
      
      d = ceiling((row %% (4*n_dim))/4)
      d = d + n_dim*(d == 0)
      
      m = ceiling((col %% (4*n_models))/4)
      m = m + n_models*(m == 0)
      
      poly_path = row %% 4
      poly_path = poly_path + 4*(poly_path == 0)
      
      poly_traj = col %% 4
      poly_traj = poly_traj + 4*(poly_traj == 0)
      
      
      integrand = function(t){
        
        cur_pos = getEstPathPosition(t, EstPath)
        
        VF = baseVectorFields(t, cur_pos)[d,m]
        Path_Polynomial = (4-poly_path)*(t-Path_TimeSplits[Path_Region])^(max(c(3-poly_path,0)))
        Traj_Polynomial = (t-Traj_TimeSplits[Traj_Region])^(4-poly_traj)
        
        if(is.na(Path_Polynomial) | is.na(Traj_Polynomial)){
          return(0)
        } else{
          return(VF*Path_Polynomial*Traj_Polynomial)
        }
        
      }
      
      err = tryCatch(expr = stats::integrate(f = Vectorize(integrand), lower = lower, upper = upper, subdivisions = 10000)$value,
                     error = function(e){NULL})
      
      if(is.null(err)){
        err = cubature::cubintegrate(integrand, lower = lower, upper = upper, method = 'pcubature')$integral
      }
      
      err
      
    }
    
    c(rep(0, 4*n_dim*(Overlaps[1]-1)), apply(matrix((4*n_dim*(Overlaps[1]-1)+1):(4*n_dim*(Overlaps[length(Overlaps)]))), MARGIN = 1, FUN = integrateMatrixElement), rep(0,Total_Row - (4*n_dim*(Overlaps[length(Overlaps)]))))
    
  }
  
  
  first_Matrix = apply(matrix(1:(4*n_models) + 4*n_models*(worst_split-1)), MARGIN = 1, FUN = integrateMatrixCol, simplify = T)
  
  second_Matrix = apply(matrix(1:(4*n_models) + 4*n_models*(worst_split)), MARGIN = 1, FUN = integrateMatrixCol, simplify = T)
  
  cur_integral_list[[worst_split]] = first_Matrix
  
  integral_list = append(cur_integral_list, list(second_Matrix), after = worst_split)
  
  integral_list
  
}


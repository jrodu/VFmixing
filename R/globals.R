# Column names used through non-standard evaluation in ggplot2::aes(),
# dplyr verbs, tidyr::pivot_longer() and plotly/leaflet formulas.
# Declaring them here silences R CMD check's "no visible binding" NOTE.
utils::globalVariables(c(
  # ggplot2 plotting functions
  "X1", "X2", "X1v", "X2v", "upper", "lower",
  "Trajectory", "Model", "Angle", "Sample", "x", "y", "id",
  # drifterVisualizationTool
  "time", "time_id", "time_diff", "type", "lat", "lon",
  "true_lat", "true_lon", "point_id", "method_index", "method_name"
))

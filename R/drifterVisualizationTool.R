#' drifterVisualizationTool
#' 
#' @description
#' A shiny application to take a deep dive into both paths and trajectories.
#' 
#'
#' @param path an EstimatedPath object.
#' @param traj_list a list of EstimatedTrajectory objects.
#' @param baseVectorFields a baseVectorFields function.
#' @param propagation_list a Propagation List. 
#' @param path_t_grid_size a numeric scalar. The time grid size used to plot the path.
#' @param traj_t_grid_size a numeric scalar. The time grid size used to plot the trajectories.
#' @param model_names a string vector. The names of the mixing vector fields
#' @param method_names a string vector. The names of the methods used to compute the trajectories.
#' @param model_colors a string vector. The names or hexcodes of colors to plot the trajectories.
#' @param method_colors a string vector. The names or hexcodes of colors to plot the methods.
#' @param vector_field_arrow_scale a numeric scalar. The scale used to draw the arrows in Vector Field Mode.
#' @param vector_field_grid_resolution a numeric scalar. The number of arrows in either dimension to display the vector fields in Vector Field Mode.
#'
#' @returns a shiny app.
#' @export
#'
drifterVisualizationTool = function(path, traj_list, baseVectorFields, propagation_list, 
                                    path_t_grid_size, traj_t_grid_size, 
                                    model_names, method_names,
                                    model_colors, method_colors, 
                                    vector_field_arrow_scale, vector_field_grid_resolution){
  
  if (!requireNamespace("leaflet","leaflet.extras2","plotly","RColorBrewer",
                        "shiny","shinyjs","tidyr","waiter", quietly = TRUE)) {
    stop(
      "Packages \"leaflet\",\"leaflet.extras2\",\"plotly\",\"RColorBrewer\",
                        \"shiny\",\"shinyjs\",\"tidyr\", and \"waiter\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  path_t_seq = seq(min(path$TimeSplits), max(path$TimeSplits), by = path_t_grid_size)
  path_pos = getEstPathPosition(path_t_seq, path)
  
  path_highres_df <- data.frame(
    time = as.POSIXct("1993-01-01 00:00:00", tz = "UTC") + path_t_seq*24*60*60,
    lat = path_pos[,2],
    lon = path_pos[,1]
  )
  
  # --- DATA 2: High-resolution data for the subplot functions ---
  
  min_traj_time = min(sapply(traj_list, function(x){min(x$TimeSplits)}))
  max_traj_time = min(sapply(traj_list, function(x){max(x$TimeSplits)}))
  traj_t_seq = seq(min_traj_time, max_traj_time, traj_t_grid_size)
  
  Traj_Vals = do.call(cbind, lapply(traj_list, FUN = function(x){getEstTrajValue(traj_t_seq,x)}))
  
  subplot_highres_df = data.frame(time = as.POSIXct("1993-01-01 00:00:00", tz = "UTC") + traj_t_seq*24*60*60, 
                                  Traj_Vals)
  traj_names_grid = expand.grid(model_names, method_names)
  
  names(subplot_highres_df) = c('time', paste(traj_names_grid[,2], traj_names_grid[,1], sep = "_"))
  
  # --- DATA 3: Sparse "Error Points" - The link between map and plots ---
  
  propagation_list_finalpos =  propagation_list[[1]]
  propagation_inter_locations = propagation_list[[2]]
  
  error_times = propagation_list_finalpos[,1]
  true_pos = getEstPathPosition(error_times, path)
  
  n_dim = ncol(true_pos)
  n_methods = length(method_names)
  n_models = length(model_names)
  
  propagation_final_locations = do.call(cbind, apply(matrix(1:n_methods), 1, function(i){true_pos - propagation_list_finalpos[,1:n_dim + 1 + (i-1)*n_dim]}))
  
  # Generate the error point data
  
  error_points_df = data.frame(
    
    time = as.POSIXct("1993-01-01 00:00:00", tz = "UTC") + error_times*24*60*60,
    time_id = 1:length(error_times),
    true_lon = true_pos[,1],
    true_lat = true_pos[,2],
    propagation_final_locations
  )
  
  prop_names_grid = as.matrix(expand.grid(method_names, c('lon','lat')))
  
  names(error_points_df) = c('time','time_id','true_lon','true_lat',paste(prop_names_grid[,1], prop_names_grid[,2], sep = "_"))
  
  # --- Combine Data for Plotting on the Map ---
  # This data is used for all the interactive markers on the map
  
  method_names_w_true = c("Anchor Position", method_names)
  
  combined_map_data <- dplyr::bind_rows(apply(matrix(1:(n_methods+1)), MARGIN = 1, FUN = function(i){error_points_df %>% dplyr::select(time_id, time, lon = !!dplyr::sym(names(error_points_df[3 + (i-1)*n_dim])), lat = !!dplyr::sym(names(error_points_df[4 + (i-1)*n_dim]))) %>% ggpubr::mutate(type = method_names_w_true[i])})) %>%
    ggpubr::mutate(unique_marker_id = paste(type, time_id, sep = "-"))
  
  vis_method_names = c(method_names, model_names)
  vis_traj_list = traj_list
  
  for(i in 1:n_models){
    
    cur_traj = list(matrix(c(rep(0,(i-1)*4),
                             rep(0,3),1,
                             rep(0, (n_models-i)*4)), nrow = n_models, byrow = T))
    vis_traj_list[[length(vis_traj_list)+1]] = list(Traj = cur_traj, TimeSplits = c(0,1))
    
    
  }
  
  
  all_vis_method_names <- c("All Models", vis_method_names) 
  
  
  # 4. Define the UI for the Shiny App
  
  ui <- shiny::fluidPage(
    waiter::use_waiter(), # 1. Set up the waiter resources
    shinyjs::useShinyjs(),
    shiny::titlePanel("Drifter Visualization Tool"),
    shiny::fluidRow(
      shiny::column(7,
             shiny::div(id = "map_overlay", class = "overlay-pane"),
             leaflet::leafletOutput("map", height = "85vh"),
             
             # Use absolutePanel for the button controls
             shiny::absolutePanel(
               top = 5, right = 20,
               width = "170px", # Give the panel a fixed width
               draggable = TRUE,
               
               # This is the main button that will toggle the mode
               shiny::actionButton("toggle_mode_button", "Enter Vector Field Mode", width = "100%"),
               
               # This div contains the extra buttons. It's hidden by default.
               shinyjs::hidden(
                 shiny::div(
                   id = "vector_mode_controls",
                   style = "margin-top: 5px; display: flex; justify-content: space-between;",
                   shiny::actionButton("back_button", "Back", width = "48%"),
                   shiny::actionButton("next_button", "Next", width = "48%")
                 )
               ),
               shinyjs::hidden(
                 shiny::div(
                   id = "trajectory_selector_panel",
                   style = "margin-top: 10px; background-color: rgba(255,255,255,0.8); padding: 5px; border-radius: 5px;",
                   shiny::radioButtons(
                     inputId = "trajectory_selector",
                     label = "Vector Field Source:",
                     # We'll use method_names as the source for each trajectory.
                     # In a real app, these names might come from your traj_list.
                     choices = all_vis_method_names[1:(length(all_vis_method_names) - n_models)],
                     selected = all_vis_method_names[1:(length(all_vis_method_names) - n_models)][1]
                   )
                 )
               )
             )
      ),
      shiny::column(5, 
             shiny::div(id = "plot_overlay", class = "overlay-pane"),
             plotly::plotlyOutput("subplots_combined", height = "85vh"))
    )
  )
  
  # 5. Define the Server logic for the Shiny App
  server <- function(input, output, session) {
    
    waiter::waiter_show(
      html = shiny::tagList(
        # The background div
        tags$div(
          style = "
        position: fixed; 
        top: 0; 
        left: 0; 
        width: 100%; 
        height: 100%; 
        background-image: url('https://media2.giphy.com/media/v1.Y2lkPTc5MGI3NjExa2FxaTRzNjRpaHFxZGc5bmlsOXoyMXNnZWFsYmU1anVhYjlkczY2MyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/PmXy0k63DhdDjVQBgQ/giphy.gif'); 
        background-size: cover; 
        background-position: center;
        z-index: 1; /* Push it to the back */
      "
        ),
        # Your custom text
        shiny::h1(
          "Loading Application...",
          style = "
        position: relative; /* Required for z-index to work */
        z-index: 2; /* Bring it to the front */
        color: white; 
        font-family: 'Georgia', serif;
        font-weight: bold;
        text-shadow: 2px 2px 8px #000000;
      "
        )
      )
    )
    
    shinyjs::runjs("
    $('.action-button').on('click', function() {
      var btn = $(this);
      setTimeout(function() {
        btn.blur();
      }, 250); // 1000ms = 1 second
    });
  ")
    
    # --- State Management ---
    selected_time_id <- shiny::reactiveVal(NULL)
    hovered_time <- shiny::reactiveVal(NULL)
    event_source <- shiny::reactiveVal(NULL)
    map_zoom_time_range <- shiny::reactiveVal(NULL)
    vector_mode_on <- shiny::reactiveVal(FALSE)
    subplot_hover_event <- shiny::debounce(shiny::reactive(plotly::event_data("plotly_hover")), 1)
    
    ### NEW: Reactive value to store the count of visible time IDs.
    # Initialize with the total number of points.
    visible_time_id_count <- shiny::reactiveVal(length(unique(error_points_df$time_id)))
    visible_error_points <- shiny::reactive({
      bounds <- shiny::req(input$map_bounds)
      error_points_df %>%
        dplyr::filter(
          true_lat >= bounds$south, true_lat <= bounds$north,
          true_lon >= bounds$west, true_lon <= bounds$east
        ) %>%
        dplyr::arrange(time_id)
    })
    selected_trajectory_for_field <- shiny::reactiveVal(NULL)
    precomputed_vectors <- shiny::reactiveVal(list())
    
    #Computing the Model Vector Fields
    
    computeBaseFieldsForTime <- function(time_id_to_calc) {
      bounds <- shiny::req(input$map_bounds)
      selected_time_posix <- error_points_df$time[error_points_df$time_id == time_id_to_calc]
      
      if (length(selected_time_posix) == 0) return(NULL)
      
      selected_time_numeric <- as.numeric(difftime(selected_time_posix, as.Date('1993-01-01'), units = 'days'))
      
      grid_resolution <- vector_field_grid_resolution
      lon_seq <- seq(bounds$west, bounds$east, length.out = grid_resolution)
      lat_seq <- seq(bounds$south, bounds$north, length.out = grid_resolution)
      grid_df <- expand.grid(lon = lon_seq, lat = lat_seq)
      
      # Call the base function once for each grid point and store the raw result
      base_fields_grid <- lapply(1:nrow(grid_df), function(i) {
        baseVectorFields(t = selected_time_numeric, curPos = c(grid_df[i,1], grid_df[i,2]))
      })
      
      # Return both the grid locations and the calculated base vectors
      return(list(grid_df = grid_df, base_fields_grid = base_fields_grid, selected_time_numeric = selected_time_numeric))
    }
    
    ##Calculate Vector Fields
    
    calculateVectorField <- function(base_fields_data, trajectory_index, traj_list) {
      if (is.null(base_fields_data)) return(NULL)
      
      grid_df <- base_fields_data$grid_df
      base_fields_grid <- base_fields_data$base_fields_grid
      selected_time_numeric <- base_fields_data$selected_time_numeric
      
      # Calculate the final vectors by applying the weights to the pre-computed base fields
      vectors_matrix <- t(sapply(base_fields_grid, function(base_vectors) {
        
        cur_traj_value = getEstTrajValue(t = selected_time_numeric, EstTraj = traj_list[[trajectory_index]])
        
        base_vectors %*% t(cur_traj_value)
      }))
      
      all_lines <- list()
      arrow_scale <- vector_field_arrow_scale
      
      # This part for creating arrow geometry remains unchanged
      for (i in 1:nrow(grid_df)) {
        start_lon <- grid_df$lon[i]
        start_lat <- grid_df$lat[i]
        u <- vectors_matrix[i, 1] * arrow_scale
        v <- vectors_matrix[i, 2] * arrow_scale
        end_lon <- start_lon + u
        end_lat <- start_lat + v
        
        shaft <- data.frame(lon = c(start_lon, end_lon), lat = c(start_lat, end_lat), id = paste0("arrow_", i, "_shaft"))
        all_lines[[length(all_lines) + 1]] <- shaft
        
        angle <- atan2(v, u)
        arrowhead_length <- 0.3 * sqrt(u^2 + v^2)
        arrowhead_angle <- pi / 6
        
        head1_angle <- angle + pi + arrowhead_angle
        head1 <- data.frame(lon = c(end_lon, end_lon + arrowhead_length * cos(head1_angle)), lat = c(end_lat, end_lat + arrowhead_length * sin(head1_angle)), id = paste0("arrow_", i, "_head1"))
        all_lines[[length(all_lines) + 1]] <- head1
        
        head2_angle <- angle + pi - arrowhead_angle
        head2 <- data.frame(lon = c(end_lon, end_lon + arrowhead_length * cos(head2_angle)), lat = c(end_lat, end_lat + arrowhead_length * sin(head2_angle)), id = paste0("arrow_", i, "_head2"))
        all_lines[[length(all_lines) + 1]] <- head2
      }
      
      return(all_lines)
    }
    
    ## Draw Vector Field Function
    
    drawVectorField <- function(time_id_to_draw, trajectory_name_to_draw) {
      proxy <- leaflet::leafletProxy("map", session)
      proxy %>% leaflet::clearGroup("vector_field_layer")
      
      if (trajectory_name_to_draw == "All Models") {
        
        # Define a very small jitter amount based on map bounds
        bounds <- shiny::req(input$map_bounds)
        #jitter_amount <- (bounds$east - bounds$west) / 1000 
        jitter_amount = 0
        num_grid_points <- vector_field_grid_resolution^2
        
        # --- Interleaving Loop ---
        # 1. Loop through each grid point first.
        for (grid_idx in 1:num_grid_points) {
          
          # 2. Then, loop through each model for that specific grid point.
          for (model_idx in 1:n_models) {
            model_name <- model_names[model_idx]
            model_color <- model_colors[model_idx]
            
            # Retrieve the pre-computed arrow segments for this model
            # The list contains all segments, so we need to find the ones for this grid point
            all_segments_for_model <- precomputed_vectors()[[model_name]][[as.character(time_id_to_draw)]]
            
            if (is.null(all_segments_for_model)) next # Skip if no data
            
            # Isolate the 3 segments for the current grid point's arrow
            # The segments are ordered, so we can calculate their indices
            segment_indices <- ((grid_idx - 1) * 3 + 1):(grid_idx * 3)
            arrow_segments <- all_segments_for_model[segment_indices]
            
            # 3. Apply a tiny random jitter to the arrow's coordinates
            jitter_lon <- 0
            jitter_lat <- 0
            
            # 4. Draw this single, jittered arrow
            for (segment in arrow_segments) {
              if (!is.null(segment)) {
                jittered_segment <- segment
                jittered_segment$lon <- jittered_segment$lon + jitter_lon
                jittered_segment$lat <- jittered_segment$lat + jitter_lat
                
                proxy %>% leaflet::addPolylines(
                  data = jittered_segment,
                  lng = ~lon, lat = ~lat,
                  group = "vector_field_layer",
                  layerId = paste0(model_name, "_", segment$id[1]), # Ensure unique layerId
                  color = model_color,
                  weight = 2, opacity = 0.5
                )
              }
            }
          }
        }
        
      } else {
        # This is the original logic for drawing a single vector field. It remains unchanged.
        lines_to_draw <- precomputed_vectors()[[trajectory_name_to_draw]][[as.character(time_id_to_draw)]]
        if (is.null(lines_to_draw) || length(lines_to_draw) == 0) return()
        for (line_segment in lines_to_draw) {
          proxy %>% leaflet::addPolylines(
            data = line_segment,
            lng = ~lon, lat = ~lat,
            group = "vector_field_layer", layerId = ~id,
            color = "#003366", weight = 2, opacity = 0.5
          )
        }
      }
    }
    
    drawIntermediatePoints <- function(time_id_to_draw) {
      proxy <- leaflet::leafletProxy("map", session)
      
      # Clear any previously drawn intermediate points
      proxy %>% leaflet::clearGroup("intermediate_points_layer")
      
      # Get the correct dataframe from the list
      inter_df <- propagation_inter_locations[[time_id_to_draw]]
      
      if (is.null(inter_df) || nrow(inter_df) == 0) return()
      
      # --- Step 1: Standardize Column Names (The Fix) ---
      # Programmatically create a standard set of names.
      # This assumes the columns are ordered: time, lon_meth1, lat_meth1, lon_meth2, lat_meth2, etc.
      new_colnames <- c(
        "time", 
        paste(
          rep(c("lon", "lat"), n_methods), 
          rep(1:n_methods, each = 2), 
          sep = "_"
        )
      )
      
      # Assign the new, predictable names to the dataframe
      colnames(inter_df) <- new_colnames
      
      # --- Step 2: Reshape the Data ---
      # Now pivot_longer is much simpler and more reliable.
      long_inter_df <- inter_df %>%
        # Add a row number to identify each point in the sequence
        ggpubr::mutate(point_id = dplyr::row_number()) %>%
        tidyr::pivot_longer(
          cols = -c(time, point_id),
          # Use names_sep instead of a complex pattern
          names_to = c(".value", "method_index"), 
          names_sep = "_"
        ) %>%
        # Assign the correct method name based on the index
        ggpubr::mutate(
          method_name = method_names[as.numeric(method_index)]
        )
      
      # --- Step 3: Drawing Logic (No Changes Here) ---
      # Loop through each method to assign the correct color
      for (i in 1:n_methods) {
        method_data <- dplyr::filter(long_inter_df, method_name == method_names[i])
        
        if (nrow(method_data) > 0) {
          proxy %>% leaflet::addCircleMarkers(
            data = method_data,
            lng = ~lon,
            lat = ~lat,
            group = "intermediate_points_layer",
            layerId = ~paste0("inter_", method_name, "_", point_id),
            color = method_colors[i],
            radius = 3,
            weight = 1,
            fillOpacity = 0.8, 
            stroke = TRUE
          )
        }
      }
    }
    
    updateMapAndPlotHighlights <- function(current_id, current_hover_time) {
      # This is the same logic from your existing highlight observer
      
      # 1. Prepare static and dynamic shapes for the plot
      static_shapes <- lapply(1:n_methods, function(i) {
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 1, y1 = 1, 
             yref = if(i==1) "y" else paste0("y", i), line = list(color = "grey", dash = "dash"))
      })
      dynamic_shapes <- static_shapes
      
      # 2. Clear previous highlights from the map
      proxy <- leaflet::leafletProxy("map", session)
      proxy %>% leaflet::clearGroup("highlight_layer") %>% leaflet::clearGroup("hover_layer")
      
      # 3. Draw blue selection highlight (if a point is selected)
      if (!is.null(current_id)) {
        highlight_data <- combined_map_data %>% dplyr::filter(time_id == current_id[1])
        proxy %>% leaflet::addCircleMarkers(data = highlight_data, lng = ~lon, lat = ~lat, group = "highlight_layer", 
                                   radius = 10, color = "blue", stroke = TRUE, weight = 1, fillOpacity = 0.05)
        
        selected_time <- error_points_df$time[error_points_df$time_id == current_id[1]]
        if (length(selected_time) > 0) {
          blue_line <- list(type = "line", x0 = selected_time, x1 = selected_time, y0 = 0, y1 = 1, 
                            yref = "paper", line = list(color = "blue", dash = "dash"))
          dynamic_shapes <- append(dynamic_shapes, list(blue_line))
        }
      }
      
      # 4. Draw red hover highlight (if a point is hovered)
      if (!is.null(current_hover_time)) {
        hover_point <- path_highres_df %>% dplyr::filter(time == current_hover_time)
        proxy %>% leaflet::addCircleMarkers(data = hover_point, lng = ~lon, lat = ~lat, group = "hover_layer", 
                                   radius = 10, color = "red", stroke = TRUE, weight = 1, fillOpacity = 0.8, 
                                   label = ~paste("Time:", time, "| Lon:", round(lon, 4), "| Lat:", round(lat, 4)))
        
        red_line <- list(type = "line", x0 = current_hover_time, x1 = current_hover_time, y0 = 0, y1 = 1, 
                         yref = "paper", line = list(color = "red", dash = "dash"))
        dynamic_shapes <- append(dynamic_shapes, list(red_line))
      }
      
      # 5. Update the plot with all shapes
      plotly::plotlyProxy("subplots_combined", session) %>%
        plotly::plotlyProxyInvoke("relayout", list(shapes = dynamic_shapes))
    }
    
    # --- Initial Map and Plot Rendering (No changes in this section) ---
    
    output$map <- leaflet::renderLeaflet({
      
      map_types <- sort(c("Anchor Position", method_names))
      map_colors <- c("grey", method_colors)[order(c("Anchor Position", method_names))]
      
      map_pal <- leaflet::colorFactor(
        palette = map_colors,
        domain = map_types
      )
      
      # Initial map setup remains the same
      map <- leaflet::leaflet(data = combined_map_data) %>%
        leaflet::fitBounds(lng1 = min(path_highres_df$lon), lat1 = min(path_highres_df$lat),
                  lng2 = max(path_highres_df$lon), lat2 = max(path_highres_df$lat)) %>%
        leaflet::addProviderTiles(providers$Esri.OceanBasemap, group = "Ocean Basemap") %>%
        # Add arrowheads as before
        leaflet.extras2::addArrowhead(data = path_highres_df, lng = ~lon, lat = ~lat, color = 'black', stroke = F,
                     options = leaflet.extras2::arrowheadOptions(frequency = "80px", size = "10px", fill = TRUE, opacity = 0.9)) %>%
        leaflet::addLegend("bottomright", colors = c("black", "grey", method_colors,"red"), labels = c("Estimated Path", method_names_w_true,"Hover"), title = "Legend")
      
      # This replaces your single addPolylines call for the high-res path
      for (i in 1:(nrow(path_highres_df) - 1)) {
        map <- map %>%
          leaflet::addPolylines(
            data = path_highres_df[i:(i + 1), ], # Data for just ONE segment
            lng = ~lon,
            lat = ~lat,
            color = "black",
            opacity = 0.8,
            weight = 3,
            # Set the label to the time of the segment's starting point
            label = ~paste("Time:", path_highres_df$time[i]),
            # Provide a unique layerId for each segment
            layerId = paste0("segment_", path_highres_df$time[i]),
          )
      }
      
      map = map %>% 
        leaflet::addCircleMarkers(lng = ~lon, lat = ~lat, layerId = ~unique_marker_id,
                         color = ~map_pal(type),
                         radius = 5, stroke = FALSE, fillOpacity = 0.9,
                         popup = ~paste("<b>Type:</b>", type, "<br>", "<b>Time:</b>", time,
                                        "<br>", "<b>Longitude:</b>", round(lon,4), "<br>", "<b>Latitude:</b>", round(lat, 4)))
      
      # Return the final map object
      map
    })
    
    output$subplots_combined <- plotly::renderPlotly({
      subplot_list = list()
      
      for(i in 1:n_methods){
        
        cur_plot = plotly::plot_ly(subplot_highres_df, x = ~time)
        
        for(j in 1:n_models){
          
          cur_plot = cur_plot %>% 
            
            plotly::add_lines(y = stats::as.formula(paste0("~`", names(subplot_highres_df)[1 + (i-1)*n_models + j], "`")), name = model_names[j], legendgroup = paste0("group",j), showlegend = (i == 1), line = list(color = model_colors[j]))
          
        }
        
        subplot_list[[i]] = cur_plot %>%
          plotly::layout(yaxis = list(title = method_names[i]), 
                 shapes = list(list(
                   type = "line",
                   x0 = 0, x1 = 1, xref = "paper", # Make it span the full plot width
                   y0 = 1, y1 = 1, yref = "y", # Reference the data axis
                   line = list(color = "grey", dash = "dash")
                 )))
        
      }
      
      fig <- plotly::subplot(subplot_list, nrows = length(subplot_list), shareX = TRUE, titleY = TRUE)
      
      layout_args <- list(
        title = "Trajectory Functions",
        legend = list(title = list(text = "Currents")),
        showlegend = TRUE 
      )
      
      for(i in 1:length(subplot_list)) {
        y_axis_name <- if(i == 1) "yaxis" else paste0("yaxis", i)
        
        if(i == 1){
          layout_args[[y_axis_name]] <- list(title = method_names[i], fixedrange = TRUE)
        } else{
          layout_args[[y_axis_name]] <- list(title = method_names[i], matches = 'y')
        }
        
        if(i > 1) {
          x_axis_name <- paste0("xaxis", i)
          layout_args[[x_axis_name]] <- list(matches = 'x')
        }
      }
      
      fig <- do.call(plotly::layout, c(list(p = fig), layout_args))
      
      plotly::event_register(fig, "plotly_relayout")
      plotly::event_register(fig, "plotly_click")
      plotly::event_register(fig, "plotly_hover")
      waiter::waiter_hide()
      fig
    })
    
    # --- Observers for Interactivity (No changes in this section) ---
    
    shiny::observeEvent(input$map_marker_click, {
      click_id <- input$map_marker_click$id
      time_id <- combined_map_data$time_id[combined_map_data$unique_marker_id == click_id]
      selected_time_id(time_id)
    })
    
    shiny::observeEvent(plotly::event_data("plotly_click"), {
      
      # First, check if vector mode is currently on.
      if (vector_mode_on()) {
        
        # If it is, show a warning message and stop immediately.
        shiny::showNotification(
          "Please exit Vector Field Mode to interact with the plots.",
          type = "warning",
          duration = 3
        )
        return()
        
      }
      
      # If vector mode is off, run the original logic to select a time point.
      event <- plotly::event_data("plotly_click")
      if (!is.null(event$x)) {
        clicked_time <- as.POSIXct(event$x, origin = "1970-01-01", tz = "UTC")
        closest_point <- error_points_df %>%
          ggpubr::mutate(time_diff = abs(difftime(time, clicked_time, units = "secs"))) %>%
          dplyr::filter(time_diff == min(time_diff)) %>%
          dplyr::slice(1)
        selected_time_id(closest_point$time_id)
      }
    })
    
    shiny::observeEvent(input$map_click, {
      # First, check if vector mode is currently on.
      if (vector_mode_on()) {
        
        # If it is, show a warning message and stop immediately.
        shiny::showNotification(
          "Please exit Vector Field Mode to interact with the map.",
          type = "warning",
          duration = 3 
        )
        return()
        
      } else {
        
        # If vector mode is off, run the original code to clear selections.
        selected_time_id(NULL)
        hovered_time(NULL)
        
      }
    })
    
    shiny::observeEvent(input$map_shape_mouseover, {
      info <- input$map_shape_mouseover
      
      if (startsWith(info$id, "segment_")) {
        time_str <- sub("segment_", "", info$id)
        hovered_time(as.POSIXct(time_str, tz = "UTC"))
      }
    })
    
    shiny::observeEvent(subplot_hover_event(), {
      event <- subplot_hover_event()
      if (!is.null(event$x)) {
        hover_time <- as.POSIXct(event$x, origin = "1970-01-01", tz = "UTC")
        closest_time <- path_highres_df %>%
          ggpubr::mutate(time_diff = abs(difftime(time, hover_time, units = "secs"))) %>%
          dplyr::filter(time_diff == min(time_diff)) %>%
          dplyr::slice(1)
        
        hovered_time(closest_time$time)
      }
    })
    
    
    ### NEW: Observer to enable/disable the button based on the count.
    shiny::observe({
      count <- visible_time_id_count()
      
      # Add a title to the button to inform the user why it's disabled.
      title_text <- if (count > 10) "Zoom in to less than 10 time points to enable" else ""
      
      if (count > 10) {
        shinyjs::disable("toggle_mode_button")
      } else {
        shinyjs::enable("toggle_mode_button")
      }
      
      # Update the button's title attribute
      shinyjs::runjs(sprintf("$('#toggle_mode_button').prop('title', '%s');", title_text))
    })
    
    ### MODIFIED: This observer now also updates the count of visible points.
    shiny::observeEvent(input$map_bounds, {
      if (!is.null(event_source()) && event_source() == "subplot") {
        event_source(NULL)
        return()
      }
      
      bounds <- input$map_bounds
      
      # --- Logic for counting visible points ---
      visible_error_points <- error_points_df %>%
        dplyr::filter(
          true_lat >= bounds$south, true_lat <= bounds$north,
          true_lon >= bounds$west, true_lon <= bounds$east
        )
      count <- length(unique(visible_error_points$time_id))
      visible_time_id_count(count)
      
      # --- Existing logic for zooming subplots ---
      filtered_data <- path_highres_df %>% dplyr::filter(lat >= bounds$south, lat <= bounds$north, lon >= bounds$west, lon <= bounds$east)
      if(nrow(filtered_data) > 0) {
        min_time <- min(filtered_data$time)
        max_time <- max(filtered_data$time)
        
        map_zoom_time_range(c(min_time, max_time))
        
        event_source("map")
      }
    })
    
    
    shiny::observeEvent(list(selected_time_id(), hovered_time()), {
      
      # Add this condition to prevent it from running in vector mode
      if (!vector_mode_on()) {
        updateMapAndPlotHighlights(selected_time_id(), hovered_time())
      }
      
    }, ignoreNULL = FALSE)    
    
    # --- Zooming logic (No changes in this section) ---
    shiny::observeEvent(map_zoom_time_range(), {
      if (!is.null(event_source()) && event_source() == "subplot") {
        event_source(NULL)
        return()
      }
      
      current_zoom_time_range <- map_zoom_time_range()
      
      if(!is.null(current_zoom_time_range)){
        
        zoomed_subplot_df = subplot_highres_df[subplot_highres_df$time >= current_zoom_time_range[1] & subplot_highres_df$time <= current_zoom_time_range[2],]
        current_traj_range = c(min(zoomed_subplot_df[,-1]),max(zoomed_subplot_df[,-1]))
        
        plotly::plotlyProxy("subplots_combined", session) %>% plotly::plotlyProxyInvoke("relayout", list(xaxis = list(range = current_zoom_time_range, autorange = FALSE)))
        
      }
      
    }, ignoreNULL = F)
    
    shiny::observeEvent(plotly::event_data("plotly_relayout"), {
      if (!is.null(event_source()) && event_source() == "map") {
        event_source(NULL)
        return()
      }
      
      event <- plotly::event_data("plotly_relayout")
      
      if (!is.null(event[["xaxis.range[0]"]])) {
        min_time <- as.POSIXct(event[["xaxis.range[0]"]], tz = "UTC")
        max_time <- as.POSIXct(event[["xaxis.range[1]"]], tz = "UTC")
        
        map_data_in_view <- path_highres_df %>%
          dplyr::filter(time >= min_time, time <= max_time)
        
        if (nrow(map_data_in_view) > 0) {
          event_source("subplot")
          leaflet::leafletProxy("map", session) %>%
            leaflet::flyToBounds(
              lng1 = min(map_data_in_view$lon), lat1 = min(map_data_in_view$lat),
              lng2 = max(map_data_in_view$lon), lat2 = max(map_data_in_view$lat)
            )
        }
      } else if (!is.null(event[["xaxis.autorange"]]) && event[["xaxis.autorange"]] == TRUE) {
        event_source("subplot")
        leaflet::leafletProxy("map", session) %>%
          leaflet::flyToBounds(
            lng1 = min(path_highres_df$lon), lat1 = min(path_highres_df$lat),
            lng2 = max(path_highres_df$lon), lat2 = max(path_highres_df$lat)
          )
      }
    })
    
    shiny::observeEvent(hovered_time(), {
      # We want hover to work in both modes, so no req(vector_mode_on())
      
      current_hover_time <- hovered_time()
      
      # --- 1. Update the Map (simple) ---
      proxy <- leaflet::leafletProxy("map", session)
      proxy %>% leaflet::clearGroup("hover_layer") # Clear previous red circle
      
      if (!is.null(current_hover_time)) {
        hover_point <- path_highres_df %>% dplyr::filter(time == current_hover_time)
        proxy %>% leaflet::addCircleMarkers(data = hover_point, lng = ~lon, lat = ~lat, group = "hover_layer", 
                                   radius = 10, color = "red", stroke = TRUE, weight = 1, fillOpacity = 0.8, 
                                   label = ~paste("Time:", time, "| Lon:", round(lon, 4), "| Lat:", round(lat, 4)))
      }
      
      # --- 2. Update the Plot (rebuilding all shapes) ---
      # Start with the static horizontal lines
      static_shapes <- lapply(1:n_methods, function(i) {
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 1, y1 = 1, 
             yref = if(i==1) "y" else paste0("y", i), line = list(color = "grey", dash = "dash"))
      })
      
      shapes_to_draw <- static_shapes
      
      # Add the blue selection line, if a point is currently selected
      current_id <- selected_time_id()
      if (!is.null(current_id)) {
        selected_time <- error_points_df$time[error_points_df$time_id == current_id[1]]
        if (length(selected_time) > 0) {
          blue_line <- list(type = "line", x0 = selected_time, x1 = selected_time, y0 = 0, y1 = 1, 
                            yref = "paper", line = list(color = "blue", dash = "dash"))
          shapes_to_draw <- append(shapes_to_draw, list(blue_line))
        }
      }
      
      # Add the new red hover line
      if (!is.null(current_hover_time)) {
        red_line <- list(type = "line", x0 = current_hover_time, x1 = current_hover_time, y0 = 0, y1 = 1, 
                         yref = "paper", line = list(color = "red", dash = "dash"))
        shapes_to_draw <- append(shapes_to_draw, list(red_line))
      }
      
      # Send the complete list of shapes to the plot
      plotly::plotlyProxy("subplots_combined", session) %>%
        plotly::plotlyProxyInvoke("relayout", list(shapes = shapes_to_draw))
      
    }, ignoreNULL = FALSE)
    
    # --- Vector Field Mode Button Logic (No changes in this section) ---
    
    shiny::observeEvent(input$toggle_mode_button, {
      vector_mode_on(!vector_mode_on())
      is_on <- vector_mode_on()
      
      # Your existing code to lock/unlock Plotly subplots
      plotly_layout_args <- list('xaxis.fixedrange' = is_on)
      if (n_methods > 1) {
        for (i in 2:n_methods) {
          axis_name <- paste0("yaxis", i, ".fixedrange")
          plotly_layout_args[[axis_name]] <- is_on
        }
      }
      
      if (is_on) {
        # --- ENTERING VECTOR MODE ---
        shiny::updateActionButton(session, "toggle_mode_button", label = "Exit Vector Field Mode")
        shinyjs::show("vector_mode_controls", anim = TRUE)
        shinyjs::show("trajectory_selector_panel", anim = TRUE)
        
        # --- MODIFIED JAVASCRIPT ---
        # Forcefully disable scroll wheel zoom by setting the map's internal option as well.
        shinyjs::runjs("
        var map = $('#map').data('leaflet-map'); 
        if (map) { 
          map.dragging.disable(); 
          map.touchZoom.disable(); 
          map.doubleClickZoom.disable(); 
          
          map.options.scrollWheelZoom = false; 
          map.scrollWheelZoom.disable(); 
          
          $('.leaflet-control-zoom').hide(); 
        }
      ")
        plotly::plotlyProxy("subplots_combined", session) %>% plotly::plotlyProxyInvoke("relayout", plotly_layout_args)
        
        # --- PRE-COMPUTATION WITH PROGRESS BAR ---
        visible_ids <- visible_error_points()$time_id
        num_trajectories <- length(method_names) + n_models
        
        if (length(visible_ids) > 0 && num_trajectories > 0) {
          progress <- shiny::Progress$new(session, min = 0, max = length(visible_ids))
          progress$set(message = "Pre-computing all vector fields...", value = 0)
          on.exit(progress$close())
          
          temp_vector_cache <- list()
          
          # Loop through each visible time step
          for (time_idx in 1:length(visible_ids)) {
            id <- visible_ids[time_idx]
            progress$set(value = time_idx, detail = paste("Computing base fields for time point", time_idx))
            
            # STAGE 1: Compute the base fields ONCE for this time_id
            base_fields_data <- computeBaseFieldsForTime(id)
            
            # STAGE 2: Loop through each trajectory/method to apply weights
            for (traj_idx in 1:num_trajectories) {
              traj_name <- vis_method_names[traj_idx]
              
              # The first time we see a trajectory, initialize its list
              if (is.null(temp_vector_cache[[traj_name]])) {
                temp_vector_cache[[traj_name]] <- list()
              }
              
              # Pass the pre-computed base fields to the simplified calculation function
              temp_vector_cache[[traj_name]][[as.character(id)]] <- calculateVectorField(base_fields_data, traj_idx, vis_traj_list)
            }
          }
          
          precomputed_vectors(temp_vector_cache)
          
          # (Set initial selections code remains the same)
          selected_trajectory_for_field(all_vis_method_names[1])
          selected_time_id(visible_ids[1])
        }
        
      } else {
        # --- EXITING VECTOR MODE ---
        shiny::updateActionButton(session, "toggle_mode_button", label = "Enter Vector Field Mode")
        shinyjs::hide("vector_mode_controls", anim = TRUE)
        shinyjs::hide("trajectory_selector_panel", anim = TRUE) # Hide the selector
        
        # Clear the cache and the map layer
        precomputed_vectors(list())
        selected_trajectory_for_field(NULL)
        leaflet::leafletProxy("map", session) %>% 
          leaflet::clearGroup("vector_field_layer") %>%
          leaflet::clearGroup("intermediate_points_layer")
        
        # --- MODIFIED JAVASCRIPT ---
        # Re-enable scroll wheel zoom by setting the option back to true.
        shinyjs::runjs("
        var map = $('#map').data('leaflet-map'); 
        if (map) { 
          map.dragging.enable(); 
          map.touchZoom.enable(); 
          map.doubleClickZoom.enable(); 
          
          map.options.scrollWheelZoom = true;
          map.scrollWheelZoom.enable(); 
          
          $('.leaflet-control-zoom').show(); 
        }
      ")
        
        plotly::plotlyProxy("subplots_combined", session) %>% plotly::plotlyProxyInvoke("relayout", plotly_layout_args)
        
      }
    })
    
    shiny::observeEvent(input$next_button, {
      shiny::req(vector_mode_on()) # Only run in vector mode
      
      visible_ids <- visible_error_points()$time_id
      current_id <- selected_time_id()
      
      current_index <- which(visible_ids == current_id)
      
      # If we are not at the last point, select the next one
      if (length(current_index) > 0 && current_index < length(visible_ids)) {
        selected_time_id(visible_ids[current_index + 1])
      }
    })
    
    # Observer for the "Back" button click
    shiny::observeEvent(input$back_button, {
      shiny::req(vector_mode_on()) # Only run in vector mode
      
      visible_ids <- visible_error_points()$time_id
      current_id <- selected_time_id()
      
      current_index <- which(visible_ids == current_id)
      
      # If we are not at the first point, select the previous one
      if (length(current_index) > 0 && current_index > 1) {
        selected_time_id(visible_ids[current_index - 1])
      }
    })
    
    shiny::observeEvent(input$trajectory_selector, {
      selected_trajectory_for_field(input$trajectory_selector)
    })
    
    shiny::observeEvent(list(selected_time_id(), selected_trajectory_for_field()), {
      shiny::req(vector_mode_on())
      
      current_id <- selected_time_id()
      current_traj <- selected_trajectory_for_field()
      
      # When the selection is cleared (e.g., exiting mode), do nothing here.
      if (is.null(current_id) || is.null(current_traj)) return()
      
      # --- This block now only draws the "heavy" items ---
      
      # 1. Update the main blue highlight and the static plot lines
      updateMapAndPlotHighlights(current_id, hovered_time()) # We pass hover_time() so the red line doesn't disappear on click
      
      # 2. Draw the intermediate points
      drawIntermediatePoints(current_id)
      
      # 3. Draw the vector field
      drawVectorField(current_id, current_traj)
      
    }, ignoreNULL = FALSE)   
    
    # Observer to enable/disable the Next/Back buttons
    shiny::observe({
      shiny::req(vector_mode_on()) # Only run in vector mode
      
      visible_ids <- visible_error_points()$time_id
      current_id <- selected_time_id()
      
      # If there are no points or no selection, disable both
      if (length(visible_ids) == 0 || is.null(current_id)) {
        shinyjs::disable("next_button")
        shinyjs::disable("back_button")
        return()
      }
      
      current_index <- which(visible_ids == current_id)
      
      # Enable/disable "Back"
      if (current_index == 1) {
        shinyjs::disable("back_button")
      } else {
        shinyjs::enable("back_button")
      }
      
      # Enable/disable "Next"
      if (current_index == length(visible_ids)) {
        shinyjs::disable("next_button")
      } else {
        shinyjs::enable("next_button")
      }
    })
    
  }
  
  # 7. Run the application
  shiny::shinyApp(ui, server)
  
  
  
  
}
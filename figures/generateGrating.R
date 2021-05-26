generateGrating = function(spatialFreq = 5, angle = 0, phase = 0, fixation = TRUE, pctg_radius_in = 0.05, gaussianEnvelope = TRUE, sigma = 3, export = ".eps"){
  
  if (!("ggplot2" %in% .packages())){
    library(ggplot2)
    
  }
  
  # Create x axis for generating signal
  range = seq(-8,8,0.05)
  
  # White noise
  if (spatialFreq == 0){
    
    # Number of rows
    cols = length(range)
    rows = round(0.7*length(range))
    
    # Generate random noise
    white_noise = matrix(runif(cols*rows,min = -1, 1), nrow = rows, ncol = cols)
    
    # Find center of image
    midx = ceiling((nrow(white_noise)+1)/2)
    
    midy = ceiling((ncol(white_noise)+1)/2)
    
    # Add inner circles into the matrix with the grating image
    grating_mask = white_noise
    
    grating_mask = EBImage::gblur(grating_mask,1)
    
    # Create coordinates for plot and to generate outer and inner circles
    grid_coordinates = expand.grid(1:nrow(grating_mask),1:ncol(grating_mask))
    
    # Estimate distance from center to create outer and inner circles
    distance_from_center_grid = sqrt((grid_coordinates$Var1 - midx)^2 + (grid_coordinates$Var2 - midy)^2)
    
    # Find radius of outer and inner circles based on proportion of plot
    # radius_out = 0.4*max(distance_from_center_grid)
    
    radius_in = pctg_radius_in*max(distance_from_center_grid)
    
    # grating_mask = ifelse(circle_out == 1, 0, grating_rotated)
    
    # create outer and inner circles
    circle_fill_in = ifelse(((grid_coordinates$Var1 - midx)^2 + (grid_coordinates$Var2 - midy)^2) < radius_in^2, 1, 0)
    
    circle_in = matrix(circle_fill_in, nrow(grating_mask),ncol(grating_mask))
    
    if (fixation == TRUE){
      
      grating_mask = ifelse(circle_in == 1, 0, grating_mask)
      
    }
    
    # Fixation dot position
    dot.frame = data.frame(x = midx, y = midy, grating = 0)
    
    # Apply Gaussian envelope
    if (gaussianEnvelope == TRUE){
      
      x_axis = range
      
      y_axis = seq(x_axis[round(length(x_axis)*.15)],x_axis[floor(length(x_axis)*.85)],0.05)
      
      grid_axis = expand.grid(x_axis,y_axis)
      
      densities = unlist(purrr::map2(.x = grid_axis$Var1, .y = grid_axis$Var2, .f = ~mvtnorm::dmvnorm(c(.x,.y), mean = c(0,0), sigma = diag(rep(sigma,2)))))
      
      matrix_density = matrix(densities, nrow = length(y_axis), ncol = length(x_axis), byrow = TRUE)
      
      grating_mask = grating_mask*matrix_density
      
    }
    
    # Create data.frame with all data for plotting
    grating.frame = data.frame(x = grid_coordinates$Var1, y = grid_coordinates$Var2, grating = as.vector(grating_mask))
    
    # Generate plot
    plot_grating = ggplot(grating.frame, aes(x = x, y = y, fill = grating),show.legend = FALSE) + 
      geom_raster(show.legend = FALSE) +
      scale_fill_gradient2(low = "black",mid = "gray",high = "white") + 
      coord_flip() + 
      theme_void() + 
      theme(panel.border = element_blank())
    
    if (fixation == TRUE){
      
      plot_grating = plot_grating + geom_point(data = dot.frame, aes(x = x, y = y), shape = 16, size = 2, show.legend = FALSE)
      
    }
    # Save plot
    ggplot2::ggsave(filename = paste("grating_noise",export,sep = ""),plot = plot_grating, width = 5, height = 3.5, dpi = 1000)
    
    # Gabor patch    
  }else{
    
    # Create signal with specified frequency and phase
    signal = sin(range*spatialFreq + phase)
    
    # Replicate signal over many rows
    grating = rep(signal,times = round(0.7*length(signal)))
    
    # Reorganize data into a matrix
    grating = matrix(grating,nrow = round(0.7*length(signal)),ncol = (length(signal)), byrow = TRUE)
    
    # Find center of image
    midx = ceiling((nrow(grating)+1)/2)
    
    midy = ceiling((ncol(grating)+1)/2)
    
    # Transform angle from degrees to radians
    angle_change = angle - 90
    angle_rad = (angle_change * pi) / (180)
    
    # Create empty matrix to store rotated image
    grating_rotated = matrix(0,nrow = nrow(grating), ncol = ncol(grating))
    
    # Rotate grating
    for (i in 1:nrow(grating_rotated)){
      
      for (j in 1:ncol(grating_rotated)){
        
        # x= (i-midx)*cos(angle_rad)+(j-midy)*sin(angle_rad)
        x= (i-midx)*cos(angle_rad) - (j-midy)*sin(angle_rad);
        
        y= (i-midx)*sin(angle_rad) + (j-midy)*cos(angle_rad)
        
        x=round(x)+midx
        
        y=round(y)+midy
        
        if ((x>=1 && y>=1 && x<=nrow(grating) && y<=ncol(grating))){
          
          grating_rotated[i,j]=grating[x,y]; # k degrees rotated image  
          
        }
        
      }
    }
    
    # Create coordinates for plot and to generate outer and inner circles
    grid_coordinates = expand.grid(1:nrow(grating),1:ncol(grating))
    
    # Estimate distance from center to create outer and inner circles
    distance_from_center_grid = sqrt((grid_coordinates$Var1 - midx)^2 + (grid_coordinates$Var2 - midy)^2)
    
    # Find radius of outer and inner circles based on proportion of plot
    # radius_out = 0.4*max(distance_from_center_grid)
    
    radius_in = pctg_radius_in*max(distance_from_center_grid)
    
    # # create outer and inner circles
    # circle_fill = ifelse(sqrt((grid_coordinates$Var1 - midx)^2 + (grid_coordinates$Var2 - midy)^2) >= radius_out, 1, 0)
    # 
    # circle_out = matrix(circle_fill,nrow(grating),ncol(grating), byrow = FALSE)
    
    circle_fill_in = ifelse(((grid_coordinates$Var1 - midx)^2 + (grid_coordinates$Var2 - midy)^2) < radius_in^2, 1, 0)
    
    circle_in = matrix(circle_fill_in, nrow(grating),ncol(grating))
    
    
    # Add outer and inner circles into the matrix with the grating image
    grating_mask = grating_rotated
    
    grating_mask = EBImage::gblur(grating_mask,2)
    
    # grating_mask = ifelse(circle_out == 1, 0, grating_rotated)
    
    if (fixation == TRUE){
      
      grating_mask = ifelse(circle_in == 1, 0, grating_mask)
      
      
      
    }
    
    
    # Fixation dot position
    dot.frame = data.frame(x = midx, y = midy, grating = 0)
    
    # Apply Gaussian envelope
    if (gaussianEnvelope == TRUE){
      
      x_axis = range
      
      y_axis = seq(x_axis[round(length(x_axis)*.15)],x_axis[floor(length(x_axis)*.85)],0.05)
      
      grid_axis = expand.grid(x_axis,y_axis)
      
      densities = unlist(purrr::map2(.x = grid_axis$Var1, .y = grid_axis$Var2, .f = ~mvtnorm::dmvnorm(c(.x,.y), mean = c(0,0), sigma = diag(rep(sigma,2)))))
      
      matrix_density = matrix(densities, nrow = length(y_axis), ncol = length(x_axis), byrow = TRUE)
      
      grating_mask = grating_mask*matrix_density
      
    }
    
    # Create data.frame with all data for plotting
    grating.frame = data.frame(x = grid_coordinates$Var1, y = grid_coordinates$Var2, grating = as.vector(grating_mask))
    
    # Generate plot
    plot_grating = ggplot(grating.frame, aes(x = x, y = y, fill = grating),show.legend = FALSE) + 
      geom_raster(show.legend = FALSE) +
      scale_fill_gradient2(low = "black",mid = "gray",high = "white") + 
      coord_flip() + 
      theme_void() + 
      theme(panel.border = element_blank())
    
    if (fixation == TRUE){
      
      plot_grating = plot_grating + geom_point(data = dot.frame, aes(x = x, y = y), shape = 16, size = 2, show.legend = FALSE)
      
    }
    # Save plot
    ggplot2::ggsave(filename = paste("grating_",angle,export,sep = ""),plot = plot_grating, width = 5, height = 3.5, dpi = 1000)
  }
}
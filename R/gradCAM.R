#' Generate Grad-CAM heatmap for a CNN model
#'
#' @param model A trained luz model
#' @param img_path Path to the PNG image
#' @param target_class Target class index (0-based) for which to compute Grad-CAM
#' @param target_layer_name Name of the layer to visualize (default: last conv layer)
#'
#' @return A list with the original image, heatmap, and overlaid visualization
#' @export
generate_gradcam <- function(model, img_path, target_class = NULL, target_layer_name = NULL) {
  require(torch)
  require(png)
  require(grid)
  
  # Load image
  img <- png::readPNG(img_path)
  
  # Ensure image has 3 channels (RGB)
  if (length(dim(img)) == 2) {
    # Grayscale image - convert to RGB by replicating
    img <- array(img, dim = c(dim(img), 3))
    img[,,2] <- img[,,1]
    img[,,3] <- img[,,1]
  } else if (dim(img)[3] == 4) {
    # RGBA image - drop alpha channel
    img <- img[,,1:3]
  }
  
  # Convert to tensor and add batch dimension
  # img is [H, W, C], need [1, C, H, W]
  img_tensor <- torch_tensor(img)$permute(c(3, 1, 2))$unsqueeze(1)
  img_tensor$requires_grad_(TRUE)
  
  # Get the model's features module
  features <- model$model$features
  classifier <- model$model$classifier
  
  # Variables to store activations and gradients
  conv_output <- NULL
  activations_grad <- NULL
  
  # Forward pass to get conv output BEFORE adaptive pooling
  # The model has: conv layers, then adaptive_avg_pool2d, then dropout
  x <- img_tensor
  
  # Process through all layers except the last 2 (adaptive pool and dropout)
  for (i in 1:(length(features) - 2)) {
    x <- features[[i]](x)
  }
  
  # This is the output of the last conv layer with spatial dimensions intact
  conv_output <- x$clone()
  x$retain_grad()
  activations_grad <- x
  
  # Now apply adaptive pooling and dropout
  x <- features[[length(features) - 1]](x)  # adaptive_avg_pool2d
  x <- features[[length(features)]](x)      # dropout
  
  # Complete forward pass through classifier
  output <- x$squeeze()
  output <- classifier(output)
  
  # If target_class not specified, use the predicted class
  if (is.null(target_class)) {
    target_class <- torch_argmax(output)$item()
  }
  
  # Get the score for target class
  score <- output[target_class + 1]  # +1 for R's 1-based indexing
  
  # Backward pass to get gradients
  score$backward()
  
  # Get gradients of the target class w.r.t. feature maps
  gradients <- activations_grad$grad
  
  # Global average pooling of gradients
  weights <- gradients$mean(dim = c(3, 4), keepdim = TRUE)
  
  # Weighted combination of feature maps
  cam <- (conv_output * weights)$sum(dim = 2)
  cam <- torch_relu(cam)  # Apply ReLU to focus on positive influence
  
  # Normalize to 0-1
  cam <- cam - cam$min()
  cam <- cam / (cam$max() + 1e-8)
  
  # Convert to R array and resize to original image size
  cam_array <- as.array(cam$squeeze()$cpu())
  
  # Ensure cam_array is 2D
  if (length(dim(cam_array)) == 0) {
    cam_array <- matrix(cam_array, nrow = 1, ncol = 1)
  } else if (length(dim(cam_array)) == 1) {
    # Make it a square-ish matrix
    n <- length(cam_array)
    side <- ceiling(sqrt(n))
    cam_array <- matrix(c(cam_array, rep(0, side^2 - n)), nrow = side, ncol = side)
  }
  
  # Print dimensions for debugging
  message(sprintf("Heatmap dimensions: %d x %d", nrow(cam_array), ncol(cam_array)))
  
  # Resize heatmap to match input image dimensions
  cam_resized <- resize_heatmap(cam_array, dim(img)[1:2])
  
  # Create overlay
  overlay <- create_overlay(img, cam_resized)
  
  return(list(
    original = img,
    heatmap = cam_resized,
    overlay = overlay,
    predicted_class = target_class,
    heatmap_raw_dims = dim(cam_array)
  ))
}

#' Resize heatmap to target dimensions using bilinear interpolation
resize_heatmap <- function(heatmap, target_dims) {
  # Ensure heatmap is a matrix
  if (!is.matrix(heatmap)) {
    heatmap <- as.matrix(heatmap)
  }
  
  # Handle edge cases
  if (nrow(heatmap) == 0 || ncol(heatmap) == 0) {
    return(matrix(0, nrow = target_dims[1], ncol = target_dims[2]))
  }
  
  # Use a more efficient interpolation with expand.grid
  h_old <- nrow(heatmap)
  w_old <- ncol(heatmap)
  h_new <- target_dims[1]
  w_new <- target_dims[2]
  
  # Create coordinate grids
  x_ratio <- (h_old - 1) / max(h_new - 1, 1)
  y_ratio <- (w_old - 1) / max(w_new - 1, 1)
  
  interp_result <- matrix(0, nrow = h_new, ncol = w_new)
  
  for (i in 1:h_new) {
    for (j in 1:w_new) {
      # Calculate source coordinates
      x <- (i - 1) * x_ratio + 1
      y <- (j - 1) * y_ratio + 1
      
      # Get integer parts
      x1 <- floor(x)
      x2 <- min(ceiling(x), h_old)
      y1 <- floor(y)
      y2 <- min(ceiling(y), w_old)
      
      # Get fractional parts
      x_frac <- x - x1
      y_frac <- y - y1
      
      # Bilinear interpolation
      if (x1 == x2 && y1 == y2) {
        interp_result[i, j] <- heatmap[x1, y1]
      } else if (x1 == x2) {
        interp_result[i, j] <- heatmap[x1, y1] * (1 - y_frac) + heatmap[x1, y2] * y_frac
      } else if (y1 == y2) {
        interp_result[i, j] <- heatmap[x1, y1] * (1 - x_frac) + heatmap[x2, y1] * x_frac
      } else {
        interp_result[i, j] <- 
          heatmap[x1, y1] * (1 - x_frac) * (1 - y_frac) +
          heatmap[x2, y1] * x_frac * (1 - y_frac) +
          heatmap[x1, y2] * (1 - x_frac) * y_frac +
          heatmap[x2, y2] * x_frac * y_frac
      }
    }
  }
  
  return(interp_result)
}

#' Create overlay of heatmap on original image
create_overlay <- function(img, heatmap, alpha = 0.4) {
  # Convert heatmap to RGB using a colormap (jet-like)
  heatmap_rgb <- array(0, dim = c(dim(heatmap), 3))
  
  # Simple jet colormap approximation (blue -> cyan -> green -> yellow -> red)
  for (i in 1:nrow(heatmap)) {
    for (j in 1:ncol(heatmap)) {
      val <- heatmap[i, j]
      if (val < 0.25) {
        # Blue to Cyan
        heatmap_rgb[i, j, 3] <- 1
        heatmap_rgb[i, j, 2] <- val * 4
      } else if (val < 0.5) {
        # Cyan to Green
        heatmap_rgb[i, j, 3] <- 1 - (val - 0.25) * 4
        heatmap_rgb[i, j, 2] <- 1
      } else if (val < 0.75) {
        # Green to Yellow
        heatmap_rgb[i, j, 2] <- 1
        heatmap_rgb[i, j, 1] <- (val - 0.5) * 4
      } else {
        # Yellow to Red
        heatmap_rgb[i, j, 1] <- 1
        heatmap_rgb[i, j, 2] <- 1 - (val - 0.75) * 4
      }
    }
  }
  
  # Blend with original image
  overlay <- img * (1 - alpha) + heatmap_rgb * alpha
  overlay <- pmin(pmax(overlay, 0), 1)
  
  return(overlay)
}

#' Plot Grad-CAM visualization
#'
#' @param gradcam_result Output from generate_gradcam()
#' @param save_path Optional path to save the plot
#' @export
plot_gradcam <- function(gradcam_result, save_path = NULL) {
  if (!is.null(save_path)) {
    png(save_path, width = 1600, height = 533, res = 100)
  }
  
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
  
  # Original image
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = "Original Image")
  rasterImage(gradcam_result$original, 0, 0, 1, 1, interpolate = FALSE)
  
  # Heatmap
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Grad-CAM Heatmap (%dx%d)", 
                      gradcam_result$heatmap_raw_dims[1], 
                      gradcam_result$heatmap_raw_dims[2]))
  heatmap_rgb <- array(0, dim = c(dim(gradcam_result$heatmap), 3))
  # Use jet colormap for heatmap visualization
  for (i in 1:nrow(gradcam_result$heatmap)) {
    for (j in 1:ncol(gradcam_result$heatmap)) {
      val <- gradcam_result$heatmap[i, j]
      if (val < 0.25) {
        heatmap_rgb[i, j, 3] <- 1
        heatmap_rgb[i, j, 2] <- val * 4
      } else if (val < 0.5) {
        heatmap_rgb[i, j, 3] <- 1 - (val - 0.25) * 4
        heatmap_rgb[i, j, 2] <- 1
      } else if (val < 0.75) {
        heatmap_rgb[i, j, 2] <- 1
        heatmap_rgb[i, j, 1] <- (val - 0.5) * 4
      } else {
        heatmap_rgb[i, j, 1] <- 1
        heatmap_rgb[i, j, 2] <- 1 - (val - 0.75) * 4
      }
    }
  }
  rasterImage(heatmap_rgb, 0, 0, 1, 1, interpolate = FALSE)
  
  # Overlay
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Overlay (Predicted Class: %d)", gradcam_result$predicted_class))
  rasterImage(gradcam_result$overlay, 0, 0, 1, 1, interpolate = FALSE)
  
  if (!is.null(save_path)) {
    dev.off()
  }
}
#' Generate Grad-CAM by manually capturing activations
#'
#' @param model A trained luz model
#' @param img_path Path to the PNG image
#' @param target_class Target class index (0-based) for which to compute Grad-CAM
#'
#' @return A list with the original image, heatmap, and overlaid visualization
#' @export
generate_gradcam_hooks <- function(model, img_path, target_class = NULL) {
  require(torch)
  require(png)
  
  # Load and preprocess image
  img <- png::readPNG(img_path)
  
  # Ensure image has 3 channels (RGB)
  if (length(dim(img)) == 2) {
    img <- array(img, dim = c(dim(img), 3))
    img[,,2] <- img[,,1]
    img[,,3] <- img[,,1]
  } else if (dim(img)[3] == 4) {
    img <- img[,,1:3]
  }
  
  # Convert to tensor [1, C, H, W]
  img_tensor <- torch_tensor(img)$permute(c(3, 1, 2))$unsqueeze(1)
  img_tensor$requires_grad_(TRUE)
  
  # Get model components
  features <- model$model$features
  classifier <- model$model$classifier
  
  # Forward pass through layers up to and including the 5th conv2d (layer 17)
  x <- img_tensor
  for (i in 1:17) {
    x <- features[[i]](x)
  }
  
  # Save the activations from the 5th conv layer
  activations <- x$clone()
  x$retain_grad()
  saved_x <- x
  
  # Continue forward pass through remaining feature layers
  for (i in 18:length(features)) {
    x <- features[[i]](x)
  }
  
  # Forward through classifier
  x <- x$squeeze()
  output <- classifier(x)
  
  # Get predicted class if not specified
  if (is.null(target_class)) {
    target_class <- torch_argmax(output)$item()
  }
  
  message(sprintf("Target class: %d", target_class))
  
  # Backward pass for target class
  class_score <- output[target_class + 1]  # +1 for R's 1-based indexing
  class_score$backward()
  
  # Get gradients from the saved layer
  gradients <- saved_x$grad
  
  # Check if we captured activations and gradients
  if (is.null(activations) || is.null(gradients)) {
    stop("Failed to capture activations or gradients")
  }
  
  message(sprintf("Activation shape: [%s]", paste(dim(activations), collapse=", ")))
  message(sprintf("Gradient shape: [%s]", paste(dim(gradients), collapse=", ")))
  
  # Compute Grad-CAM
  # weights = global average pooling of gradients
  weights <- gradients$mean(dim = c(3, 4), keepdim = TRUE)
  
  # Weighted combination of activation maps
  cam <- (activations * weights)$sum(dim = 2)  # Sum over channel dimension
  
  # Apply ReLU to focus on positive contributions
  cam <- torch_relu(cam)
  
  # Normalize to [0, 1]
  cam <- cam - cam$min()
  cam <- cam / (cam$max() + 1e-8)
  
  # Convert to R array
  cam_array <- as.array(cam$squeeze()$cpu())
  
  message(sprintf("CAM dimensions: %d x %d", nrow(cam_array), ncol(cam_array)))
  
  # Resize to match image dimensions
  cam_resized <- resize_cam_bilinear(cam_array, c(dim(img)[1], dim(img)[2]))
  
  # Create overlay
  overlay <- create_cam_overlay(img, cam_resized, alpha = 0.5)
  
  return(list(
    original = img,
    heatmap = cam_resized,
    overlay = overlay,
    predicted_class = target_class,
    cam_raw = cam_array,
    cam_raw_dims = dim(cam_array)
  ))
}

#' Bilinear interpolation resize for CAM
resize_cam_bilinear <- function(cam, target_size) {
  h_old <- nrow(cam)
  w_old <- ncol(cam)
  h_new <- target_size[1]
  w_new <- target_size[2]
  
  # Create output matrix
  cam_resized <- matrix(0, nrow = h_new, ncol = w_new)
  
  # Compute scaling factors
  scale_h <- h_old / h_new
  scale_w <- w_old / w_new
  
  for (i in 1:h_new) {
    for (j in 1:w_new) {
      # Map output coordinates to input coordinates
      src_y <- (i - 0.5) * scale_h + 0.5
      src_x <- (j - 0.5) * scale_w + 0.5
      
      # Clamp to valid range
      src_y <- max(1, min(h_old, src_y))
      src_x <- max(1, min(w_old, src_x))
      
      # Get surrounding pixels
      y1 <- floor(src_y)
      y2 <- min(ceiling(src_y), h_old)
      x1 <- floor(src_x)
      x2 <- min(ceiling(src_x), w_old)
      
      # Compute interpolation weights
      wy <- src_y - y1
      wx <- src_x - x1
      
      # Bilinear interpolation
      if (y1 == y2 && x1 == x2) {
        cam_resized[i, j] <- cam[y1, x1]
      } else if (y1 == y2) {
        cam_resized[i, j] <- cam[y1, x1] * (1 - wx) + cam[y1, x2] * wx
      } else if (x1 == x2) {
        cam_resized[i, j] <- cam[y1, x1] * (1 - wy) + cam[y2, x1] * wy
      } else {
        cam_resized[i, j] <- 
          cam[y1, x1] * (1 - wx) * (1 - wy) +
          cam[y1, x2] * wx * (1 - wy) +
          cam[y2, x1] * (1 - wx) * wy +
          cam[y2, x2] * wx * wy
      }
    }
  }
  
  return(cam_resized)
}

#' Create overlay with jet colormap
create_cam_overlay <- function(img, cam, alpha = 0.5) {
  # Apply jet colormap to CAM
  cam_colored <- apply_jet_colormap(cam)
  
  # Blend with original image
  overlay <- img * (1 - alpha) + cam_colored * alpha
  overlay <- pmin(pmax(overlay, 0), 1)
  
  return(overlay)
}

#' Apply jet colormap to a matrix
apply_jet_colormap <- function(values) {
  # Create RGB array
  rgb_array <- array(0, dim = c(dim(values), 3))
  
  for (i in 1:nrow(values)) {
    for (j in 1:ncol(values)) {
      val <- values[i, j]
      
      # Jet colormap: blue -> cyan -> green -> yellow -> red
      if (val < 0.25) {
        rgb_array[i, j, 3] <- 1
        rgb_array[i, j, 2] <- val * 4
      } else if (val < 0.5) {
        rgb_array[i, j, 2] <- 1
        rgb_array[i, j, 3] <- 1 - (val - 0.25) * 4
      } else if (val < 0.75) {
        rgb_array[i, j, 2] <- 1
        rgb_array[i, j, 1] <- (val - 0.5) * 4
      } else {
        rgb_array[i, j, 1] <- 1
        rgb_array[i, j, 2] <- 1 - (val - 0.75) * 4
      }
    }
  }
  
  return(rgb_array)
}

#' Plot Grad-CAM visualization with hooks method
#'
#' @param gradcam_result Output from generate_gradcam_hooks()
#' @param save_path Optional path to save the plot
#' @export
plot_gradcam_hooks <- function(gradcam_result, save_path = NULL) {
  if (!is.null(save_path)) {
    png(save_path, width = 1600, height = 533, res = 100)
  }
  
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
  
  # Original image
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = "Original Image", asp = 1)
  rasterImage(gradcam_result$original, 0, 0, 1, 1, interpolate = FALSE)
  box()
  
  # Heatmap
  heatmap_colored <- apply_jet_colormap(gradcam_result$heatmap)
  
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Grad-CAM (%dx%d)", 
                      gradcam_result$cam_raw_dims[1],
                      gradcam_result$cam_raw_dims[2]), asp = 1)
  rasterImage(heatmap_colored, 0, 0, 1, 1, interpolate = TRUE)
  box()
  
  # Overlay
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Overlay (Class: %d)", gradcam_result$predicted_class), asp = 1)
  rasterImage(gradcam_result$overlay, 0, 0, 1, 1, interpolate = TRUE)
  box()
  
  if (!is.null(save_path)) {
    dev.off()
  }
}
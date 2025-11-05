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
  
  # Forward pass through features to get activations
  activations <- features(img_tensor)
  
  # Store the last conv layer output (before adaptive pooling removes it)
  # We need to hook into the layer before adaptive pooling
  # Let's get the output before the last few layers
  
  # Re-run with hook to capture the last conv output
  conv_output <- NULL
  activations_grad <- NULL
  
  # Forward pass to get conv output (before adaptive pooling)
  x <- img_tensor
  for (i in 1:length(features)) {
    x <- features[[i]](x)
    # Capture output after last conv2d (before adaptive pooling)
    if (i == length(features) - 2) {  # Before adaptive pool and last dropout
      conv_output <- x$clone()
      x$retain_grad()
      activations_grad <- x
    }
  }
  
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
  
  # Resize heatmap to match input image dimensions
  cam_resized <- resize_heatmap(cam_array, dim(img)[1:2])
  
  # Create overlay
  overlay <- create_overlay(img, cam_resized)
  
  return(list(
    original = img,
    heatmap = cam_resized,
    overlay = overlay,
    predicted_class = target_class
  ))
}

#' Resize heatmap to target dimensions
resize_heatmap <- function(heatmap, target_dims) {
  # Ensure heatmap is a matrix
  if (!is.matrix(heatmap)) {
    heatmap <- as.matrix(heatmap)
  }
  
  # Handle edge cases
  if (nrow(heatmap) == 0 || ncol(heatmap) == 0) {
    return(matrix(0, nrow = target_dims[1], ncol = target_dims[2]))
  }
  
  x_old <- seq(0, 1, length.out = nrow(heatmap))
  y_old <- seq(0, 1, length.out = ncol(heatmap))
  
  x_new <- seq(0, 1, length.out = target_dims[1])
  y_new <- seq(0, 1, length.out = target_dims[2])
  
  # Bilinear interpolation
  interp_result <- matrix(0, nrow = target_dims[1], ncol = target_dims[2])
  
  for (i in seq_along(x_new)) {
    for (j in seq_along(y_new)) {
      x_idx <- findInterval(x_new[i], x_old)
      y_idx <- findInterval(y_new[j], y_old)
      
      x_idx <- min(max(x_idx, 1), length(x_old) - 1)
      y_idx <- min(max(y_idx, 1), length(y_old) - 1)
      
      interp_result[i, j] <- heatmap[x_idx, y_idx]
    }
  }
  
  return(interp_result)
}

#' Create overlay of heatmap on original image
create_overlay <- function(img, heatmap, alpha = 0.4) {
  # Convert heatmap to RGB using a colormap (jet-like)
  heatmap_rgb <- array(0, dim = c(dim(heatmap), 3))
  
  # Simple jet colormap approximation
  for (i in 1:nrow(heatmap)) {
    for (j in 1:ncol(heatmap)) {
      val <- heatmap[i, j]
      if (val < 0.25) {
        heatmap_rgb[i, j, 3] <- 1
        heatmap_rgb[i, j, 1] <- val * 4
      } else if (val < 0.5) {
        heatmap_rgb[i, j, 3] <- 1 - (val - 0.25) * 4
        heatmap_rgb[i, j, 1] <- 1
      } else if (val < 0.75) {
        heatmap_rgb[i, j, 1] <- 1
        heatmap_rgb[i, j, 2] <- (val - 0.5) * 4
      } else {
        heatmap_rgb[i, j, 1] <- 1 - (val - 0.75) * 4
        heatmap_rgb[i, j, 2] <- 1
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
    png(save_path, width = 1200, height = 400)
  }
  
  par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
  
  # Original image
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Original")
  rasterImage(gradcam_result$original, 0, 0, 1, 1)
  
  # Heatmap
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Grad-CAM Heatmap")
  heatmap_rgb <- array(0, dim = c(dim(gradcam_result$heatmap), 3))
  heatmap_rgb[,,1] <- gradcam_result$heatmap
  rasterImage(heatmap_rgb, 0, 0, 1, 1)
  
  # Overlay
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Overlay (Class: %d)", gradcam_result$predicted_class))
  rasterImage(gradcam_result$overlay, 0, 0, 1, 1)
  
  if (!is.null(save_path)) {
    dev.off()
  }
}
#' Generate Grad-CAM heatmap for a CNN model
#'
#' @param model A trained luz model
#' @param img_path Path to the PNG image
#' @param target_class Target class index (0-based) for which to compute Grad-CAM
#'
#' @return A list with the original image, heatmap, and overlaid visualization
#' @export
generate_gradcam <- function(model, img_path, target_class = NULL) {
  require(torch)
  require(png)
  
  # Load image
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
  
  # Forward pass through features (stop before adaptive pooling)
  x <- img_tensor
  for (i in 1:(length(features) - 2)) {
    x <- features[[i]](x)
  }
  
  # Save activations from last conv layer
  conv_output <- x$clone()
  x$retain_grad()
  activations <- x
  
  # Complete forward pass
  x <- features[[length(features) - 1]](x)  # adaptive_avg_pool2d
  x <- features[[length(features)]](x)      # dropout
  output <- x$squeeze()
  output <- classifier(output)
  
  # Get target class
  if (is.null(target_class)) {
    target_class <- torch_argmax(output)$item()
  }
  
  # Backward pass
  score <- output[target_class + 1]
  score$backward()
  
  # Compute Grad-CAM
  gradients <- activations$grad
  weights <- gradients$mean(dim = c(3, 4), keepdim = TRUE)
  cam <- (conv_output * weights)$sum(dim = 2)
  cam <- torch_relu(cam)
  
  # Normalize
  cam <- cam - cam$min()
  cam <- cam / (cam$max() + 1e-8)
  
  # Convert to array
  cam_array <- as.array(cam$squeeze()$cpu())
  
  message(sprintf("Raw heatmap size: %d x %d", nrow(cam_array), ncol(cam_array)))
  
  # Resize to match image dimensions using simple nearest neighbor for clarity
  cam_resized <- resize_heatmap_simple(cam_array, dim(img)[1:2])
  
  # Create overlay
  overlay <- create_overlay_simple(img, cam_resized)
  
  return(list(
    original = img,
    heatmap = cam_resized,
    overlay = overlay,
    predicted_class = target_class,
    heatmap_raw_size = dim(cam_array)
  ))
}

#' Simple resize using nearest neighbor interpolation
resize_heatmap_simple <- function(heatmap, target_dims) {
  if (!is.matrix(heatmap)) {
    heatmap <- as.matrix(heatmap)
  }
  
  h_old <- nrow(heatmap)
  w_old <- ncol(heatmap)
  h_new <- target_dims[1]
  w_new <- target_dims[2]
  
  result <- matrix(0, nrow = h_new, ncol = w_new)
  
  # Calculate which source pixel each target pixel maps to
  for (i in 1:h_new) {
    for (j in 1:w_new) {
      # Simple mapping: divide target space into equal blocks
      src_i <- ceiling(i * h_old / h_new)
      src_j <- ceiling(j * w_old / w_new)
      
      # Clamp to valid range
      src_i <- min(max(src_i, 1), h_old)
      src_j <- min(max(src_j, 1), w_old)
      
      result[i, j] <- heatmap[src_i, src_j]
    }
  }
  
  return(result)
}

#' Create overlay with jet colormap
create_overlay_simple <- function(img, heatmap, alpha = 0.5) {
  # Create RGB heatmap with jet colormap
  heatmap_rgb <- array(0, dim = c(dim(heatmap), 3))
  
  for (i in 1:nrow(heatmap)) {
    for (j in 1:ncol(heatmap)) {
      val <- heatmap[i, j]
      
      # Jet colormap: blue -> cyan -> green -> yellow -> red
      if (val < 0.25) {
        # Blue to cyan
        heatmap_rgb[i, j, 3] <- 1
        heatmap_rgb[i, j, 2] <- val * 4
      } else if (val < 0.5) {
        # Cyan to green
        heatmap_rgb[i, j, 2] <- 1
        heatmap_rgb[i, j, 3] <- 1 - (val - 0.25) * 4
      } else if (val < 0.75) {
        # Green to yellow
        heatmap_rgb[i, j, 2] <- 1
        heatmap_rgb[i, j, 1] <- (val - 0.5) * 4
      } else {
        # Yellow to red
        heatmap_rgb[i, j, 1] <- 1
        heatmap_rgb[i, j, 2] <- 1 - (val - 0.75) * 4
      }
    }
  }
  
  # Blend with original
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
       main = "Original Image", asp = 1)
  rasterImage(gradcam_result$original, 0, 0, 1, 1, interpolate = FALSE)
  box()
  
  # Heatmap
  heatmap_vis <- gradcam_result$heatmap
  heatmap_rgb <- array(0, dim = c(dim(heatmap_vis), 3))
  
  for (i in 1:nrow(heatmap_vis)) {
    for (j in 1:ncol(heatmap_vis)) {
      val <- heatmap_vis[i, j]
      if (val < 0.25) {
        heatmap_rgb[i, j, 3] <- 1
        heatmap_rgb[i, j, 2] <- val * 4
      } else if (val < 0.5) {
        heatmap_rgb[i, j, 2] <- 1
        heatmap_rgb[i, j, 3] <- 1 - (val - 0.25) * 4
      } else if (val < 0.75) {
        heatmap_rgb[i, j, 2] <- 1
        heatmap_rgb[i, j, 1] <- (val - 0.5) * 4
      } else {
        heatmap_rgb[i, j, 1] <- 1
        heatmap_rgb[i, j, 2] <- 1 - (val - 0.75) * 4
      }
    }
  }
  
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Grad-CAM Heatmap (%dx%d)", 
                      gradcam_result$heatmap_raw_size[1],
                      gradcam_result$heatmap_raw_size[2]), asp = 1)
  rasterImage(heatmap_rgb, 0, 0, 1, 1, interpolate = FALSE)
  box()
  
  # Overlay
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       main = sprintf("Overlay (Class: %d)", gradcam_result$predicted_class), asp = 1)
  rasterImage(gradcam_result$overlay, 0, 0, 1, 1, interpolate = FALSE)
  box()
  
  if (!is.null(save_path)) {
    dev.off()
  }
}
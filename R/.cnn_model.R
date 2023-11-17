require(luz)
require(torch)
require(torchvision)


convnet_dropout <- nn_module(
  "convnet",
  initialize = function() {
    self$features <- nn_sequential(
      nn_conv2d(3, npx, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx, npx*2, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*2, npx*4, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*4, npx*8, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*8, npx*16, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_adaptive_avg_pool2d(c(1, 1)),
      nn_dropout2d(p = 0.05),
    )
    self$classifier <- nn_sequential(
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, classes)
    )
  },
  forward = function(x) {
    x <- self$features(x)$squeeze()
    x <- self$classifier(x)
    x
  }
)

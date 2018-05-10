#' Calculate a bootstrapped assessment of height specific RMSE and bias:
#'
#' @param x The reference data
#' @param y The predicted data
#' @param height_window The tree height window
#' @param n_boot The number of times to draw a random sample within the height window
#' @param n_sample The number of samples to draw at each iteration
#' @return A dataframe containing the estimates for each height window
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-09-01
#'


boot_height_rmse<-function(x, y, height_window, n_boot, n_sample){

  # remove incomplete cases !!!!
  xy_cc_ind<-complete.cases(cbind(x, y))
  x<-x[xy_cc_ind]
  y<-y[xy_cc_ind]

  x_min<-0
  x_max<-max(x, na.rm = TRUE)

  height_rmse <- lapply(seq(from = x_min, to = x_max, by = 0.1), function(i) {
    x_filter<-x[(x >= i &
                   x < (i + height_window)) ]
    y_filter<-y[(x >= i &
                   x < (i + height_window)) ]

    smplr <-
      sapply(1:n_boot, function(X)
        sample(1:length(x_filter), size = n_sample, replace = TRUE))

    # Run the bootstrap
    y_boot <- apply(smplr, 2, function(X) {
      x_smpl <-x_filter[as.numeric(X)]
      y_smpl <-y_filter[as.numeric(X)]
      m_height <-mean(x_smpl) # The mean true height for the sample
      rmse_smpl<-rmse(y = y_smpl, x = x_smpl)
      bias_smpl<-bias(y = y_smpl, x = x_smpl)

      rmse_smpl<- (rmse_smpl/m_height) *100 # convert RMSE to percentages
      bias_smpl<- (bias_smpl/m_height) *100# ... and again for bias

      output<-data.frame(height = m_height,
                         rmse = rmse_smpl,
                         bias = bias_smpl)

      return(output)
    })
    y_boot<-do.call(rbind, y_boot)

    m_height <-mean(y_boot$height)

    step_rmse <- mean(y_boot$rmse)
    lower_rmse <- quantile(y_boot$rmse, 0.05, na.rm = TRUE)
    upper_rmse <- quantile(y_boot$rmse, 0.95, na.rm = TRUE)

    step_bias <- mean(y_boot$bias)
    lower_bias <- quantile(y_boot$bias, 0.05, na.rm = TRUE)
    upper_bias <- quantile(y_boot$bias, 0.95, na.rm = TRUE)

    y <-
      data.frame(
        height = m_height,
        rmse = step_rmse,
        lower_rmse = lower_rmse,
        upper_rmse = upper_rmse,
        bias = step_bias,
        lower_bias = lower_bias,
        upper_bias = upper_bias
      )
    return(y)
  })
  height_rmse <- do.call(rbind, height_rmse)
  return(height_rmse)
}

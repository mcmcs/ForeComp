
#' Visualizes the size distortion maximum power loss tradeoff
#'
#' @description `Plot_Tradeoff` creates a plot to show sensitivity of statistical significance to the choice of bandwidth and how size distortion and maximum power loss vary.
#'
#' @param data Data frame
#' @param f1 Column containing forecaster 1's predictions. Should be a string.
#' @param f2 Column containing forecaster 2's predictions. Should be a string.
#' @param y Column containing the realized value for the outcome variable. Should be a string.
#' @param loss_function The transformation applied to the forecast error. Defaults to squared error loss. The user supplied function should take two inputs and a scalar output. For example, loss = loss_function(f, y).
#' @param n_sim The number of simulations used to generate the ARIMA model. Defaults to 1,000.
#' @param m_set The truncation parameter. Defaults to c(1:10, seq( 11, floor(nrow(data)/2), 10)). For a standard long-run variance calculation (for example, using Bartlett kernel), it controls the number of terms used in estimating the autocovariance matrix. It should be a vector of integers with the values of M you would like to plot.
#' @param no_m_label TRUE if plot without m labels. Defaults to FALSE.
#' @return A list of length 2. The first element is a ggplot2 object of the size-power tradeoff. The second element is the underlying data used to construct the plot in element 1.
#' @importFrom forecast auto.arima
#' @importFrom stats acf
#' @importFrom stats pnorm
#' @importFrom stats approx
#' @importFrom astsa arma.spec
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' \donttest{
#' ## A typical example
#' set.seed(1234);
#' output = Plot_Tradeoff(
#'   data = TBILL,
#'   f1   = "SPFfor_Step1",
#'   f2   = "NCfor_Step1",
#'   y    = "Realiz1"
#' )
#' output[[1]] # The first element is a ggplot2 object of the size-power tradeoff.
#' output[[2]] # The second element is the underlying data used to construct the plot in element 1.
#'
#'
#' ## An example with a user supplied M values (with a larger set of M values)
#' set.seed(1234);
#' Plot_Tradeoff(
#'   data = TBILL,
#'   f1 = "SPFfor_Step1",
#'   f2 = "NCfor_Step1",
#'   y  = "Realiz1",
#'   m_set = c(1:10, seq(from = 11, to = nrow(TBILL) - 20, by = 10))
#' )
#'
#' ## An example without (f1, f2, y). The function will take the first three columns and use them
#' set.seed(1234);
#' tmpdata = TBILL[, c("SPFfor_Step1", "NCfor_Step1", "Realiz1")]; # data with [f1, f2, y]
#' Plot_Tradeoff(
#'   data = tmpdata
#' )
#'
#' }



Plot_Tradeoff <- function(data,
                          f1 = NULL,
                          f2 = NULL,
                          y  = NULL,
                          loss_function = NULL,
                          n_sim = 1000,
                          m_set = NULL,
                          no_m_label = FALSE) {

  # Handling options

  # Data
  if (is.null(f1)& is.null(f2) & is.null(y)){
    f1 = data[[1]];
    f2 = data[[2]];
    y  = data[[3]];
  } else if (!is.null(f1) & !is.null(f2) & !is.null(y)) {
    f1 <- data[[f1]]
    f2 <- data[[f2]]
    y <- data[[y]]
  } else {
    stop("f1, f2, y have to be supplied altogether.")
  }

  # If the user does not supply a loss_function, we use a quadratic loss
  if (is.null(loss_function)) {
    loss_function = function(f, y){ return( (f-y)^2 ); };
  }

  # calculating loss
  loss1 = loss_function(f1, y);
  loss2 = loss_function(f2, y);
  d     = loss1 - loss2;

  # mset - default value
  if (is.null(m_set)){
    if ( floor(nrow(data)/2) > 10){
      m_set = c(1:10, seq( 11, floor(nrow(data)/2), 10));
    } else {
      m_set = seq(1, floor(nrow(data)/2), 1);
    }
  }

  if (!is.logical(no_m_label)) {stop("Argument 'no_m_label' should be either TRUE or FALSE.")}
  if (!is.numeric(n_sim) | n_sim %% 1 != 0 | n_sim <= 0) {stop("Argument 'n_sim' should be a natural number.")}


  # other info
  conf_level <- 0.05
  series_length <- nrow(data)
  m_set_length <- length(m_set)

  # matrix to store
  mat_size_distortion_dm <- matrix(NA, m_set_length, 1)
  mat_power_loss_dm <- matrix(NA, m_set_length, 1)

  mat_size_distortion_b <- matrix(NA, m_set_length, 1)
  mat_power_loss_b <- matrix(NA, m_set_length, 1)

  mat_size_distortion_d <- matrix(NA, m_set_length, 1)
  mat_power_loss_d <- matrix(NA, m_set_length, 1)

  mut <- mean(d); # demean
  dt_tilde <- d - mut


  # --- Estimate ARIMA and extract information
  a = forecast::auto.arima(y=dt_tilde, max.p = 12, max.q=12, stationary=T, ic="aic", seasonal=F, allowmean=F);
  i_ar = grep("ar", names(a$coef));
  i_ma = grep("ma", names(a$coef));
  m_sim = list("ar"=a$coef[i_ar], "ma"=a$coef[i_ma]);

  if (length(m_sim$ar)==0){m_sim$ar=0.0}
  if (length(m_sim$ma)==0){m_sim$ma=0.0}

  # --- Generate data for size and power computation (this loop alos pre-computes long-run variance)
  Mmax = min(c(series_length-1,200)); # maximum possible M considered in this experimen
  mat_dt   = matrix(NA, series_length, n_sim); # matrix that stores data (for size calculation)
  mat_dtm  = matrix(NA, 1, n_sim);    # matrix that stores mean of data
  mat_acf  = matrix(NA, (Mmax+1), n_sim); # matrix that stores acf
  for (irep in 1:n_sim){

    # data
    dt_sim = arima.sim(m_sim, n=series_length, innov = rnorm(series_length, 0, sqrt(a$sigma2)));
    mat_dt[,irep] = dt_sim;
    mat_dtm[,irep] = mean(dt_sim);
    # autocovariance matrix
    d.cov = stats::acf(dt_sim, lag.max = (Mmax), type="covariance", plot=FALSE, demean=TRUE)$acf[,,1]; #should I compute ACF under the null regardless? Then, this should be deman = FALSE
    mat_acf[,irep] = d.cov;
  }


  # --- Setting for trade-off figure
  del_grid = seq(from=0, to=10, by=0.25);
  ndel = length(del_grid);

  v_M <- vector(mode = "integer", length = m_set_length)
  v_hypothesis_test_b <- vector(mode = "logical", length = m_set_length)
  v_test_statistic_b <- vector(mode = "logical", length = m_set_length)
  v_hypothesis_test_dm <- vector(mode = "logical", length = m_set_length)
  v_test_statistic_dm <- vector(mode = "logical", length = m_set_length)

  # --- Loop over M set
  for (iM in 1:m_set_length){

    Mchoice = m_set[iM]; #our choice of M for this iteration

    testres_b  = dm.test.bt.fb(d, cl = conf_level, M = Mchoice);
    wce_b_rej  = testres_b$rej
    wce_b_stat = testres_b$stat

    testres_dm  = dm.test.bt(d, cl = conf_level, M = Mchoice);
    wce_dm_rej  = testres_dm$rej  # CHANGE
    wce_dm_stat = testres_dm$stat # CHANGE

    v_M[iM] <- Mchoice
    v_hypothesis_test_b[iM] <- wce_b_rej
    v_test_statistic_b[iM] <- wce_b_stat
    v_hypothesis_test_dm[iM] <- wce_dm_rej
    v_test_statistic_dm[iM] <- wce_dm_stat

    # --- Oracle ---

    # --- Power for oracle

    # long-run variance
    ss = astsa::arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100);
    Om = ss$spec[1]; #this is 2*pi*f(0), spectrum at zero rather than a spectral density at zero

    # --- Standard test


    # --- Size computation for DM-WCE-dm, DM-WCE-b, DM-WPE-d
    M = Mchoice;
    mat_stat =    matrix(NA, n_sim, 1);
    mat_d.var   = matrix(NA, n_sim, 1);
    mat_rej_dm  = matrix(NA, n_sim, 1);
    mat_rej_b   = matrix(NA, n_sim, 1);
    for (irep in 1:n_sim){
      d.cov  = mat_acf[1:(M+1),irep];
      d.var  = ( d.cov[1] + 2*sum( (1 - ((1:M)/M) ) * d.cov[-1] ) ) / series_length;

      # Traditional NW, DM (DM-WCE-dm)
      dmstat = mean(mat_dtm[,irep]) / sqrt(d.var);
      pval   = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
      rej_dm    = pval < conf_level; #reject decision

      # Fixed-b NW (DM-WCE-b)
      b = M / series_length;
      if (conf_level == 0.05){
        crit = 1.9600 + 2.9694*b + 0.4160*b^2 -0.5324*b^3; #0.975 quantile
      } else if (conf_level == 0.10){
        crit = 1.6449 + 2.1859*b + 0.3142*b^2 -0.3427*b^3; #0.950 quantile
      }
      # rejection decision
      rej_b = abs(dmstat) > crit; #reject decision


      # ***Minchul DM-WPE-d
      # do be done


      # store results
      mat_stat[irep,] = dmstat;
      mat_d.var[irep,] = d.var;
      mat_rej_dm[irep, ] = rej_dm;
      mat_rej_b[irep, ] = rej_b;
    }
    size_distortion_dm = mean(mat_rej_dm) - conf_level;
    size_distortion_b = mean(mat_rej_b) - conf_level;


    # --- Power of a standard test (DM-NW)
    # Note that WCE-DM and WCE-B have the same power property when size-corrected
    # --- size-corrected crit val
    c05_star = quantile(abs(mat_stat), (1-conf_level));

    # --- size-corrected power
    mat_stat = matrix(NA, n_sim, ndel);
    mat_rej  = matrix(NA, n_sim, ndel);
    mat_rej2 = matrix(NA, n_sim, ndel);
    for (irep in 1:n_sim){

      m_i = mat_dtm[,irep];
      s_i = sqrt(mat_d.var[irep,]);

      for (idel in 1:ndel){

        # DM test with NW
        dmstat = (m_i + (1/sqrt(series_length)) *sqrt(Om)*del_grid[idel]) / s_i; #statistic
        pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
        rej = pval < conf_level; #reject decision

        mat_rej[irep, idel] = rej; #for raw power
        mat_stat[irep, idel] = dmstat;
        mat_rej2[irep, idel] = abs(dmstat) > c05_star; # for size-corrected power
      }
    }


    #======================================

    # --- another Oracle ---
    grid = seq(from=0, to=5, by=0.05);
    samp = del_grid;
    pow_gau = stats::pnorm(-1.96+grid) + stats::pnorm(-1.96-grid);

    # --- Maximum power loss
    powinterp = stats::approx(samp, apply(mat_rej2, 2, mean), grid)$y; #interpolated power
    max_power_loss_dm = max(pow_gau-powinterp);
    max_power_loss_b = max_power_loss_dm; #WCE-DM and WCE-B have the same power property

    # --- Collect results
    print(paste0("M = ", iM, " / ", m_set_length));
    mat_size_distortion_dm[iM] = size_distortion_dm;
    mat_power_loss_dm[iM] = max_power_loss_dm;

    mat_size_distortion_b[iM] = size_distortion_b;
    mat_power_loss_b[iM] = max_power_loss_b;

  } #end of iM iteration

  df_hypoth_testing <- data.frame(
    v_M,
    v_hypothesis_test_b,
    v_test_statistic_b,
    v_hypothesis_test_dm,
    v_test_statistic_dm
  )

  df_size_power <- data.frame(
    M = m_set,
    b_size_distortion = mat_size_distortion_b[,1],
    b_power_loss = mat_power_loss_b[,1]
  )

  plotting_data <- merge(df_size_power, df_hypoth_testing, by.x = "M", by.y = "v_M")
  plotting_data["v_hypothesis_test_b"] <- ifelse(v_hypothesis_test_b == TRUE, "cross", "circle")



  plot <- ggplot(plotting_data, aes(x = b_size_distortion, y = b_power_loss)) +
    geom_path(linewidth = 1, linetype = "dashed") +
    geom_point(aes(shape = v_hypothesis_test_b), size = 4.5, color = "red", stroke = 1.5) +
    scale_shape_identity() +
    labs(
      x = "Size Distortion",
      y = "Maximum Power Loss"
    ) +
    theme_minimal()

  # add m_set annotation
  if (no_m_label == FALSE){
    plot <- plot + geom_text(aes(label = M), nudge_y = .005)
  }



  return(list(plot, plotting_data))
}


#' Visualizes the size distortion maximum power loss tradeoff
#'
#' @description `Plot_Tradeoff` creates a plot to show sensitivity of statistical significance to the choice of bandwidth and how size distortion and maximum power loss vary.
#'
#' @param data
#' @param f1 Column containing forecaster 1's predictions. Should be a string.
#' @param f2 Column containing forecaster 2's predictions. Should be a string.
#' @param y Column containing the realized value for the outcome variable. Should be a string.
#' @param loss_function The transformation applied to the forecast error. Defaults to squared error loss. Should be a string.
#' @param n_sim The number of simulations used to generate the ARIMA model. Defaults to 1,000. Should be a string.
#' @param m_set The truncation parameter that controls the number of terms used in estimating the autocovariance matrix. Defaults to M = c(1:10, seq(11, nrow(data) - 1, 10)). Should be a vector of integers with the values of M you would like to plot.
#'
#' @return A list of length 2. The first element is a ggplot2 object of the size-power tradeoff. The second element is the underlying data used to construct the plot in element 1.
#' @export
#'
#' @examples

Plot_Tradeoff <- function(data, f1 = NULL, f2 = NULL, y = NULL, loss_function = NULL, n_sim = 1000, m_set =  c(1:10, seq(11, nrow(data) - 1, 10))) {

  # ==================================================
  # conf_level is currently hard-coded, and set to 0.05
  conf_level = 0.05;
  # ==================================================

  # If the user does not supply a loss_function, we use a quadratic loss
  if (is.null(loss_function)) {
    loss_function = function(f, y){ return( (f-y)^2 ); };
  }

  # Housekeeping
  data          = d_t;            # dt is a loss differential series
  series_length = nrows(data);   # the number of observations
  M_set_n       = length(M_set); # length of M grid
  number_simulations = N_sim;     # number of simulations

  # matrix to store results
  # _dm = WCE-DM (traditional NW)
  # _b  = WCE-B  (fixed-b NW)
  # _d  = WPE-D  (fixed-b)
  mat_size_distortion_dm = matrix(NA, M_set_n, 1);
  mat_power_loss_dm = matrix(NA, M_set_n, 1);

  mat_size_distortion_b = matrix(NA, M_set_n, 1);
  mat_power_loss_b = matrix(NA, M_set_n, 1);

  mat_size_distortion_d = matrix(NA, M_set_n, 1);
  mat_power_loss_d = matrix(NA, M_set_n, 1);

  mut = mean(data); # demean
  dt_tilde = data - mut;


  # --- Estimate ARIMA and extract information
  a = auto.arima(y=dt_tilde, max.p = 12, max.q=12, stationary=T, ic="aic", seasonal=F, allowmean=F);
  i_ar = grep("ar", names(a$coef));
  i_ma = grep("ma", names(a$coef));
  m_sim = list("ar"=a$coef[i_ar], "ma"=a$coef[i_ma]);

  if (length(m_sim$ar)==0){m_sim$ar=0.0}
  if (length(m_sim$ma)==0){m_sim$ma=0.0}

  # --- Generate data for size and power computation (this loop alos pre-computes long-run variance)
  Mmax = min(c(series_length-1,200)); # maximum possible M considered in this experimen
  mat_dt   = matrix(NA, series_length, number_simulations); # matrix that stores data (for size calculation)
  mat_dtm  = matrix(NA, 1, number_simulations);    # matrix that stores mean of data
  mat_acf  = matrix(NA, (Mmax+1), number_simulations); # matrix that stores acf
  for (irep in 1:number_simulations){

    # data
    dt_sim = arima.sim(m_sim, n=series_length, innov = rnorm(series_length, 0, sqrt(a$sigma2)));
    mat_dt[,irep] = dt_sim;
    mat_dtm[,irep] = mean(dt_sim);
    # autocovariance matrix
    d.cov = stats::acf(dt_sim, lag.max = (Mmax), type="covariance", plot=FALSE, demean=TRUE)$acf[,,1]; #should I compute ACF under the null regardless? Then, this should be deman = FALSE
    mat_acf[,irep] = d.cov;
  }


  # --- Setting for trade-off figure
  # ndel = 50; # number of deltas
  # del_tilde = 10; # largest delta
  # del_grid = seq(from=-del_tilde, to=del_tilde, length.out=ndel);
  del_grid = seq(from=0, to=10, by=0.25);
  ndel = length(del_grid);

  v_M <- vector(mode = "integer", length = M_set_n)
  v_hypothesis_test_b <- vector(mode = "logical", length = M_set_n)
  v_test_statistic_b <- vector(mode = "logical", length = M_set_n)
  v_hypothesis_test_dm <- vector(mode = "logical", length = M_set_n)
  v_test_statistic_dm <- vector(mode = "logical", length = M_set_n)

  # --- Loop over M set
  for (iM in 1:M_set_n){

    Mchoice = M_set[iM]; #our choice of M for this iteration

    testres_b  = dm.test.bt.fb(data, conf_level = .05, M = Mchoice);
    wce_b_rej  = testres_b$rej
    wce_b_stat = testres_b$stat

    testres_dm  = dm.test.bt(data, conf_level = .05, M = Mchoice);
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
    ss = arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100);
    Om = ss$spec[1]; #this is 2*pi*f(0), spectrum at zero rather than a spectral density at zero

    # # Oracle test (we will do a simpler calculation, see below)
    # mat_stat_o = matrix(NA, nsim, 1);
    # mat_rej_o = matrix(NA, nsim, 1);
    # for (irep in 1:nsim){
    #   dmstat = mean(mat_dtm[,irep]) / sqrt(Om / nlen); #Oracle's statistic
    #   pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
    #   rej = pval < cl; #reject decision
    #   mat_stat_o[irep, 1] = dmstat;
    #   mat_rej_o[irep,1] = rej;
    # }
    # # print("empirical size of oracle");
    # # print(mean(mat_rej_o));
    #
    #
    # # --- find a cut-off (c*) to get size-corrected critical value
    # c05_star_o = quantile(abs(mat_stat_o), (1-cl));
    #
    # # size corrected power
    # mat_rej_o  = matrix(NA, nsim, ndel);
    # mat_rej2_o = matrix(NA, nsim, ndel);
    # for (idel in 1:ndel){
    #   for (irep in 1:nsim){
    #     dmstat = (mat_dtm[,irep] + (1/sqrt(nlen)) *sqrt(Om)*del_grid[idel]) / sqrt(Om / nlen); #statistic
    #     pval = 2*stats::pnorm(-abs(dmstat), mean=0, sd=1); # p-val from normal approximation
    #     rej = pval < cl; # reject decision
    #     mat_rej_o[irep, idel] = rej; # this is for raw power
    #     mat_rej2_o[irep, idel] = abs(dmstat) > c05_star_o; # for size-corrected power
    #   }
    # }



    # --- Standard test


    # --- Size computation for DM-WCE-dm, DM-WCE-b, DM-WPE-d
    M = Mchoice;
    mat_stat =    matrix(NA, number_simulations, 1);
    mat_d.var   = matrix(NA, number_simulations, 1);
    mat_rej_dm  = matrix(NA, number_simulations, 1);
    mat_rej_b   = matrix(NA, number_simulations, 1);
    for (irep in 1:number_simulations){
      d.cov  = mat_acf[1:(M+1),irep];
      d.var  = ( d.cov[1] + 2*sum( (1 - ((1:M)/M) ) * d.cov[-1] ) ) / series_length;

      # dt_sim = mat_dt[,irep];
      # d.cov = stats::acf(dt_sim, lag.max = (M), type="covariance", plot=FALSE, demean=FALSE)$acf[,,1]; #should I compute ACF under the null regardless? Then, this should be deman = FALSE
      # d.var  = ( d.cov[1] + 2*sum( (1 - ((1:M)/M) ) * d.cov[-1] ) ) / nlen;

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
    mat_stat = matrix(NA, number_simulations, ndel);
    mat_rej  = matrix(NA, number_simulations, ndel);
    mat_rej2 = matrix(NA, number_simulations, ndel);
    for (irep in 1:number_simulations){

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
    pow_gau = pnorm(-1.96+grid) + pnorm(-1.96-grid);
    # plot(grid, approx(samp, apply(mat_rej2, 2, mean), grid)$y)
    # lines(grid, pow_gau)


    # --- Maximum power loss
    powinterp = approx(samp, apply(mat_rej2, 2, mean), grid)$y; #interpolated power
    max_power_loss_dm = max(pow_gau-powinterp);
    max_power_loss_b = max_power_loss_dm; #WCE-DM and WCE-B have the same power property

    # max_power_loss = max( apply(mat_rej2_o, 2, mean) - apply(mat_rej2, 2, mean) );

    # --- Collect results
    print(paste0("M = ", iM, " / ", M_set_n));
    mat_size_distortion_dm[iM] = size_distortion_dm;
    mat_power_loss_dm[iM] = max_power_loss_dm;

    mat_size_distortion_b[iM] = size_distortion_b;
    mat_power_loss_b[iM] = max_power_loss_b;

  } #end of iM iteration

  df_hypoth_testing <- tibble(
    v_M,
    v_hypothesis_test_b,
    v_test_statistic_b,
    v_hypothesis_test_dm,
    v_test_statistic_dm
  )

  plotting_data <- tibble(
    M = M_set,
    b_size_distortion = mat_size_distortion_b[,1],
    b_power_loss = mat_power_loss_b[,1]
  ) %>%
    left_join(., df_hypoth_testing, by = c("M" = "v_M")) %>%
    mutate(
      v_hypothesis_test_b = case_when(
        v_hypothesis_test_b == TRUE ~ "cross",
        v_hypothesis_test_b == FALSE ~ "circle"
      )
    )

  plot <- ggplot(plotting_data, aes(x = b_size_distortion, y = b_power_loss)) +
    geom_path(size = 1, linetype = "dashed") +
    geom_point(aes(shape = v_hypothesis_test_b), size = 4.5, color = "red", stroke = 1.5) +
    scale_shape_identity() +
    geom_text(aes(label = M), nudge_y = .005) +
    labs(
      x = "Size Distortion",
      y = "Maximum Power Loss",
      title = series,
      caption = str_glue("Data from {starting_year} to {ending_year}\nHorizon {horizon}")
    ) +
    theme_minimal()

  ggsave(plot = plot, filename = str_glue("spf_tradeoff_plots/{series}_{horizon}_{starting_year}_{ending_year}.png"),
         width = 7, height = 7, units = "in")
}

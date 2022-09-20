
pacman::p_load(
  tidyverse,
  lubridate,
  kableExtra
)

Load_Data <- function(dataset) {

   read_csv(paste0(dataset, "_extended.csv")) %>%
    select(X1, starts_with("SPF"), starts_with("NC"), Realiz1) %>%
    mutate(
      year_quarter  = str_replace(X1, ":", "-") %>% yq(.),
      across(-year_quarter, ~ as.double(.))
    ) %>%
    select(-X1) %>%
    select(year_quarter, everything())

}

gdp <- Load_Data("RDGP") %>% glimpse()
inflation <- Load_Data("PGDP") %>% glimpse()
unemp <- Load_Data("UNEMP") %>% glimpse()
tbill <- Load_Data("TBILL") %>% glimpse()

Replicate_Table <- function(df, starting_year_quarter, ending_year_quarter, horizon, bandwidth = 3) {

  start <- as.Date(starting_year_quarter)
  end <- as.Date(ending_year_quarter)

  eval_period <- seq.Date(start, end, by = "quarter")

  filtered_data <- df %>%
    filter(year_quarter %in% eval_period)

  NC_name <- paste0("NCfor_Step", horizon)
  SPF_name <- paste0("SPFfor_Step", horizon)

  e1 = filtered_data$Realiz1 - (filtered_data %>% pull(NC_name))  #forecast error from forecaster 1
  e1[is.na(e1)] <- 0

  e2 = filtered_data$Realiz1 - (filtered_data %>% pull(SPF_name)) #forecast error from forecaster 2
  e2[is.na(e2)] <- 0

  d  = e1^2 - e2^2  #squared loss differential
  n = length(d);

  wce_dm = dm.test.r(d, cl = 0.05, h=horizon) #dm original (WCE-DM)
  wce_b = dm.test.bt.fb(d, cl = 0.05, M = floor(sqrt(n))); #fixed-b (WCE-B)
  wpe_d = dm.test.wpe.fb(d, cl=0.05, M = floor(n^(1/bandwidth))); # (WPE-D)

  wce_dm_ten = dm.test.r(d, cl = .10, h=horizon) #dm original (WCE-DM)
  wce_b_ten = dm.test.bt.fb(d, cl = .10, M = floor(sqrt(n))); #fixed-b (WCE-B)
  wpe_d_ten = dm.test.wpe.fb(d, cl=.10, M = floor(n^(1/bandwidth))); # (WPE-D)

  df_stats <- tibble("Forecast Horizon" = horizon - 1,
         "WCD-DM" = wce_dm$stat,
         "WCE-B" = wce_b$stat,
         "WPE-D" = wpe_d$stat) %>%
    mutate(
      across(
        everything(), ~ round(., digits = 2)
      ),
      "Forecast Horizon" = as.character(`Forecast Horizon`),

    )

  df_reject <- tibble("Forecast Horizon" = horizon - 1,
         "WCD-DM_05" = wce_dm$rej,
         "WCD-DM_10" = wce_dm_ten$rej,
         "WCE-B_05" = wce_b$rej,
         "WCE-B_10" = wce_b_ten$rej,
         "WPE-D_05" = wpe_d$rej,
         "WPE-D_10" = wpe_d_ten$rej
  ) %>%
    mutate(
      "Forecast Horizon" = as.character(`Forecast Horizon`),
    )

  full_join(df_stats, df_reject, by = "Forecast Horizon")

}

Add_Asterisks <- function(df) {

  df %>%
    mutate(
      "WCE-B" = case_when(`WCE-B_05` == TRUE ~ paste0(`WCE-B`, "**"),
                          `WCE-B_10` == TRUE ~ paste0(`WCE-B`, "*"),
                          TRUE ~ as.character(`WCE-B`)),
      "WPE-D" = case_when(`WPE-D_05` == TRUE ~ paste0(`WPE-D`, "**"),
                          `WPE-D_10` == TRUE ~ paste0(`WPE-D`, "*"),
                          TRUE ~ as.character(`WPE-D`)),
      "WCD-DM" = case_when(`WCD-DM_05` == TRUE ~ paste0(`WCD-DM`, "**"),
                          `WCD-DM_10` == TRUE ~ paste0(`WCD-DM`, "*"),
                          TRUE ~ as.character(`WCD-DM`)),
    ) %>%
  select("Forecast Horizon",
         "WCD-DM",
         "WCE-B, M = $T^{1/2}$" = "WCE-B",
         "WPE-D, m = $T^{1/3}$" = "WPE-D"
         ) %>%
  t()

}

Save_Kable <- function(df, file_name) {

  caption_ <- switch (file_name,
    "gdp" = "Real output growth: SPF versus random walk",
    "inf" = "GNP/GDP inflation: SPF versus random walk",
    "unemp" = "Unemployment rate: SPF versus random walk",
    "tbill" = "Three-month Treasury bill: SPF versus random walk"
  )

  kable(df, format = "latex", booktabs = TRUE, caption = caption_, escape = FALSE) %>%
    # row_spec(c(1, 4, 7, 10, 13)) %>%
    # kable_material("striped") %>%
    pack_rows("Evaluation Period: 1987Q1 - 2016Q4, T = 120", 2, 4) %>%
    pack_rows("Evaluation Period: 1987Q1 - 2021Q4, T = 140", 5, 7) %>%
    pack_rows("Evaluation Period: 1987Q1 - 1996, T = 40", 8, 10) %>%
    pack_rows("Evaluation Period: 1997Q1 - 2006Q4, T = 40", 11, 13) %>%
    pack_rows("Evaluation Period: 2007Q1 - 2016Q4, T = 40", 14, 16) %>%
    pack_rows("Evaluation Period: 2017Q1 - 2021Q1, T = 20", 17, 19) %>%
    footnote(general = "Asterisks ** and * indicate two-sided significance at respectively, the 5% and 10% level using fixed-basymptotics for WCE-B and fixed-m asymptotics for WPE-D.",
             threeparttable = TRUE) %>%
    kable_styling(font_size = 4) %>%
    save_kable(., paste0("../../note/note02_slides/tables/", file_name, ".tex"))

}

v_rows_to_keep <- c(1:4, 6:8, 10:12, 14:16, 18:20, 22:24)

# GDP ---------------------------------------------------------------------

gdp_1 <- map_dfr(1:5, ~ Replicate_Table(gdp, "1987-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
gdp_2 <- map_dfr(1:5, ~ Replicate_Table(gdp, "1987-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)
gdp_3 <- map_dfr(1:5, ~ Replicate_Table(gdp, "1987-01-01", "1996-10-01", .)) %>% Add_Asterisks(.)
gdp_4 <- map_dfr(1:5, ~ Replicate_Table(gdp, "1997-01-01", "2006-10-01", .)) %>% Add_Asterisks(.)
gdp_5 <- map_dfr(1:5, ~ Replicate_Table(gdp, "2007-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
gdp_6 <- map_dfr(1:5, ~ Replicate_Table(gdp, "2017-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)

mat_gdp <- rbind(gdp_1, gdp_2, gdp_3, gdp_4, gdp_5, gdp_6)[v_rows_to_keep, ]

Save_Kable(mat_gdp, "gdp")

# Inflation ---------------------------------------------------------------------

inf_1 <- map_dfr(1:5, ~ Replicate_Table(inflation, "1987-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
inf_2 <- map_dfr(1:5, ~ Replicate_Table(inflation, "1987-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)
inf_3 <- map_dfr(1:5, ~ Replicate_Table(inflation, "1987-01-01", "1996-10-01", .)) %>% Add_Asterisks(.)
inf_4 <- map_dfr(1:5, ~ Replicate_Table(inflation, "1997-01-01", "2006-10-01", .)) %>% Add_Asterisks(.)
inf_5 <- map_dfr(1:5, ~ Replicate_Table(inflation, "2007-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
inf_6 <- map_dfr(1:5, ~ Replicate_Table(inflation, "2017-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)

mat_inf <- rbind(inf_1, inf_2, inf_3, inf_4, inf_5, inf_6)[v_rows_to_keep, ]

Save_Kable(mat_inf, "inf")


# Unemployment ------------------------------------------------------------

unemp_1 <- map_dfr(1:5, ~ Replicate_Table(unemp, "1987-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
unemp_2 <- map_dfr(1:5, ~ Replicate_Table(unemp, "1987-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)
unemp_3 <- map_dfr(1:5, ~ Replicate_Table(unemp, "1987-01-01", "1996-10-01", .)) %>% Add_Asterisks(.)
unemp_4 <- map_dfr(1:5, ~ Replicate_Table(unemp, "1997-01-01", "2006-10-01", .)) %>% Add_Asterisks(.)
unemp_5 <- map_dfr(1:5, ~ Replicate_Table(unemp, "2007-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
unemp_6 <- map_dfr(1:5, ~ Replicate_Table(unemp, "2017-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)

mat_unemp <- rbind(unemp_1, unemp_2, unemp_3, unemp_4, unemp_5, unemp_6)[v_rows_to_keep, ]

Save_Kable(mat_unemp, "unemp")

# TBill -------------------------------------------------------------------

tbill_1 <- map_dfr(1:5, ~ Replicate_Table(tbill, "1987-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
tbill_2 <- map_dfr(1:5, ~ Replicate_Table(tbill, "1987-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)
tbill_3 <- map_dfr(1:5, ~ Replicate_Table(tbill, "1987-01-01", "1996-10-01", .)) %>% Add_Asterisks(.)
tbill_4 <- map_dfr(1:5, ~ Replicate_Table(tbill, "1997-01-01", "2006-10-01", .)) %>% Add_Asterisks(.)
tbill_5 <- map_dfr(1:5, ~ Replicate_Table(tbill, "2007-01-01", "2016-10-01", .)) %>% Add_Asterisks(.)
tbill_6 <- map_dfr(1:5, ~ Replicate_Table(tbill, "2017-01-01", "2021-10-01", .)) %>% Add_Asterisks(.)

mat_tbill <- rbind(tbill_1, tbill_2, tbill_3, tbill_4, tbill_5, tbill_6)[v_rows_to_keep, ]

Save_Kable(mat_tbill, "tbill")

# Robustness Checks -------------------------------------------------------


if (type == "Robustness Check") {

wce_b = dm.test.bt.fb(d, cl = 0.05, M = bandwidth); #fixed-b (WCE-B)
wpe_d = dm.test.wpe.fb(d, cl=0.05, M = bandwidth); # (WPE-D)
wce_dm = dm.test.bt(d, cl = .05, M = bandwidth)

wce_b_ten = dm.test.bt.fb(d, cl = .10, M = bandwidth); #fixed-b (WCE-B)
wpe_d_ten = dm.test.wpe.fb(d, cl=.10, M = bandwidth); # (WPE-D)
wce_dm_ten = dm.test.bt(d, cl = .05, M = bandwidth)


}

num_obs <- nrow(gdp)
b <- seq(.1, 1, .1)
m <- floor(b * num_obs)

v_bandwidth <- c(1:5, m)
horizon <- 1:5

df_horizon_bandwidths <- expand_grid(horizon, v_bandwidth)

Run_Robustness_Check <- function(df_string, function_name) {

  df <- get(df_string)

  caption_ <- switch (df_string,
    "gdp" = str_glue("{function_name}, Real Output Growth"),
    "inflation" = str_glue("{function_name}, GDP/GNP Inflation"),
    "unemp" = str_glue("{function_name}, Unemploymnet Rate"),
    "tbill" = str_glue("{function_name}, Three Month Treasury Bill")
  )

  file_name <- paste0(df_string, "_", function_name, "_robust.tex")

  df_table <- map2_dfr(df_horizon_bandwidths$horizon, df_horizon_bandwidths$v_bandwidth,
         ~ Replicate_Table(df,
                           "1987-01-01",
                           "2016-10-01",
                           horizon = .x,
                           bandwidth = .y,
                           type = "Robustness Check") %>%
  mutate(
    m = .y
    )
  ) %>%
  as_tibble() %>%
  select(`Forecast Horizon`, m, starts_with("WPE"), starts_with("WCE-B"), starts_with("DM")) %>%
  mutate(
    "WPE-D" = case_when(`WPE-D_05` == TRUE ~ paste0(`WPE-D`, "**"),
                        `WPE-D_10` == TRUE ~ paste0(`WPE-D`, "*"),
                        TRUE ~ as.character(`WPE-D`)),
    "WCE-B" = case_when(`WCE-B_05` == TRUE ~ paste0(`WCE-B`, "**"),
                        `WCE-B_10` == TRUE ~ paste0(`WCE-B`, "*"),
                        TRUE ~ as.character(`WCE-B`)),
    "DM-BT" = case_when(`WCD-DM_05` == TRUE ~ paste0(`DM`, "**"),
                          `WCD-DM_10` == TRUE ~ paste0(`DM`, "*"),
                          TRUE ~ as.character(`DM`))
  ) %>%
  select("Forecast Horizon", "m", all_of(function_name))

  list(robustness_table = df_table,
       table_caption = caption_,
       output_file_name = file_name)

}

Save_Robustness_Table <- function(list_table_caption_file) {

  df <- list_table_caption_file %>%
    pluck("robustness_table")

  caption_ <- list_table_caption_file %>%
    pluck("table_caption")

  file_name <- list_table_caption_file %>%
    pluck("output_file_name")

  df_to_save <- pivot_wider(df, names_from = `Forecast Horizon`, values_from = 3) %>%
    mutate(
      m = as.character(m),
      across(everything(.), ~ as.character(.))
    ) %>%
    rename("Forecast Horizon" = m) %>%
    add_row(
      "Forecast Horizon" = "m = ",
      `0` = "",
      `1` = "",
      `2` = "",
      `3` = "",
      `4` = "",
      .before = 1)

  kable(df_to_save, format = "latex", escape = FALSE, caption = caption_, booktabs = TRUE) %>%
    column_spec(., 1, border_right = TRUE) %>%
    kable_styling(font_size = 6) %>%
    footnote(
      general = "Asterisks ** and * indicate two-sided significance at, respectively, the 5% and 10% level using fixed-m asymptotics. m is chosen as $m \\in{1,2,3,4,5,floor(b * t)}$ with $t = 219$ and $b \\in{.1, .2, \\dots,.9, 1.0}$",
             threeparttable = TRUE) %>%
    save_kable(., paste0("../../note/note02_slides/tables/", file_name))
}

df_robust_combinations <- expand_grid(dataset = c("gdp", "inflation", "unemp", "tbill"),
                                      dm = c("WPE-D", "WCE-B", "DM-BT"))

Run_Robustness_Check("gdp", "WPE-D")

walk2(df_robust_combinations$dataset, df_robust_combinations$dm, ~ Run_Robustness_Check(.x, .y) %>% Save_Robustness_Table(.))




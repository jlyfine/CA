#------load packages------
rm(list = ls());gc()
require(extRemes);require(VineCopula);require(abind);require(segmented);
require(lmom);require(foreach);require(doParallel);require(parallel);
require(snow);require(ncdf4);require(countrycode)
#-----------------------------
# file path input
#-----------------------------

base_dir = "/share/home/junli2025/24Attribution/0data/0base_data"
#base_dir = "D:/23Attribution_onlyco2/0data"

co2_loo_file = file.path(
  base_dir,
  "CO2_ppm_from_each_country_due to fossil fuel from 1850-2023")

obs_co2_file = file.path(
  base_dir,
  "CO2_ppm_historical_period_1850_2023")

grid_file = file.path(
  base_dir,
  "global_grid_2.5degree.nc")

#-----------------------------
# Read country-level CO2 data
#-----------------------------

co2_loo = readRDS(co2_loo_file)

country_names = colnames(co2_loo)[1:214]

country_names_iso3 = gsub("0", " ", country_names)
country_names_iso3 = countrycode(country_names_iso3, "country.name", "iso3c")  
colnames(co2_loo)[1:214] = country_names_iso3

# If needed:
# country_iso3[is.na(country_iso3)] <- "XKX"  # Kosovo

#-----------------------------
# Read observed historical CO2
#-----------------------------

obs_co2_raw = readRDS(obs_co2_file)
obs_co2 = obs_co2_raw$co2 - obs_co2_raw$co2[1]

co2_without_group = co2_loo
co2_without_group[] = NA
for (ik in 1:ncol(co2_loo)) {#ik=215
  co2_without_group[,ik] = obs_co2 - co2_loo[,ik]
}

#-----------------------------
# Present-day CO2 means
#-----------------------------

present_window = 20

present_rows = (nrow(co2_without_group) - present_window + 1):
  nrow(co2_without_group)

co2_model_all_group_loo = mean(
  co2_without_group[present_rows, ncol(co2_without_group)],
  na.rm = TRUE)

co2_model_each_group_loo = colMeans(
  co2_without_group[present_rows, 1:(ncol(co2_without_group) - 1)],
  na.rm = TRUE)

#-----------------------------
# Observed present-day CO2
#-----------------------------
# 2004–2023 corresponds to rows 155–174 if the series starts in 1850.

present_start_year = 2004
present_end_year   = 2023

present_year_index = (present_start_year - 1850 + 1):
  (present_end_year - 1850 + 1)

co2_obs_present_day = mean(
  obs_co2[present_year_index], 
  na.rm = TRUE
)

co2_model_present_day = rep(
  co2_obs_present_day,
  length(co2_model_each_group_loo)
)

#-----------------------------
# Read grid information
#-----------------------------

grid_nc = nc_open(grid_file)
nc_close(grid_nc)

#-----------------------------
# CO2 during 1950–2023
#-----------------------------

analysis_start_year = 1950
analysis_end_year   = 2023

analysis_year_index = (
  analysis_start_year - 1850 + 1
):(
  analysis_end_year - 1850 + 1
)

co2_obs_1950_2023 = obs_co2[analysis_year_index]

message("co2_obs_1950_2023 in")

#===========================================================
# 1.1 Check fitted GEV model
#============================================================

check_gev_fit = function(mydata, par) {
  
  y = mydata$y
  X.loc = model.matrix(~x, data = mydata)
  X.sc  = model.matrix(~x, data = mydata)
  
  p.loc = ncol(X.loc)
  p.sc  = ncol(X.sc)
  
  beta.mu = par[1:p.loc]
  beta.sc = par[(p.loc + 1):(p.loc + p.sc)]
  xi      = par[length(par)]
  
  mu = as.vector(X.loc %*% beta.mu)
  sigma = as.vector(X.sc %*% beta.sc)
  
  z = (y - mu) / sigma
  
  if (abs(xi) < 1e-6) {
    p = exp(-exp(-z))
  } else {
    t = 1 + xi * z
    p = rep(NA, length(y))
    ok = t > 0
    p[ok] = exp(-t[ok]^(-1 / xi))
    p[!ok & xi > 0] <- 0
    p[!ok & xi < 0] <- 1
  }
  
  return(list(
    mu = mu,
    sigma = sigma,
    xi = xi,
    u = p
  ))
}

#============================================================
# 1.2 Fit GEV models for one grid cell
#============================================================

fit_gev_one_cell = function(y, x) {
  #tada = readRDS("E:/23Attribution/0data/0obs_standardized_index/sti_1month_scale_era5_1950-2023_2.5degree_summer")
  #y = tada[120,40,]; x = as.vector(co2_summer)
  mydata = data.frame(y = y, x = x)
  mydata = na.omit(mydata)
  
  fitx = fevd(y, data = mydata, location.fun = ~ x, scale.fun = ~ x)
  
  par_summary = summary.fevd(fitx, silent = TRUE)
  
  fit_check = check_gev_fit(mydata, par_summary$par)

  u = sort(fit_check$u)
  
  q_theory = qnorm((1:length(u)) / (length(u) + 1))
  q_emp    = qnorm(u)
  
  fitness_p = ks.test(q_theory, q_emp)$p.value
  
  list(
    fitness_p = fitness_p,
    par       = par_summary$par
  )
}

#============================================================
# 1.2 Fit GEV models for all grid cells
#============================================================

#tada = readRDS("E:/23Attribution/0data/0obs_standardized_index/sti_1month_scale_era5_1950-2023_2.5degree_summer")
#y = tada[120,40,]; x = as.vector(co2_summer)
cal_gev_i = function(tada, co2_obs_1950_2023, whether_sample = "no") {
  
  co2_monthly = matrix(
    rep(co2_obs_1950_2023, each = 12),
    ncol = length(co2_obs_1950_2023)
  )
  co2_summer = as.vector(co2_monthly[c(6, 7, 8), ])
  
  n_time = length(co2_summer)
  n_lat  = dim(tada)[1]
  n_lon  = dim(tada)[2]
  
  if (whether_sample == "yes") {
    sample_index = sample(seq_len(n_time), replace = TRUE)
  } else {
    sample_index = seq_len(n_time)
  }
  
  outda = array(NA, dim = c(n_lat, n_lon, n_time))
  template = matrix(NA, nrow = n_lat, ncol = n_lon)
  
  fitness_gev = template
  par_mu0     = template
  par_mu1     = template
  par_sigma0  = template
  par_sigma1  = template
  par_shape   = template
  
  min_valid_n = n_time / 12
  
  for (i in seq_len(n_lat)) {
    for (j in seq_len(n_lon)) {
      if (is.na(tada[i, j, 1])) next#i=120;j=40
      
      y_full = tada[i, j, ]
      x_full = co2_summer
      
      mydata_full = data.frame(y = y_full, x = x_full)
      
      if (nrow(na.omit(mydata_full)) < min_valid_n) next
      
      y_sample = y_full[sample_index]
      x_sample = x_full[sample_index]
      
      outda[i, j, ] = y_sample
      
      fit_result = try(
        fit_gev_one_cell(y = y_sample, x = x_sample),
        silent = TRUE
      )
      
      if (inherits(fit_result, "try-error")) next
      
      fitness_gev[i, j] = fit_result$fitness_p
      
      par_mu0   [i, j] = fit_result$par[1]
      par_mu1   [i, j] = fit_result$par[2]
      par_sigma0[i, j] = fit_result$par[3]
      par_sigma1[i, j] = fit_result$par[4]
      par_shape [i, j] = fit_result$par[5]
    }
  }
  
  list(
    outda = outda,
    
    fitness_gev = fitness_gev,
    
    par_mu0    = par_mu0,
    par_mu1    = par_mu1,
    par_sigma0 = par_sigma0,
    par_sigma1 = par_sigma1,
    par_shape  = par_shape,
    
    sample_index0 = sample_index
  )
}

#============================================================
# 1.4 Run repeated GEV fitting
#============================================================

cal_gev_fitness = function(
    fpdain = "E:/23Attribution/0data/0obs_standardized_index/sti_1month_scale_era5_1950-2023_2.5degree_summer",
    fpout  = "E:/23Attribution/0results/",
    rep_times = 500,
    co2_obs_1950_2023 = co2_obs_1950_2023,
    core = 6
) {
  
  tada = readRDS(fpdain)
  
  input_name = basename(fpdain)
  
  sample_types = c("no", rep("yes", rep_times))
  
  start_time = proc.time()
  
  cl = makeCluster(core)
  registerDoParallel(cl)
  
  #"check_gev_fit"; "make_summer_co2_series"
  #"fit_gev_one_cell"; "cal_gev_i"
  foreach(
    k = seq_along(sample_types),
    .packages = c("extRemes", "ncdf4"),
    .export = c(
      "check_gev_fit",
      "fit_gev_one_cell",
      "cal_gev_i"
    )
  ) %dopar% {#k=1
    
    sample_type = sample_types[k]
    
    gev_result = cal_gev_i(
      tada = tada,
      co2_obs_1950_2023 = co2_obs_1950_2023,
      whether_sample = sample_type
    )

    output_file = file.path(
      fpout,
      paste0(
        input_name,
        "_gev_fitness_times",
        k,
        "_sample_",
        sample_type
      )
    )
    #gev_result = ''
    saveRDS(gev_result, file = output_file, compress = "xz")
  }
  
  stopCluster(cl)
  
  end_time = proc.time()
  
  message(
    "Time used: ",
    round((end_time - start_time)[3] / 60, 2),
    " minutes"
  )
}

#============================================================
# 2.1 Calculate compound-event probability ratios
#============================================================

calpr_compound = function(par_ta, par_sm, threshold, copula_future, tas_i) {
  
  get_gev_par = function(par, tas_value) {
    c(
      location = par[1] + par[2] * tas_value,
      scale    = par[3] + par[4] * tas_value,
      shape    = par[5]
    )
  }
 #parx_his = c(parx[1] + parx[2]*0    , parx[3] + parx[4]*0    , parx[5])
 #pary_his = c(pary[1] + pary[2]*0    , pary[3] + pary[4]*0    , pary[5])
 #parx_fut = c(parx[1] + parx[2]*tas_i, parx[3] + parx[4]*tas_i, parx[5])
 #pary_fut = c(pary[1] + pary[2]*tas_i, pary[3] + pary[4]*tas_i, pary[5])
  
  safe_prob = function(p, eps = 1e-9) {
    ifelse(p == 0, eps, p)
  }

  par_ta_his = get_gev_par(par_ta, tas_value = 0)
  par_sm_his = get_gev_par(par_sm, tas_value = 0)
  
  par_ta_fut = get_gev_par(par_ta, tas_value = tas_i)
  par_sm_fut = get_gev_par(par_sm, tas_value = tas_i)
  
  invalid_scale = any(c(
    par_ta_his["scale"],
    par_sm_his["scale"],
    par_ta_fut["scale"],
    par_sm_fut["scale"]
  ) <= 0)
  
  if (invalid_scale) {
    return(c(rep(NA, 5), tas_i))
  }

  ta_mean_scale = pevd(
    qevd(threshold, par_ta_his["location"], par_ta_his["scale"], par_ta_his["shape"]),
    par_ta_fut["location"], par_ta_fut["scale"], par_ta_fut["shape"]
    )
  
  sm_mean_scale = pevd(
    qevd(threshold, par_sm_his["location"], par_sm_his["scale"], par_sm_his["shape"]),
    par_sm_fut["location"], par_sm_fut["scale"], par_sm_fut["shape"]
    )
  
  ta_scale_only = pevd(
    qevd(threshold, par_ta_his["location"], par_ta_his["scale"], par_ta_his["shape"]),
    par_ta_his["location"], par_ta_fut["scale"], par_ta_his["shape"]
  )
  
  sm_scale_only = pevd(
    qevd(threshold, par_sm_his["location"], par_sm_his["scale"], par_sm_his["shape"]),
    par_sm_his["location"], par_sm_fut["scale"], par_sm_his["shape"]
  )
  
  p_joint_future = 1 - sm_mean_scale - ta_mean_scale +
    BiCopCDF(sm_mean_scale, ta_mean_scale, copula_future)
  
  p_joint_baseline = 1 - threshold - threshold +
    BiCopCDF(threshold, threshold, copula_future)
  
  p_sm_mean_scale = 1 - sm_mean_scale - threshold +
    BiCopCDF(sm_mean_scale, threshold, copula_future)
  
  p_ta_mean_scale = 1 - threshold - ta_mean_scale +
    BiCopCDF(threshold, ta_mean_scale, copula_future)
  
  p_sm_scale_only = 1 - sm_scale_only - threshold +
    BiCopCDF(sm_scale_only, threshold, copula_future)
  
  p_ta_scale_only = 1 - threshold - ta_scale_only +
    BiCopCDF(threshold, ta_scale_only, copula_future)
  
  p_joint_future   = safe_prob(p_joint_future)
  p_joint_baseline = safe_prob(p_joint_baseline)
  p_sm_mean_scale  = safe_prob(p_sm_mean_scale)
  p_ta_mean_scale  = safe_prob(p_ta_mean_scale)
  p_sm_scale_only  = safe_prob(p_sm_scale_only)
  p_ta_scale_only  = safe_prob(p_ta_scale_only)
  
  c(
    joint_mean_scale_ratio = p_joint_future / p_joint_baseline,
    sm_mean_ratio          = p_sm_mean_scale / p_sm_scale_only,
    ta_mean_ratio          = p_ta_mean_scale / p_ta_scale_only,
    sm_scale_ratio         = p_sm_scale_only / p_joint_baseline,
    ta_scale_ratio         = p_ta_scale_only / p_joint_baseline,
    tas                    = tas_i
  )
}

#============================================================
# 2.2 Extract GEV parameters at one grid cell
#============================================================

get_cell_gev_par = function(gev_fit, i, j) {
  
  c(
    mu0    = gev_fit$par_mu0[i, j],
    mu1    = gev_fit$par_mu1[i, j],
    sigma0 = gev_fit$par_sigma0[i, j],
    sigma1 = gev_fit$par_sigma1[i, j],
    shape  = gev_fit$par_shape[i, j]
  )
}

#============================================================
# 2.3 Fit compound LMF for one pair of GEV-fit files
#============================================================

cal_fit_compound <- function(
    ta_gev_fit,
    sm_gev_fit,
    threshold,
    co2_obs_present_day,
    output_length,
    tas_sequence,
    family_use
) {
  #sm_gev_fit = readRDS("D:/23Attribution_onlyco2/0results/1hot_dry_gev_fitness_only_summer/spi/spi_opposite_1month_scale_era5_1950-2023_2.5degree_summer_gev_fitness_times1_sample_no")
  #ta_gev_fit = readRDS("D:/23Attribution_onlyco2/0results/1hot_dry_gev_fitness_only_summer/sti/sti_1month_scale_era5_1950-2023_2.5degree_summer_gev_fitness_times1_sample_no")
  #output_length = 150
  #family_use = 5
  #tas_sequence = seq(0, 110, length.out=150)
  ta_data = ta_gev_fit$outda
  sm_data = sm_gev_fit$outda
  
  n_lat = dim(ta_data)[1]
  n_lon = dim(ta_data)[2]
  
  lmf_all = array(NA, dim = c(n_lat, n_lon, output_length))
  
  for (i in seq_len(n_lat)) {
    for (j in seq_len(n_lon)) {
      
      ta_cell = ta_data[i, j, ]
      sm_cell = sm_data[i, j, ]
      
      if (all(is.na(ta_cell)) || all(is.na(sm_cell))) next
      
      par_ta = c(ta_gev_fit$par_mu0   [i,j],
                 ta_gev_fit$par_mu1   [i,j],
                 ta_gev_fit$par_sigma0[i,j], 
                 ta_gev_fit$par_sigma1[i,j], 
                 ta_gev_fit$par_shape [i,j])  
      par_sm = c(sm_gev_fit$par_mu0   [i,j],
                 sm_gev_fit$par_mu1   [i,j],
                 sm_gev_fit$par_sigma0[i,j], 
                 sm_gev_fit$par_sigma1[i,j], 
                 sm_gev_fit$par_shape [i,j]) 

      if (anyNA(par_ta) || anyNA(par_sm)) next
      
      copula_par = BiCopEst(
        pobs(ta_cell),
        pobs(sm_cell),
        family = family_use
      )
      
      copula_future = BiCop(
        family = family_use,
        par    = copula_par$par,
        par2   = copula_par$par2
      )
      
      lmf_by_tas = c()
      for (k in 1:length(tas_sequence)) {#k=50
        lmf_by_tas = c(lmf_by_tas, 
                       calpr_compound(par_ta        = par_ta,
                                      par_sm        = par_sm,
                                      threshold     = threshold,
                                      copula_future = copula_future,
                                      tas_i         = tas_sequence[k]))
      }
      ncolda = length(lmf_by_tas)/length(tas_sequence)
      lmf_by_tas = matrix(lmf_by_tas, ncol=ncolda, byrow=T)
      
      lmf_present_day = calpr_compound(
        par_ta        = par_ta,
        par_sm        = par_sm,
        threshold     = threshold,
        copula_future = copula_future,
        tas_i         = co2_obs_present_day
      )
      
      lmf_all[i, j, ] = c(
        matrix(lmf_by_tas, ncol=1),
        lmf_present_day[1:5],
        copula_par$par,
        copula_par$par2
      )
    }
  }
  
  lmf_all
}

#============================================================
# 2.4 Write NetCDF output
#============================================================

write_lmf_netcdf = function(
    lon,
    lat,
    var_name,
    var_unit,
    long_name,
    time_unit,
    output_file,
    data,
    output_dir
) {
  
  lon_dim = ncdim_def(
    name     = "lon",
    units    = "degrees_east",
    vals     = lon,
    longname = "longitude"
  )
  
  lat_dim = ncdim_def(
    name     = "lat",
    units    = "degrees_north",
    vals     = lat,
    longname = "latitude"
  )
  
  time_dim = ncdim_def(
    name     = "time",
    units    = time_unit,
    calendar = "standard",
    vals     = seq_len(dim(data)[3]) - 1
  )
  
  var_def = ncvar_def(
    name        = var_name,
    units       = var_unit,
    dim         = list(lon_dim, lat_dim, time_dim),
    longname    = long_name,
    compression = 9
  )
  
  output_path = file.path(output_dir, output_file)
  
  nc = nc_create(output_path, var_def)
 #on.exit(nc_close(nc), add = TRUE)
  
  ncvar_put(nc, var_def, data)
  ncatt_put(nc, varid = 0, attname = "description", attval = long_name)
  nc_close(nc)
  
}

#============================================================
# 2.5 Main LMF calculation
#============================================================

cal_LMF = function(
    fpin_sti = "D:/23Attribution_onlyco2/0results/1hot_dry_gev_fitness_only_summer/sti",
    fpin_spi = "D:/23Attribution_onlyco2/0results/1hot_dry_gev_fitness_only_summer/spi",
    fpout = "E:/23Attribution/0results/",
    modoutna = "fit_hot_wet_era5",
    family_use = 5,
    dasta = 1,
    daend = 5,
    threshold = 0.9,
    core = 5,
    co2_dis = 150,
    grid_nc,
    co2_obs_present_day
) {
  
  tas_sequence = seq(
    from = 0,
    to   = co2_obs_present_day,
    length.out = co2_dis
  )
  
  n_output_per_tas = 6
  
  output_length = length(tas_sequence) * n_output_per_tas + 7
  
  sti_files = file.path(fpin_sti, list.files(fpin_sti))[dasta:daend]
  spi_files = file.path(fpin_spi, list.files(fpin_spi))[dasta:daend]
  
  ta_gev_fit_all = lapply(sti_files, readRDS)
  sm_gev_fit_all = lapply(spi_files, readRDS)

  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_output_all = foreach(
    k = seq_along(ta_gev_fit_all),
    .packages = c("VineCopula", "extRemes", "foreach", "ncdf4", "lmom"),
    .export = c(
      "calpr_compound",
      "get_cell_gev_par",
      "cal_fit_compound",
      "write_lmf_netcdf"
    )
  ) %dopar% {#k=1
    
    lmf_output = cal_fit_compound(
      ta_gev_fit          = ta_gev_fit_all[[k]],
      sm_gev_fit          = sm_gev_fit_all[[k]],
      threshold           = threshold,
      co2_obs_present_day = co2_obs_present_day,
      output_length       = output_length,
      tas_sequence        = tas_sequence,
      family_use          = family_use)
    
    lmf_output[is.infinite(lmf_output)] = 9999

    lmf_output
  }
  
  stopCluster(cl)
  
  for (k in seq_along(lmf_output_all)) {#k=1
    
    file_parts = strsplit(basename(sti_files[k]), "_")[[1]]
    
    start_id = which(grepl("times", file_parts))
    end_id   = which(grepl("yes|no", file_parts))
    
    sample_name = paste(file_parts[start_id:end_id], collapse = "_")
    
    output_file = paste0(
      modoutna,
      "_1950_2023_thre",
      threshold,
      "_",
      sample_name,
      ".nc"
    )
    
    write_lmf_netcdf(
      lon         = grid_nc$dim$lon$vals,
      lat         = rev(grid_nc$dim$lat$vals),
      var_name    = "NA",
      var_unit    = "NA",
      long_name   = "NA",
      time_unit   = "months since 1850-01",
      output_file = output_file,
      data        = lmf_output_all[[k]],
      output_dir  = fpout
    )
  
  }
}

#============================================================
# 3.1 Fit segmented model and make leave-one-out predictions
#============================================================

fit_segmented_prediction = function(
    response_data,
    co2_model_each_group_loo,
    co2_model_present_day,   
    co2_model_all_group_loo, 
    n_breakpoints) {
  
  seg_model = tryCatch(
    segmented(
      lm(y ~ x, data = response_data),
      seg.Z = ~ x,
      npsi = n_breakpoints
    ),
    error = function(e) NULL
  )
  
  if (is.null(seg_model)) return(NULL)
  
  slope_values = tryCatch(slope(seg_model)$x[, 1], error = function(e) NULL)
  intercept_values = tryCatch(intercept(seg_model)$x[, 1], error = function(e) NULL)
  model_summary = tryCatch(summary(seg_model), error = function(e) NULL)
  
  if (
    is.null(slope_values) ||
    is.null(intercept_values) ||
    is.null(model_summary)
  ) {
    return(NULL)
  }
  
  co2_model_each_group_loo = as.numeric(co2_model_each_group_loo)
  
  lmf_group_loo = predict(seg_model, 
                          newdata = data.frame(x = co2_model_present_day)) - 
    predict(seg_model,
            newdata = data.frame(x = co2_model_each_group_loo))
  
  lmf_all_group_loo = predict(seg_model,
                              newdata = data.frame(x = unique(co2_model_present_day))) - 
    predict(seg_model, 
            newdata = data.frame(x = co2_model_all_group_loo))
  
  obs_mod_last_value = c(
    observed_last = tail(response_data$y, 1),
    fitted_last   = tail(seg_model$fitted.values, 1)
  )
  
  #lmf_group_loo[lmf_group_loo < 0] = 0
  #lmf_all_group_loo[lmf_all_group_loo < 0] = 0
  
  list(
    lmf_group_loo    = lmf_group_loo,
    lmf_all_group_loo = lmf_all_group_loo,
    slope            = slope_values,
    intercept        = intercept_values,
    r_squared        = model_summary$r.squared,
    obs_mod_last_value = obs_mod_last_value
  )
}

#============================================================
# 3.2 Apply segmented prediction to all grid cells
#============================================================

run_segmented_prediction_grid = function(
    co2_model_each_group_loo,
    co2_model_present_day,   
    co2_model_all_group_loo, 
    lmf_fit,
    n_breakpoints
) {
  #lmf_fit = lmf_fit_all[[k]]
  lmf_total  = lmf_fit$lmf_tot
  tas_change = lmf_fit$tas_change
  
  n_lon = dim(lmf_total)[1]
  n_lat = dim(lmf_total)[2]
  n_group = length(co2_model_each_group_loo)
  
  lmf_group_loo = array(NA, dim = c(n_lon, n_lat, n_group))
  lmf_all_group_loo = matrix(NA, nrow = n_lon, ncol = n_lat)
  
  slope_segment = array(NA, dim = c(n_lon, n_lat, n_breakpoints + 1))
  intercept_segment = array(NA, dim = c(n_lon, n_lat, n_breakpoints + 1))
  
  r_squared = matrix(NA, nrow = n_lon, ncol = n_lat)
  obs_mod_last_value = array(NA, dim = c(n_lon, n_lat, 2))
  
  for (i in seq_len(n_lon)) {
    for (j in seq_len(n_lat)) {
      if (is.na(lmf_total[i, j, 1])) next#i=120;j=40
      
      response_data = data.frame(
        y = lmf_total[i, j, ],
        x = tas_change[i, j, ])
      
      pred_result = fit_segmented_prediction(
        response_data            = response_data,
        co2_model_each_group_loo = co2_model_each_group_loo,
        co2_model_present_day    = co2_model_present_day,
        co2_model_all_group_loo  = co2_model_all_group_loo,
        n_breakpoints            = n_breakpoints
      )
      
      if (is.null(pred_result)) next
      
      lmf_group_loo     [i, j, ] = pred_result$lmf_group_loo
      lmf_all_group_loo [i, j]   = pred_result$lmf_all_group_loo
      slope_segment     [i, j, ] = pred_result$slope
      intercept_segment [i, j, ] = pred_result$intercept
      r_squared         [i, j]   = pred_result$r_squared
      obs_mod_last_value[i, j, ] = pred_result$obs_mod_last_value
    }
  }
  
  list(
    lmf_group_loo_i    = lmf_group_loo,
    lmf_all_group_loo  = lmf_all_group_loo,
    slope_segment      = slope_segment,
    intercept_segment  = intercept_segment,
    mod_r_validation   = r_squared,
    obs_mod_last_value = obs_mod_last_value,
    lmf_sm             = lmf_fit$lmf_sm,
    lmf_ta             = lmf_fit$lmf_ta
  )
}

#============================================================
# 3.3 Main prediction function
#============================================================

cal_prediction_all = function(
    fp_fit_datain = "D:/23Attribution_onlyco2/0results/2calculate_LMF_only_summer/dry",
    fpout = "/public1/home/m8s001859/0jkl/22ac/0results/0fit_compound_dry_hot/cmip6_mpi",
    dasta = 1,
    daend = 5,
    core = 5,
    co2_model_each_group_loo,
    co2_model_present_day,
    co2_model_all_group_loo,
    co2_dis = 150,
    n_breakpoints = 1
) {
  
  input_files = file.path(
    fp_fit_datain,
    list.files(fp_fit_datain)[dasta:daend]
  )
  
  if (length(input_files) == 0) {
    stop("No input NetCDF files found.")
  }
  
  lmf_fit_all = list()
  for (k in 1:length(input_files)) {#k=2
    fcon = nc_open(input_files[k])
    fda = ncvar_get(fcon, 'NA')
    nc_close(fcon)
    
    fda[fda == 9999] = NA
    #lmf_sm_m   = fda[,,(co2_dis*1+1):(co2_dis*2)];lmf_ta_m   = fda[,,(co2_dis*2+1):(co2_dis*3)];lmf_sm_v   = fda[,,(co2_dis*3+1):(co2_dis*4)];lmf_ta_v   = fda[,,(co2_dis*4+1):(co2_dis*5)]
    
    lmf_tot    = fda[ , , (co2_dis * 0 + 1):(co2_dis * 1)]
    tas_change = fda[ , , (co2_dis * 5 + 1):(co2_dis * 6)]    
    
    lmf_sm = fda[ , , (co2_dis * 2)] * fda[ , , (co2_dis * 4)]
    lmf_ta = fda[ , , (co2_dis * 3)] * fda[ , , (co2_dis * 5)]
    
    lmf_fit_all[[k]] = list(lmf_tot = lmf_tot, 
                            tas_change = tas_change,
                            lmf_sm = lmf_sm,
                            lmf_ta = lmf_ta)
  }
  gc()
  
  
  cl = makeCluster(core)
  registerDoParallel(cl)
  prediction_all <- foreach(
    k = seq_along(input_files),
    .packages = c("segmented", "foreach", "ncdf4"),
    .export = c(
      "fit_segmented_prediction",
      "run_segmented_prediction_grid"
    )
  ) %dopar% {
    
    run_segmented_prediction_grid(
      co2_model_each_group_loo = co2_model_each_group_loo,
      co2_model_present_day    = co2_model_present_day,
      co2_model_all_group_loo  = co2_model_all_group_loo,
      lmf_fit                  = lmf_fit_all[[k]],
      n_breakpoints            = n_breakpoints
    )
  }
  
  stopCluster(cl)
  
  for (k in seq_along(input_files)) {#k = 1
    
    file_parts = strsplit(basename(input_files[k]), "_")[[1]]
    start_id = which(grepl("times", file_parts))
    end_id   = which(grepl("yes|no", file_parts))
    sample_name = paste(file_parts[start_id:end_id], collapse = "_")
    modoutna = paste('pre_', paste(file_parts[3:5], collapse = "_"), sep = '')
    threshold = file_parts[8]
    
    output_file = paste0(
      modoutna,
      "_1950_2023_",
      threshold,
      "_segmented",
      n_breakpoints,
      "_",
      sample_name
    )
    
    output_path = file.path(fpout, output_file)
    saveRDS(
      prediction_all[[k]],
      file = output_path,
      compress = "xz"
    )
  }
}


#
#		would it be better to load libraries here
#					or in the "_run"-file?
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#			functions to read meteorological data 
#				of a number popular sources
# some of them (Roshydromet-Obnins archive, e.g.) require registration;
# that is why data are supposed to be stored on a local server before to be processed
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#********************************************
# designed to read and transform KNMI-Climate Explorer data
ReadInd <- function(dir_name, file_name) {	
	read_df <- read.table(file.path(dir_name, file_name), 
							stringsAsFactors = FALSE, header = FALSE,
							na.strings = "NA",
							comment.char = "#")

	if (ncol(read_df) == 2) {
		res_df <- TrsfIndLng(df_inp = read_df)
	} else {
		 res_df <- TrsfIndWd(df_inp = read_df)
	}

	return(res_df)
}	

# transform long-formatted monthly data to wide format
TrsfIndLng <- function(df_inp) {	

	colnames(df_inp) <- c("year", "dat")
	df_inp[df_inp < -900] <- NA 	# -999.9 may mean NA, too!

	months_names <- paste("month", 1:12, sep = "_")

	# month may be coded as a fraction part of a year
	# check-up of triviality https://stackoverflow.com/a/45273181/8465924
	if (all(round(df_inp$year) == df_inp$year)){ # all years are integers
		years_set <- unique(df_inp$year)
		months_col <- rep(months_names, times = length(years_set))
	} else { 									# non-trivial case: years have a fractional part
		month_ns <- round((df_inp$year %% 1) * 12, 1)
		if (any(month_ns < 1)) {
			month_ns <- month_ns + 1
		}
		months_col <- paste("month", month_ns, sep = "_")
		df_inp <- df_inp %>% mutate(year = floor(year))
	}
	df_inp <- cbind(df_inp, month_name = months_col)

	df_inp <- spread(data = df_inp, key = month_name, value = dat)
	df_inp <- df_inp %>% select_(.dots = c("year", months_names))

	return(df_inp)
}

# add colnames and calculate ann-av of the wide-formatted data
TrsfIndWd <- function(df_inp) {	
	colnames(df_inp) <- c("year", paste("month", 1:12, sep = "_"))
	df_inp[df_inp < -900] <- NA 	# -999.9 may mean NA, too!
	df_inp$annual <- rowSums(df_inp[, -1])
	return(df_inp)
}
#********************************************

#********************************************
#	work with Obninsk-Roshydromet data	

# read reference info 
# (that's a table with a header and "" sep)
ReadRefInfo <- function(dir_name, file_name) {
	ref_info_file_name <- paste0(dir_name, file_name)
	res_df <- read.csv(ref_info_file_name, stringsAsFactors = FALSE, sep = "")
	return(res_df)
}

# read monthly-aggregated data from Obninsk archive
ReadObninsk <- function(dir_name, file_name) {

	ref_info_file_name <- paste0(dir_name, file_name)
	res_df <- read.csv(ref_info_file_name, stringsAsFactors = FALSE,
		sep = "", header = FALSE)
	colnames(res_df) <- c("st_id", "year", paste("month", 1:12, sep = "_"))

	res_df$annual <- rowSums(res_df[, -c(1, 2)])

	return(res_df)
}

# read corrected presipitation data from Obninsk archive
ReadObninskCorrtdPrecip <- function(dir_name, file_name) {

	ref_info_file_name <- paste0(dir_name, file_name)
	res_df <- read.csv(ref_info_file_name, stringsAsFactors = FALSE,
		sep = ";", header = FALSE, na.strings = "9999.9")
	colnames(res_df) <- c("st_id", "year", 
		"month", "fullness_flag",
		"raw_sum", "crrtd_sum", "crrtd_fluid", "crrtd_solid",
		"crrtd_mix")

	# res_df$annual <- rowSums(res_df[, -c(1, 2)])

	return(res_df)
}

# @geo_string_vct is a vector of characters formatted like "62ТА 01?"
# @sep_mark_1, @sep_mark_2, @tail_mark are special symbols
# causing troubles by transforming @geo_string_vct to numeric
ExtractRoshydromCoord <- function(geo_string_vct,
	sep_mark_1, sep_mark_2, tail_mark)
{
	temp_lat_string <- unlist(strsplit(
		geo_string_vct,
		# regular expression with "[", "]" and fixed = FALSE  
		split = paste0("[", sep_mark_1, sep_mark_2, "|", 
			sep_mark_2, sep_mark_1, "]")
			)
		)
	temp_lat_string_2 <- unlist(strsplit(temp_lat_string, split = tail_mark,
			fixed = TRUE))
	str_vct_temp <- gsub(" ", "", temp_lat_string_2, fixed = TRUE)
	# odd elements are degrees, even ones are minutes
	return(as.numeric(str_vct_temp[seq(from = 1, to = length(str_vct_temp), by = 2)]) +
	as.numeric(str_vct_temp[seq(from = 2, to = length(str_vct_temp), by = 2)])/60)
}

# read Roshydromet meta-data
ExtractRoshydromMetaData <- function(dir_name, file_name) {

	meta_df <- read.csv2(paste0(dir_name, file_name), 
		header = FALSE, skip = 3,
		dec = ".", stringsAsFactors = FALSE)
	meta_df <- meta_df[, -1]
	names(meta_df) <- c("st_id", "st_name", "lat", "lon",
		"alt", "obs_begin", "comment")

	meta_df <- meta_df %>% filter(!is.na(meta_df$st_id))

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#					transform degrees to decimals
	#	the problem is specific symbols inside the coordinates' strings

	# extract specific symbols from strings of meta-information
	# !!! TODO: is data-specific as far
	smbl_splt <- c(rawToChar(charToRaw(meta_df[1, ]$lat)[3]),
		rawToChar(charToRaw(meta_df[1, ]$lat)[4]),
		rawToChar(charToRaw(meta_df[1, ]$lat)[7]))

	# meta_df$lat_dec
	meta_df$lat_dec <- sapply(
		FUN = function(i) {
			ExtractRoshydromCoord(geo_string_vct = meta_df[i, ]$lat,
				sep_mark_1 = smbl_splt[1], sep_mark_2 = smbl_splt[2], 
				tail_mark = smbl_splt[3])
		}, X = seq_along(meta_df$lat)
	)

	meta_df$lon_dec <- sapply(
		FUN = function(i) {
			ExtractRoshydromCoord(geo_string_vct = meta_df[i, ]$lon,
				sep_mark_1 = smbl_splt[1], sep_mark_2 = smbl_splt[2], 
				tail_mark = smbl_splt[3])
		}, X = seq_along(meta_df$lon)
	) 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	return(meta_df)
}

# TODO: pass a column name for calculations -> lazyeval:: could be helpful
AddAnnAv <- function(meteo_df){

	test_ann <- meteo_df %>% group_by(st_id, year) %>%
		summarise(annual = sum(crrtd_sum))

}
#********************************************

# quick abd dirty approach for geo-clustering
TransfCoordToRegion <- function(lat, lon) {
	if ((lat > 50) & (lon < 55)) {return("Center")}
	if ((lat < 50) & (lon < 60)) {return("South")}
	if ((lat < 65) & ((lon > 60) & (lon < 90))) {
		return("West Siberia")
	}
	if ((lat < 65) & ((lon > 90) & (lon < 130))) {
		return("Soth East Siberia")
	}
	if ((lat < 55) & (lon > 130)) {
		return("Soth Far East")
	}
	if ((lat > 55) & (lon > 130)) {
		return("Northern Far East")
	} 
	else {return(NA)}
}

# vectorised quick abd dirty approach for geo-clustering
# @meteo_df should have columns lat_dec & lon_dec
TransfCoordToRegion_Vct <- function(meteo_df) {
	res_df <- 	meteo_df %>% mutate(region = case_when(
				((lat_dec > 50) & (lon_dec < 55))~ "Center",
				((lat_dec < 50) & (lon_dec < 60))~ "South",	
				((lat_dec < 65) & ((lon_dec > 60) & (lon_dec < 90))) ~ "West Siberia",			
				((lat_dec < 65) & ((lon_dec > 90) & (lon_dec < 130))) ~ "Soth East Siberia",		
				((lat_dec < 55) & (lon_dec > 130)) ~ "Soth Far East", 		
				((lat_dec > 55) & (lon_dec > 130)) ~ "Northern Far East"
			)
		)
	return(res_df)
}

# @ref_df is the data frame as returned by ReadRefInfo()
# !!! special structure of the data is kept in mind (groupped by 10 days)
ReadRefData <- function(dir_name, info_df, i_ref_st) {

	data_df <- read.csv(paste0(dir_name, info_df$name[i_ref_st], ".csv"),
		skip = 1, stringsAsFactors = FALSE, sep = ";")
	data_df <- data_df[, 1:(1 + 3 * 12)]
	colnames(data_df) <- c("year", sapply(
		function(i) {paste("month", i, 1:3, sep = "_")}, X = 1:12
		)
	)
	data_df <- data_df[!is.na(data_df$year), ]
	data_df$annual <- sapply(function(i) {sum(as.numeric(data_df[i, -1]))}, 
		X = seq_along(data_df$year))
	data_df <- data_df %>% mutate(ann_normal = (annual - mean(annual)) /mean(annual))

	return(data_df)
}

## read atm indices in KNMI format (https://climexp.knmi.nl/selectindex.cgi?id=someone@somewhere)
# ReadInd <- function(dir_name, file_name, n_skip) {	
# 	res_df <- read.table(paste0(dir_name, file_name), 
# 							stringsAsFactors = FALSE, header = FALSE,
# 							skip = n_skip, na.strings = "-999.9000")
# 	colnames(res_df) <- c("year", paste("month", 1:12, sep = "_"))
# 	res_df[which(res_df < -900)] <- NA
# 	res_df$annual <- rowSums(res_df[, -1])
# 	return(res_df)
# }	
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# @info_df has structure as returned by ReadRefInfo()
# @rad_infl radius of influence, degrees
SelectInAreaInfo <- function(dir_name = dd_name, file_name = df_name,
	info_df, rad_infl, st_name, n_avbl_years) {
	st_info_file_name <- paste0(dir_name, file_name)
	# read data of GRDC stations
	st_info <- read.csv(st_info_file_name, stringsAsFactors = FALSE, sep =";",
		quote="", header = TRUE, row.names = NULL)

	i_ref_st <- which(info_df$name %in% st_name)
	if (length(i_ref_st) < 1) {stop(cat
		("No such station in the reference info: ",
			"\n\r", "given station name:", st_name, 
		"\n\r", sep = "")
		)
	}

	st_info$distance <- sqrt((info_df$Lat[i_ref_st] - st_info$lat)^2 +
		(info_df$Lon[i_ref_st] - st_info$lon)^2)

	st_in_area_info <- st_info %>% filter (m_yrs > n_avbl_years) %>%
		filter(distance < rad_infl)	

	if (length(st_in_area_info[, 1]) < 1) {stop(cat("No stations into the area",
		" with a given number of observation years", "\n\r", "Reference station ", 
		info_df$name[i_ref_st], "\b\r",
		"rad_infl = ", rad_infl, " , n_avbl_years = ", n_avbl_years, "\n\r",
		sep = ""))}

	return(st_in_area_info)
}

# @ meteo_df is the data frame with the "annual" column
# and the "st_id" column
LastMissed <- function(meteo_df) {

	res <- meteo_df %>% filter(year < 2016) %>%
		group_by(st_id) %>%
		summarise(
			last_missed_ind = last(which(is.na(annual))),
			last_missed_year = year[last(which(is.na(annual)))]
		) 

	return(res)	
}

ExtractSt <- function(meteo_df, station_id){

	missed_df <- LastMissed(meteo_df) %>%
		filter(st_id %in% station_id)
	# return(missed_df)	

	res <- meteo_df %>% filter(st_id %in% station_id) %>%
		filter(year > missed_df$last_missed_year)

	return(res)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		functions to work with hydrologic basins
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# @dir_name is the name of the data dir (with a shape-file)
# @file_name is the name of the data file (the shape-file, actually)
ReadBasin <- function(dir_name, file_name){
	map_res <- readOGR(paste0(sh_dir_name, sh_fl_name))	
	return(map_res)
}

# map_obj is the geo-object representing he basin
# @x and @y are the coordinates of the spatial points under consideration
InBasinCheck <- function (map_obj, x, y){

	# SpatialPoints() isn't happy with vectors -> matrix or data frame is needed
	test_point <- SpatialPoints(coords = cbind(x, y), 
		proj4string = CRS(proj4string(map_obj)))
	res <- sapply(FUN = function(i) gContains(map_Volga, test_point[i]), 
		X = seq_along(test_point))

	if (sum(res) < 1) {stop(cat("No stations into the basin",
		" with a given number of observation years",
		sep = ""))}

	return(res)
}

# plots a series of the points that may be in basin or may be not
# plots interactively; it seems, ggplot2 is a better way to store plot objs
# @x and @y are the coordinates of the spatial points under consideration
InBasinPlot <- function(map_obj, x, y){

	test_point <- SpatialPoints(coords = cbind(x, y), 
		proj4string = CRS(proj4string(map_obj)))
	inside_test <- sapply(FUN = function(i) gContains(map_Volga, test_point[i]), 
		X = seq_along(test_point))

	for (i in seq_along(test_point)) {
		# wrapper for the name of the plot
		plot_name <- paste0("Point ", i)
		dev.new()
		if (inside_test[i]) {
			plot(map_obj, col = "gray", b = "gray30", main = plot_name)
			points(test_point[i], pch = 21, col = "darkgreen", bg = "green")		
		} else {
			plot(map_obj, col = "gray", b = "gray30", main = plot_name)
			points(test_point[i], pch = 21, col = "darkred", bg = "red")
		}
	}
}

# plots conseuqently the reference point (in blue) and the cloud of points around them
# the around points may be in (in green) or out (in red) of the basin
# plots interactively; it seems, ggplot2 is a better way to store plot objs
# @map_obj is an geo-object representing the considered basin
# @x and @y are the coordinates of the spatial points under consideration
RelToBasinPlot <- function(map_obj, x_ref, y_ref,
	x, y, ref_name = "", sub_text = ""){

	ref_point <- SpatialPoints(coords = cbind(x_ref, y_ref), 
		proj4string = CRS(proj4string(map_obj)))	

	test_points <- SpatialPoints(coords = cbind(x, y), 
		proj4string = CRS(proj4string(map_obj)))

	inside_test <- sapply(FUN = function(i) gContains(map_Volga, test_points[i]), 
		X = seq_along(test_points))
	
	plot(map_obj, col = "gray", b = "gray30", main = ref_name, 
		sub = sub_text, cex.sub = 0.5)
	points(ref_point, pch = 21, col = "royalblue", bg = "darkblue")	
	for (i in seq_along(test_points)) {
		# wrapper for the name of the plot
		if (inside_test[i]) {
			points(test_points[i], pch = 21, col = "darkgreen", bg = "green")		
		} else {
			points(test_points[i], pch = 21, col = "darkred", bg = "red")
		}
	}
}	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# @i_obs_st is an index of the observation station to compare with the reference station
ExtractObsData <- function(dir_name, info_df, i_obs_st) {

	# data-specific consts
	dd_name_pref <- c("europe//", "asia//")
	dd_cont_code <- c(6, 2)
	df_name_pref <- "_Q_YVM"
	# 4-digit year; one or more whitespaces digits (possibly with a point inside); end of input which may ve precede with whitespaces
	match_pattern <- "^\\d\\d\\d\\d;[[:space:]]+\\d*\\.*\\d*;[[:space:]]*\\z"

	contin_code <- substr(info_df$grdc_no[i_obs_st], 
		start = 1, stop = 1)
	prefx <- dd_name_pref[which(dd_cont_code %in% contin_code)]

	file_to_read <- paste0(dir_name, prefx, info_df$grdc_no[i_obs_st], 
		df_name_pref, ".csv")	

	test_con <- file(file_to_read, "r", blocking = FALSE)
	test_data_file <- readLines(con = test_con, n = 1000)
	close(test_con)
	unlink(test_con)

	ann_ts_data_text <- test_data_file[str_detect(test_data_file, 
		match_pattern)]
	data_ann <- read.table(text = ann_ts_data_text, sep = ";", 
		stringsAsFactors = FALSE, na.strings = -999.000)

	data_ann <- data_ann[, -3] # last NA column appears due to ; before the end of the input
	colnames(data_ann) <- c("year", "MQ") # for consistency with GRDC 2013
	data_ann <- data_ann %>% mutate(MQ_norm = (MQ - mean(MQ))/mean(MQ))

	return(data_ann)
}

# @yrs_vct is a vector of observed years
# @yrs_to_match is a vector of years in reference data
# range(obs_data[[i]]$year
AllYearsCheck <- function(yrs_vct, yrs_to_match){
	# there is no problem if all years present
	super_cond <- (range(yrs_vct)[2] - range(yrs_vct)[1]) == (length(yrs_vct) - 1)
	if (super_cond) {
		return(yrs_vct)
	} else {
		cmmn_years <- dplyr::intersect(yrs_vct, yrs_to_match)
		yrs_diff_sqv <- yrs_vct[-length(yrs_vct)] - stats::lag(yrs_vct, 1)[-1]
		# as the simplest solution just cut all years before the last lost one
		i_last_lost <- max(which(yrs_diff_sqv < -1))
		return(yrs_vct[-seq(from = 1, to = i_last_lost, by = 1)])
	}
}

SmoothObsData <- function(data_df, col_name, n_ma = 11){

	dat <- data_df[, col_name]
	ma_col_name <- paste(col_name, 11, "ma", sep = "_")

	fn <- rep(1/n_ma, n_ma)
	# result of filter is of the ts type -> force to numerical
	data_df[, ma_col_name] <- as.numeric(stats::filter(dat, fn, sides = 1))

	return(data_df)
}

# @year_min_for_annot is the xmin to ajust text across the plot field
PlotTimeSerCompar <- function(ref_df, obs_list, 
	ref_col_name = "ann_normal",
	obs_col_name = "MQ_norm",
	plot_name = "", sub_name = "",
	annot_text = "",
	year_min_for_annot = 1850){

	dark_col <- c(ref = "darkred", obs = "darkgreen")
	line_col <- c(ref = "red", obs = "forestgreen")
	fill_col <- c(ref = "coral", obs = "green")

	# delete rows with NAs for smoothed values
	ref_df_end <- ref_df[!is.na(ref_df[, ref_col_name]), ]

	plot_res <- ggplot(data = ref_df_end, aes_string(
			x = "year", y = ref_col_name)) +
		geom_line(col = line_col["ref"]) + 
		geom_point(shape = 21, fill = fill_col["ref"], colour = dark_col["ref"]) +
		labs(ggtitle(label = plot_name, subtitle = sub_name))

		alpha_seq <- seq(from = 0.25, to = 0.75, length.out = length(obs_list))	
	for (i in seq_along(obs_list)) {

		obs_dat <- obs_list[[i]][!is.na(obs_list[[i]][, obs_col_name]), ]
		# return(obs_dat)

		plot_res <- plot_res +
			geom_line(inherit.aes = FALSE, data = obs_dat,
				aes_string(x = "year", y = obs_col_name), colour = line_col["obs"],
				alpha = alpha_seq[i]) +
			geom_point(inherit.aes = FALSE, data = obs_dat,
				aes_string(x = "year", y = obs_col_name),
				shape = 21, fill = fill_col["obs"], colour = dark_col["obs"]) 
	}

	plot_res <- plot_res +
			annotate(geom = "text", xmin = year_min_for_annot - 20, x = year_min_for_annot, 
				y = 0.5, label = annot_text, size = 2.5) +
			theme(plot.subtitle = element_text(size = 5, face="italic", color = "gray10"))

	return(plot_res)
}

# TODO How to set discrete color scales?
PlotTimeSerCompar_AdvCol <- function(ref_df, obs_list, plot_name = "",
	sub_name = ""){

	dark_col <- c(ref = "darkred", obs = "darkgreen")
	line_col <- c(ref = "red", obs = "forestgreen")
	fill_col <- c(ref = "coral", obs = "green")

	colour_count <- length(obs_list)
	getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

	plot_res <- ggplot(data = ref_df, aes(x = year, y = ann_normal)) +
		geom_line(col = line_col["ref"]) + 
		geom_point(shape = 21, fill = fill_col["ref"], colour = dark_col["ref"]) +
		labs(ggtitle(label = plot_name, subtitle = sub_name))

	for (i in seq_along(obs_list)) {
		alpha_seq <- seq(from = 0.25, to = 0.75, length.out = length(obs_list))	
		plot_res <- plot_res +
			geom_line(inherit.aes = FALSE, data = obs_list[[i]],
				aes(x = year, y = MQ_norm),
				alpha = alpha_seq[i], col = getPalette(i)) +
			scale_colour_manual(values = c(getPalette)) +
			geom_point(inherit.aes = FALSE, data = obs_list[[i]],
				aes(x = year, y = MQ_norm),
				shape = 21, fill = getPalette(i)) +
			theme(legend.position="right")
	}

	return(plot_res)
}

Correl <- function(ref_df, obs_df) {	
	years_for_cor <- dplyr::intersect(ref_df[, "year"], obs_df[, "year"])
	Q_ref <- ref_df[which(ref_df[, "year"] %in% years_for_cor), "ann_normal"]
	Q_data <- obs_df[which(obs_df[, "year"] %in% years_for_cor), "MQ_norm"]
	# return(list(Qr = Q_ref, Q_d = Q_data))
	if (length(Q_ref) == 0| length(Q_data) == 0) {
		stop(paste0("Not enough data for correlation", "\n\r"))
	}
	return(list(correl_value = cor(Q_ref, Q_data), 
				n_years = length(years_for_cor)))
}

ConservCheck <- function(meteo_df,
	slctd_st_id){
	slctd_st_data <- meteo_df %>% filter(st_id %in% slctd_st_id) %>%
			filter(!is.na(annual)) %>% 
			mutate(annual_norm = annual^(1/3)) %>%
			select(st_id, year, annual, annual_norm)
		
		if(length(slctd_st_data$annual_norm) > 1){
			# kr_res <- kruskal.test(test$annual_norm)
			# brt_res <- bartlett.test(test$annual_norm)
			neu_res <- bartels.test(slctd_st_data$annual_norm)
			snh_res <- snh.test(slctd_st_data$annual_norm)
			pt_res <- pettitt.test(slctd_st_data$annual_norm)
			br_res <- br.test(slctd_st_data$annual_norm)
		} else {
			neu_res <- NA
			snh_res <- NA
			pt_res <- NA 
			br_res <- NA
		}
		
		p_conserv <- min(neu_res[["p.value"]], snh_res[["p.value"]],
				pt_res[["p.value"]], br_res[["p.value"]])

		res <- (c(st_id = slctd_st_id,
			p_neu = neu_res[["p.value"]], p_snh = snh_res[["p.value"]],
			p_pt = pt_res[["p.value"]], p_br = br_res[["p.value"]],
			p_min = p_conserv))
	
		return(res)
}

ConstructDataForCorrel <- function(obs_df, atm_df, 
	year_col = "year", col_to_cogen = "annual") {

	# NB different defaults of tibble and dataframes for @drop
	# reqire to define @drop explicitly
	years_for_cor <- unlist(dplyr::intersect(obs_df[, year_col, drop = FALSE], 
			atm_df[, year_col, drop = FALSE]))

	# not very elegant playing with @drop argument is needed to
	# overcome difficulties of tibbles vs dataframes
	atm_data <- unlist(atm_df[which(unlist(atm_df[, year_col]) %in% years_for_cor), 
			col_to_cogen])
	param_data <- unlist(obs_df[which(unlist(obs_df[, year_col]) %in% years_for_cor), 
			col_to_cogen])
	
	return(data.frame(year = years_for_cor,
		atm_ind = atm_data, obs_data = param_data))

}
# write in a conbinient format under a reasonable file name
WriteData <- function(dat_df){
	file_name <- paste0(substitute(dat_df), "_dat.csv")
	write.table(file = file_name, x = dat_df,
	quote = FALSE, col.names = TRUE, row.names = FALSE,
	sep = ";")
}

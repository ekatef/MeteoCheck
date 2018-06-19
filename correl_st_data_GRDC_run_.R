# structure of the data files to be correlated with
# is out of date (corresponds to GRDC 2013)

rm(list = ls())

Sys.setenv(LANGUAGE='en')
# library("dplyr")
# library("readr")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

library(rgdal)	# for readOGR()
library(sp)		# for Spatial classes
library(rgeos)	# for gContains

library(stringr)
library(lubridate)

# should be set here, not in the config file
# the reason is rm() in the begin of the file
wd_name <- "..."
setwd(wd_name)

source("correl_st_data_GRDC_config.R")
source("correl_st_data_GRDC_fun.R")

names_of_months <- paste("month", 1:12, sep = "_")

ref_info <- ReadRefInfo()
print(ref_info)

i_st_testing <- 1				# for testing

res_file_name <- paste0(res_dir_name, ref_info$name, "_monthly", ".csv")
for (j in seq_along(ref_info$name)){
# for (j in i_st_testing){

	ref_data <- ReadRefData(info_df = ref_info, i_ref_st = j)
	res_df <- ref_data %>% select(year)
	
	for (i in 1:12) {
		month_name <- names_of_months[i]

		ann_av <- ref_data %>% select(starts_with(paste0(month_name, "_"))) %>%
			rowSums(.)

		month_norm <- mean(ann_av)
		ann_norm <- (ann_av - month_norm)/month_norm

		monthly_names <- paste0(month_name, c("_sum", "_norm"))
		res_df[, monthly_names] <- cbind(ann_av, ann_norm)

	}

		res_ann <- res_df %>% select(ends_with("_sum")) %>%
			mutate(ann = rowSums(.))

	dat_df_long <- res_df %>% select(year, ends_with("_sum")) %>%
			gather(key = date, value = runoff_value, position = - year) %>%
			mutate(month_n = str_extract(date, "\\d\\d*")) %>%
			mutate(full_date = ymd(paste(year, month_n, "01", sep = "-"))) %>%
			arrange(rev(desc(full_date))) %>%
			mutate(month = month(full_date), year = year(full_date)) %>%
			select(full_date, month, year, runoff_value)

	dat_long_check <- dat_df_long %>% mutate(year = year(full_date)) %>%
		group_by(year) %>% summarise(ann = sum(runoff_value))

	# print(str(dat_long_check))
	# print(paste(ref_info$name[j]))
	# print(identical(res_ann$ann, dat_long_check$ann))
	# print(all.equal(res_ann$ann, dat_long_check$ann))

	long_res_file_name <- paste0(res_dir_name, ref_info$name[j], "_long_ts.csv")
	write.table(file = long_res_file_name, dat_df_long,
		quote = FALSE, row.names = FALSE, sep = ";")

	dev.new()
	plot(x = dat_long_check$year, y = dat_long_check$ann, 
		type = "l", col = "red", lwd = 2,
		main = ref_info$name[j])

	write.table(file = long_res_file_name[j], dat_df_long,
		quote = FALSE, row.names = FALSE, sep = ";")
}

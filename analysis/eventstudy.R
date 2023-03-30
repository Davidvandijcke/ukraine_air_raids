

# Ukraine air raids
# DVD + ALW
# Ridgeline plots for diagnostics

# Load packages needed
packages_load <- c("tidyr", "tidyverse", "stringr",
      "readr", "dplyr", "ggridges", "ggplot2", "here", "data.table", "ggpubr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = packages_load, character.only = TRUE)
 
# set relative directory 
cdir <- here::here() # need to open from analysis folder
setwd(gsub("/code/analysis", "", cdir))

data_dir <- "data/final/estimates_filtered" # change t
data_dir_s <- "data/final/estimates_sample" # change t
data_dir_inter <- "data/final/estimates_interaction_filtered" # change t

out_dir <- "output/figs/estimates_filtered" # change to output directory
out_dir_s <- "output/figs/estimates_sample" # change to output directory
out_dir_inter <- "output/figs/estimates_interaction" # change to output directory
final_dir <- "data/final"

basicPrepData <- function(data) {
  data$month  <- substr(data$id, 6, 7)

  data_long <- gather(data, period, movement, e_TYPE_300:e_TYPE_1800, factor_key=TRUE)

  data_long$period <- gsub('e_TYPE_', '', data_long$period)
  data_long$period <- gsub('[.]', '-', data_long$period)

  data_long$period_num <- as.numeric(data_long$period)

  data_long$period_num <- (data_long$period_num)/60

  data_long$period_str <- as.character(data_long$period_num)
  return(data_long)
}

#======================
prepData <- function(data, by_time = TRUE) {

  data_long <- basicPrepData(data)


  sort(data_long$period, decreasing = FALSE)

  if (by_time) { 
    data_long$time_period = 0
    data_long$time_period[data_long$month=="03"|data_long$month=="04"] = 1
    data_long$time_period[data_long$month=="05"|data_long$month=="06"] = 2
    data_long$time_period[data_long$month=="07"|data_long$month=="08"|data_long$month=="09"] = 3
  } else {
    data_long$time_period <- 0
  }

  plot_data <- data_long %>% 
  group_by(period_num, time_period) %>%
  summarise(x_mean = mean(movement), x_sd = sd(movement))

  plot_data$max_tick <- plot_data$x_mean + 2*plot_data$x_sd    


  plot_data$min_tick <- plot_data$x_mean - 2*plot_data$x_sd     

  return(plot_data)        

  return(data_long)
}


#======================
mergeRegionVariables <- function(fn){


  # load data
  econ2 <- fread(file.path("data", "final", "econ_normalized.csv"), drop = "dwells_shops_jan")
  econ2[, dwells_shops := dwells_shops * caid_unique] # undo sample scaling so we can aggregate


  # get name_1 names
  names <- fread(file.path("data", "final", "econ.csv"), select = c("name_1", "name_2"))
  econ2 <- merge(econ2, names[!duplicated(name_2)], by = "name_2")

  econ2[, day := as.Date(day)]

  # balance econ2 data by day
  econ2 <- econ2[day %in% unique(econ2$day),]


  # merge level 1 and level 2 regions with data sequentially
  sdcols <- c("dwells_shops", "caid_unique")
  # aggregate econ2 sdcols columns by region_english and day 
  econ1 <- econ2[, lapply(.SD, sum), by = c("name_1", "day"), .SDcols = sdcols]  

  # normalize dwells_shops_ma30 with respect to january region-specific average
  normalizeShops <- function(econ, region_name) {
    # calculate january average by region_name
    econ[, dwells_shops := dwells_shops / caid_unique]
    econ[is.na(dwells_shops), dwells_shops := 0]
    temp <- econ[month(day) == 1, .(dwells_shops_jan = mean(dwells_shops, na.rm = TRUE)), by = region_name]
    econ <- merge(econ, temp, by = region_name)
    econ[, dwells_shops_ma30 := frollmean(dwells_shops, 15, fill = NA, align = "right", na.rm = TRUE), by = region_name]
    econ[, dwells_shops_ma30_jan_norm := dwells_shops_ma30/dwells_shops_jan]
    return(econ)
  }
  econ1 <- normalizeShops(econ1, "name_1")
  econ2 <- normalizeShops(econ2, "name_2")
  

  # load air_raid_alerts.csv to match name_1 and name_2 with region_english
  alerts <- fread(file.path("data", "final", "air_raid_alerts.csv"), select = c("NAME_1", "region_english"))
  alerts <- alerts[!duplicated(NAME_1)]
  setnames(alerts, "NAME_1", "name_1")
  econ1 <- merge(econ1, alerts, by = "name_1", all.x = TRUE)
  setnames(alerts, "name_1", "name_2")
  econ2 <- merge(econ2, alerts, by = "name_2", all.x = TRUE)

  # drop region_english rows with NAs in econ1 and econ2
  econ1 <- econ1[!is.na(region_english)]
  econ2 <- econ2[!is.na(region_english)]


  # load coefficients 
  data <- fread(file.path(data_dir,fn))

  # assign identifier variables
  data$month  <- substr(data$id, 6, 7)
  data[, region_english := substr(id, 20, nchar(id))]
  data[, day := as.Date(substr(id, 1, 20))]



  # merge data first with econ1, then with econ2, and fill in missing values from each
  data <- merge(data, econ1, by = c("day", "region_english"), all.x = TRUE)
  data <- merge(data, econ2, by = c("day", "region_english"), all.x = TRUE, suffixes = c("", ".x"))
  data[, dwells_shops := ifelse(is.na(dwells_shops), dwells_shops.x, dwells_shops)]
  data[, dwells_shops_ma30 := ifelse(is.na(dwells_shops_ma30), dwells_shops_ma30.x, dwells_shops_ma30)]
  data[, dwells_shops_ma30_jan_norm := ifelse(is.na(dwells_shops_ma30_jan_norm), dwells_shops_ma30_jan_norm.x, dwells_shops_ma30_jan_norm)]
  data <- data[, -c("dwells_shops.x", "dwells_shops_ma30.x", "dwells_shops_ma30_jan_norm.x")]


  # merge in regional variables (21 indicates 3-week moving average)
  region <- fread(file.path("data", "final", "region_variables_21.csv"))
  data <- merge(data, region, by = c("day", "region_english"))


  data_long <- gather(data, period, movement, e_TYPE_300:e_TYPE_1800, factor_key=TRUE)

  data_long$period <- gsub('e_TYPE_', '', data_long$period)
  data_long$period <- gsub('[.]', '-', data_long$period)

  data_long$period_num <- as.numeric(data_long$period)

  data_long$period_num <- (data_long$period_num)/60

  data_long$period_str <- as.character(data_long$period_num)


  data_long$time_period = 0
  data_long$time_period[data_long$month=="03"|data_long$month=="04"] = 1
  data_long$time_period[data_long$month=="05"|data_long$month=="06"] = 2
  data_long$time_period[data_long$month=="07"|data_long$month=="08"|data_long$month=="09"] = 3


  # join false alarms
  panel_join_wlags <- fread(file.path(final_dir, "panel_join_wlags.csv"))
  panel_join_wlags[, false_alarm_rate := false_alarm/(false_alarm + true_alarm)]
  setnames(panel_join_wlags, "date_ymd", "day")

  # moving average of false alarm share
  data_long <- merge(data_long, panel_join_wlags, by = c("region_english", "day")) %>% setDT()
    data_long$abv_med_falserate <- ifelse(data_long$false_alarm_rate>.769, 1, 0)
  
  
  return(data_long)

}


#======================
# econ by time period 
splitByTime <- function(data_long, y_var, split_var, legend, facet = NULL, 
                facet_labels =  c("March/April", "May/June", "July/August/September"), 
                labels = c("< Mdn", "> Mdn"), grp = NULL) {

  if (!is.null(facet)) {
    data_long[, time_id := get(facet)]
  } else {
    data_long[, time_id := 1] # all observations same group for calculating average
  }

  data_long[, split_var_mean := median(get(split_var), na.rm = TRUE), by = facet] # TODO
  data_long[, split_var_above_mean := get(split_var) > split_var_mean, by = c("time_period")]

  data_long[, split_var_above_mean := get(split_var) > split_var_mean, by = c("time_period")]

  data_long[, group := split_var_above_mean ]

  if (!is.null(grp)){
    data_long[, group := get(grp)]
  }


  # econ by exposure
  plot_data <- data_long %>% 
    group_by(period_num, group, time_id) %>%
    summarise(x_mean = mean(movement), x_sd = sd(movement))

  plot_data$max_tick <- plot_data$x_mean + 2*plot_data$x_sd            
  plot_data$min_tick <- plot_data$x_mean - 2*plot_data$x_sd          


  # remove rows with missing observations for group
  plot_data <- plot_data[!is.na(plot_data$group),]

  p <- ggplot(plot_data, aes(x=period_num, y=x_mean, group=group, colour = factor(group, labels = labels))) + 
    geom_point(alpha = 0.5) +
    theme_minimal() +theme(text=element_text(size=21)) +
    ylab(paste0("Movement ", y_var)) + 
    xlab("Minutes after alert") + 
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(color = legend) + 
    stat_smooth(aes(fill=factor(group)), data = plot_data[plot_data$period_num<0,], method = "loess", formula = 'y ~ x')  +
    stat_smooth(aes(fill=factor(group)), data = plot_data[plot_data$period_num>=0,], method = "loess", formula = 'y ~ x') +
    scale_color_manual(values=c("#009392", "#d0587e"))+
    scale_fill_manual(values=c("#009392", "#d0587e"), guide="none")+
    guides(color=guide_legend(override.aes=list(fill=c("#009392",  "#d0587e")))) + 
    theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) 
    # ggtitle("Changes in population movement after Ukrainean air raid alert") +

  if (!is.null(facet)){
    p <- p + facet_wrap(~ factor(time_id, labels = facet_labels))
  }

  return(p)
}





basicPlot <- function(plot_data, y_var){
  p <- ggplot(plot_data, aes(x=period_num, y=x_mean, group=time_period, color = factor(time_period, labels=c("March, April", "May, June", "July, Aug, Sept")))) + 
  geom_point(alpha = 1/2,)  +
  scale_color_manual(values=c("#009392", "#e5b9ad", "#d0587e"))+
  scale_fill_manual(values=c("#009392", "#e5b9ad", "#d0587e"), guide="none")+
  theme_minimal() +theme(text=element_text(size=21)) +
  ylab(paste0("Movement, ", y_var)) + 
  xlab("Minutes after alert") + 
  # ggtitle("Changes in population movement after Ukrainean air raid alert") +
  labs(color = "Period of conflict")+
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # decrease axis label size
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  stat_smooth(aes(fill=factor(time_period)), data = plot_data[plot_data$period_num<0,], method = "loess", formula = 'y ~ x')  +
  stat_smooth(aes(fill=factor(time_period)), data = plot_data[plot_data$period_num>=0,], method = "loess", formula = 'y ~ x') +
  # remove the color background of the legend color lines
  guides(color=guide_legend(override.aes=list(fill=c("#009392", "#e5b9ad", "#d0587e")))) 

  return(p)

}

pooledPlot <- function(plot_data, y_var){
  clr <-  '#2A7CAD'
  p <- ggplot(plot_data, aes(x=period_num, y=x_mean), color = clr) + 
  geom_point(alpha = 1/2, color = clr)  +
  theme_minimal() +theme(text=element_text(size=21)) +
  ylab(paste0("Movement, ", y_var)) + 
  xlab("Minutes after alert") + 
  # ggtitle("Changes in population movement after Ukrainean air raid alert") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # decrease axis label size
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  stat_smooth(fill = clr, color = clr, data = plot_data[plot_data$period_num<0,], method = "loess", formula = 'y ~ x')  +
  stat_smooth(fill = clr, color = clr, data = plot_data[plot_data$period_num>=0,], method = "loess", formula = 'y ~ x') +
  theme(legend.position="none") 
  # remove the color background of the legend color lines
  return(p)

}




###################################################################################

##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##
#### Base plots #### 
##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~##

library(ggpubr)

#==================================================
# pooling periods: xy movement
#==================================================

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path(data_dir, paste0("coeffs_", coeff, ".csv"))), by_time = FALSE)
  # plot
  p <- pooledPlot(plot_data, coeffs_list[[coeff]])

}

coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_pooled_mvmt <- lapply(names(coeffs_list), makePlot)

ggarrange(plotlist = plots_pooled_mvmt, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_distance_sum.pdf")), width = 15, height = 5)

#==================================================
# split by period: xy movement
#==================================================

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path(data_dir, paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}

coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

ggarrange(plotlist = plots_split_mvmt, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_period_distance_sum.pdf")), width = 15, height = 5)

#==================================================
# split by period: z movement
#==================================================

coeffs_list <- list("altitude_avg_" = "altitude (h)",
                    "altitude_avg_log" = "altitude (asinh)")
plots_split_altitude <- lapply(names(coeffs_list), makePlot)

ggarrange(plotlist = plots_split_altitude, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_period_altitude.pdf")), width = 15, height = 5)

#==================================================
## Split exposure by time
#==================================================

makePlot <- function(coeff, split_var){
  # plot
  data_long <- mergeRegionVariables(paste0("coeffs_", coeff, ".csv"))
  p <- splitByTime(data_long, coeffs_list[[coeff]], split_var = split_var,
      legend = "Alert exposure (moving average)", facet = "time_period")
  # plot

}
# distance sum 
coeffs_list <- list("distance_sum_" = "distance sum (m)")
plots <- lapply(names(coeffs_list), function(i) lapply(c("duration_ma", "duration_cum"), function(j) makePlot(i,j)))

  plots_split_mvmt_duration <- unlist(plots,recursive=FALSE)

ggarrange(plotlist = plots_split_mvmt_duration, nrow =2, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_exposure_sum.pdf")), width = 15, height = 9)

#==================================================
## Export figure fragments
#==================================================

plots_pooled_mvmt[1] 
  ggsave(file.path(out_dir, paste0("coeffs_pooled_mvmt_levels.pdf")), width = 12, height = 6)
  plots_pooled_mvmt[2] 
  ggsave(file.path(out_dir, paste0("coeffs_pooled_mvmt_asinh.pdf")), width = 12, height = 6)
  
plots_split_mvmt[1]
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_levels.pdf")), width = 12, height = 6)
plots_split_mvmt[2]
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_asinh.pdf")), width = 12, height = 6)

plots_split_altitude[1] 
ggsave(file.path(out_dir, paste0("coeffs_split_altitude_levels.pdf")), width = 12, height = 6)
plots_split_altitude[2] 
ggsave(file.path(out_dir, paste0("coeffs_split_altitude_asinh.pdf")), width = 12, height = 6)

#plots_split_mvmt_duration[1]
ggarrange(plotlist = plots_split_mvmt_duration[1], nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_duration_ma.pdf")), width = 15, height = 5)

#plots_split_mvmt_duration[2]
ggarrange(plotlist = plots_split_mvmt_duration[2], nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_duration_cum.pdf")), width = 15, height = 5)

#ggarrange(plotlist = c(plots_pooled_mvmt[1], plots_split_mvmt[1], plots_split_altitude[1], plots_split_mvmt_duration[1]), ncol = 3, nrow =2, common.legend = TRUE, legend = "bottom")



#==================================================
# urban vs. rural: split by period
#==================================================


## urban

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path("data/final/estimates_filtered_urban", paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}


coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

plots_split_mvmt[1]
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_levels_urban.pdf")), width = 12, height = 6)


## rural
makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path("data/final/estimates_filtered_rural", paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}


coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

plots_split_mvmt[1]
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_levels_rural.pdf")), width = 12, height = 6)



#==================================================
# extended time window
#==================================================

## 60

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path("data/final/estimates_filtered_60", paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}


coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

p1 <- plots_split_mvmt[1]
p1[[1]] <- p1[[1]] + ggtitle("60 minutes") + theme(plot.title = element_text(size = 18, hjust = 0.5))



## 240

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path("data/final/estimates_filtered_240", paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}


coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

p2 <- plots_split_mvmt[1]
p2[[1]] <- p2[[1]] + ylim(-500,2500) + ggtitle("240 minutes") + theme(plot.title = element_text(size = 18, hjust = 0.5))
ggarrange(plotlist = c(p1,p2), nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_levels_extended_window.pdf")), width = 12, height = 6)



#==================================================
# two-wall rule: away from home
#==================================================

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path("data/final/estimates_filtered_away_from_home", paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}


coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

p1 <- plots_split_mvmt[1] 
p1[[1]] <- p1[[1]] + ggtitle("Away from home") + theme(plot.title = element_text(size = 18, hjust = 0.5))

#==================================================
# two-wall rule: on the move (300m / hour)
#==================================================

makePlot <- function(coeff){
  # plot
  plot_data <- prepData(read.csv(file.path("data/final/estimates_filtered_moving", paste0("coeffs_", coeff, ".csv"))))
  # plot
  p <- basicPlot(plot_data, coeffs_list[[coeff]])
  
}


coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")
plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

p2 <- plots_split_mvmt[1]
p2[[1]] <- p2[[1]] + ggtitle("On the move") + theme(plot.title = element_text(size = 18, hjust = 0.5))


ggarrange(plotlist = c(p1,p2), nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_split_mvmt_levels_twowall.pdf")), width = 12, height = 6)



#==================================================
#### Split false positive rate by time #### 
#==================================================

table(data_long$fa)

makePlot <- function(coeff, split_var){
  # plot
  data_long <- mergeRegionVariables(paste0("coeffs_", coeff, ".csv"))
  data_long[, false_3 := frollmean(false_alarm_rate, align = "right", n = 3)]
  data_long[, false_7 := frollmean(false_alarm_rate, align = "right", n = 7)]
  
  temp <- data_long[day <= as.Date("2022-04-15"), .(false_first = mean(false_alarm_rate, na.rm = TRUE)), by = "region_english"]
  data_long <- temp[data_long, on = "region_english"]
  
  p <- splitByTime(data_long, coeffs_list[[coeff]], split_var = split_var,
                   legend = "False alarm rate", facet = "time_period")
  # plot
  
}
# distance sum 
coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")

plots <- lapply(names(coeffs_list), function(i) lapply(c("false_alarm_rate"), function(j) makePlot(i,j)))

plots_split_mvmt_falsealarms <- unlist(plots,recursive=FALSE)

ggarrange(plotlist = plots_split_mvmt_falsealarms, nrow =2, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_falsealarms_sum.pdf")), width = 15, height = 9)


ggarrange(plotlist = plots_split_mvmt_falsealarms[1], nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_period_falsepositiverate.pdf")), width = 15, height = 5)



# 3-day moving average
plots <- makePlot("distance_sum_", "false_3")
ggarrange(plots, nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_period_falsepositiverate_3.pdf")), width = 15, height = 5)


# 7-day moving average
plots <- makePlot("distance_sum_", "false_7")
ggarrange(plots, nrow =1, common.legend = TRUE, legend = "bottom")
ggsave(file.path(out_dir, paste0("coeffs_period_falsepositiverate_7.pdf")), width = 15, height = 5)


# first month of alarms
plots <- makePlot("distance_sum_", "false_first")
ggarrange(plots, nrow =1, common.legend = TRUE, legend = "bottom")



#==================================================
# split by period: xy movement w consistent sample
#==================================================


prepData_s <- function(data, by_time = TRUE) {
  
  data_long <- basicPrepData(data)
  
  sort(data_long$period, decreasing = FALSE)
  
  if (by_time) { 
    data_long$time_period = 0
    data_long$time_period[data_long$month=="03"|data_long$month=="04"] = 1
    data_long$time_period[data_long$month=="05"|data_long$month=="06"] = 2
    data_long$time_period[data_long$month=="07"|data_long$month=="08"|data_long$month=="09"] = 3
  } else {
    data_long$time_period <- 0
  }

  data_long$movement <- ifelse(data_long$movement==0,"NA",data_long$movement)
    
  plot_data <- data_long %>% 
    group_by(period_num, month) %>%
    summarise(x_mean = mean(movement), x_sd = sd(movement))
  
  plot_data$max_tick <- plot_data$x_mean + 2*plot_data$x_sd    
  
  plot_data$min_tick <- plot_data$x_mean - 2*plot_data$x_sd     
  
  data_keep <- data_long 
  
  return(plot_data)        
  
  return(data_long)
}

  makePlot <- function(coeff){
    # plot
    plot_data <- prepData_s(read.csv(file.path(data_dir, paste0("coeffs_", coeff, ".csv"))))
    # plot
    p <- basicPlot(plot_data, coeffs_list[[coeff]])
    
  }

  basicPlot_s <- function(plot_data, y_var){
    p <- ggplot(plot_data, aes(x=period_num, y=x_mean, group=month, color = factor(month, labels=c("March", "April", "May", "June", "July", "Aug", "Sept")))) + 
      geom_point(alpha = 1/2,)  +
      scale_color_manual(values=c("#009392", "#e5b9ad", "#d0587e"))+
      scale_fill_manual(values=c("#009392", "#e5b9ad", "#d0587e"), guide="none")+
      theme_minimal() +theme(text=element_text(size=21)) +
      ylab(paste0("Movement, ", y_var)) + 
      xlab("Minutes after alert") + 
      # ggtitle("Changes in population movement after Ukrainean air raid alert") +
      labs(color = "Period of conflict")+
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      # decrease legend title size and legend text size
      #theme(legend.title = element_text(size = 18), legend.text = element_text(size = 18)) +
      # decrease axis label size
      theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
      stat_smooth(aes(fill=factor(month)), data = plot_data[plot_data$period_num<0,], method = "loess", formula = 'y ~ x')  +
      stat_smooth(aes(fill=factor(month)), data = plot_data[plot_data$period_num>=0,], method = "loess", formula = 'y ~ x') +
      # remove the color background of the legend color lines
      guides(color=guide_legend(override.aes=list(fill=c("#009392", "#e5b9ad", "#d0587e")))) 
    
    return(p)
    
  }
  
coeffs_list <- list("distance_sum_" = "distance sum (m)", "distance_sum_log" = "distance sum (asinh)")

plots_split_mvmt <- lapply(names(coeffs_list), makePlot)

plots_split_mvmt[1]

ggsave(file.path(out_dir_s, paste0("coeffs_split_mvmt_levels_sample.pdf")), width = 12, height = 6)

##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##
##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##
##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##~~~~~~~~~~~##

library(lubridate)

data_long <- read.csv(file.path(data_dir_s, paste0("coeffs_", "distance_sum_", ".csv")))

#Gather data
data_long <- gather(data_long, period, movement, e_TYPE_300:e_TYPE_1800, factor_key=TRUE)

data_long$period <- gsub('e_TYPE_', '', data_long$period)
data_long$period <- gsub('[.]', '-', data_long$period)

data_long$period_num <- as.numeric(data_long$period)
data_long$period_num <- (data_long$period_num)/60
data_long$period_str <- as.character(data_long$period_num)

#manipulating dates for join
data_long$year  <- substr(data_long$id, 1, 4)
data_long$day  <- substr(data_long$id, 9, 10)
data_long$month  <- substr(data_long$id, 6, 7)
data_long$date <- paste(data_long$year, data_long$month, data_long$day)
data_long$date_ymd <- ymd(data_long$date)
data_long$region_english  <- substr(data_long$id, 20, 50)
data_long$region_english <- trimws(data_long$region_english)

data_long$time_period[data_long$month=="03"|data_long$month=="04"] = 1
data_long$time_period[data_long$month=="05"|data_long$month=="06"] = 2
data_long$time_period[data_long$month=="07"|data_long$month=="08"|data_long$month=="09"] = 3

data_long$movement <- ifelse(data_long$movement==0,NA,data_long$movement)

data_long <- na.omit(data_long)

plot_data_sample <- data_long %>% 
  group_by(period_num, time_period) %>%
  summarise(x_mean = mean(movement), x_sd = sd(movement))

ggplot(plot_data_sample, aes(x=period_num, y=x_mean, group=time_period, color = factor(time_period, labels=c("March, April", "May, June", "July, Aug, Sept")))) + 
  geom_jitter(alpha = 1/2,)  +
  scale_color_manual(values=c("#009392", "#e5b9ad", "#d0587e"))+
  scale_fill_manual(values=c("#009392", "#e5b9ad", "#d0587e"), guide="none")+
  theme_minimal() +theme(text=element_text(size=21)) +
  ylab("Movement (distance sum, m)") + 
  xlab("Minutes after alert") + 
  # ggtitle("Changes in population movement after Ukrainean air raid alert") +
  labs(color = "Period of Conflict")+
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_smooth(aes(fill=factor(time_period)), data = plot_data_sample[plot_data_sample$period_num<0,], method = "loess", formula = 'y ~ x')  +
  stat_smooth(aes(fill=factor(time_period)), data = plot_data_sample[plot_data_sample$period_num>=0,], method = "loess", formula = 'y ~ x') +
  # remove the color background of the legend color lines
  guides(color=guide_legend(override.aes=list(fill=c("#009392", "#e5b9ad", "#d0587e"))))


library(lubridate)

data_long <- read.csv(file.path(data_dir, paste0("coeffs_", "distance_sum_", ".csv")))

#Gather data
data_long <- gather(data_long, period, movement, e_TYPE_300:e_TYPE_1800, factor_key=TRUE)

data_long$period <- gsub('e_TYPE_', '', data_long$period)
data_long$period <- gsub('[.]', '-', data_long$period)

data_long$period_num <- as.numeric(data_long$period)
data_long$period_num <- (data_long$period_num)/60
data_long$period_str <- as.character(data_long$period_num)

#manipulating dates for join
data_long$year  <- substr(data_long$id, 1, 4)
data_long$day  <- substr(data_long$id, 9, 10)
data_long$month  <- substr(data_long$id, 6, 7)
data_long$date <- paste(data_long$year, data_long$month, data_long$day)
data_long$date_ymd <- ymd(data_long$date)
data_long$region_english  <- substr(data_long$id, 20, 50)
data_long$region_english <- trimws(data_long$region_english)

data_long$movement <- ifelse(data_long$movement==0,NA,data_long$movement)

data_long <- na.omit(data_long)

plot_data_sample <- data_long %>% 
  group_by(period_num, month) %>%
  summarise(x_mean = mean(movement), x_sd = sd(movement))

ggplot(plot_data_sample, aes(x=period_num, y=x_mean, group=month, color = factor(month, labels=c("March", "April", "May", "June", "July", "Aug", "Sept")))) + 
  geom_jitter(alpha = 1/2,)  +
  scale_color_manual(values=c("#009392", "#e5b9ad", "#d0587e", "#009392", "#e5b9ad", "#d0587e", "#009392"))+
  scale_fill_manual(values=c("#009392", "#e5b9ad", "#d0587e", "#009392", "#e5b9ad", "#d0587e", "#009392"), guide="none")+
  theme_minimal() +theme(text=element_text(size=21)) +
  ylab("Movement (distance sum, m)") + 
  xlab("Minutes after alert") + 
  # ggtitle("Changes in population movement after Ukrainean air raid alert") +
  labs(color = "Period of Conflict")+
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_smooth(aes(fill=factor(month)), data = plot_data_sample[plot_data_sample$period_num<0,], method = "loess", formula = 'y ~ x')  +
  stat_smooth(aes(fill=factor(month)), data = plot_data_sample[plot_data_sample$period_num>=0,], method = "loess", formula = 'y ~ x')  #+
  # remove the color background of the legend color lines
  guides(color=guide_legend(override.aes=list(fill=c("#009392", "#e5b9ad", "#d0587e", "#009392", "#e5b9ad", "#d0587e", "#009392"))))


  
  
  #==================================================
  # split by distance from frontline, simple post term
  #==================================================
  data_dir_post <- "data/final/post_filtered" # change t
    
  data_long <- read.csv(file.path(data_dir_post, paste0("coeffs_", "distance_sum_", ".csv"))) %>% setDT()
  
  median(data_long$post_inter)
  
  
  
  #==================================================
  # calculate coefficient on front line interaction
  #==================================================
  
  out_dir_inter <- "output/figs/estimates_interaction_filtered" # change to output directory
  
  
  library(lubridate)
  
  data_long <- read.csv(file.path(data_dir_inter, paste0("coeffs_", "distance_sum_", ".csv"))) %>% setDT()
  
  data_long$month  <- substr(data_long$id, 6, 7)
  
  
  cols <- colnames(data_long)
  cols <- cols[grepl("TYPE", cols)]
  distcols <- cols[grepl("distance", cols)]
  basecols <- setdiff(cols, distcols)
  
  data_long <- gather(data_long[,- ..basecols], period, movement,all_of(distcols),  factor_key=TRUE) %>% setDT()
  
  
  
  data_long[,period := as.character(period)]
  data_long$conf <- substr(data_long$period, nchar(data_long$period)-4, nchar(data_long$period))
  data_long$period <- gsub('e_TYPE_|...distance_to_front.|X.', '', data_long$period)
  data_long$period <- gsub('[.]', '-', data_long$period)
  
  data_long$period_num <- as.numeric(data_long$period)
  data_long$period_num <- (data_long$period_num)/60
  
  data_long$time_period = 0
  data_long$time_period[data_long$month=="03"|data_long$month=="04"] = 1
  data_long$time_period[data_long$month=="05"|data_long$month=="06"] = 2
  data_long$time_period[data_long$month=="07"|data_long$month=="08"|data_long$month=="09"] = 3
  
  # trim outliers in each period
  seven <- quantile(data_long$movement, 0.75, na.rm = TRUE)
  twenty <- quantile(data_long$movement, 0.25) 
  iqr <- seven - twenty
  # data_long <- data_long[!(((movement >= (seven  + 6*iqr)) |
  #                                                    (movement <  (twenty- 6 * iqr))))]

  
  plot_data <- data_long %>% 
    group_by(period_num, time_period, conf) %>%
    summarise(x_mean = mean(movement), x_sd = sd(movement))
  
  plot_data$group <- 1
  group <- "group"
  y_var <- ", distance sum (m)"
  facet_labels <- c("March/April", "May/June", "July/August/September")
  labels = NULL
  legend <- "Distance from front"
  
  p <- ggplot(plot_data, aes(x=period_num, y=x_mean,  colour = factor(group))) + 
    geom_point(alpha = 0.5) +
    theme_minimal() +theme(text=element_text(size=21)) +
    ylab(paste0("Movement ", y_var)) + 
    xlab("Minutes after alert") + 
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(color = legend) + 
    stat_smooth(aes(fill=factor(group)), data = plot_data[plot_data$period_num<0,], method = "loess", formula = 'y ~ x')  +
    stat_smooth(aes(fill=factor(group)), data = plot_data[plot_data$period_num>=0,], method = "loess", formula = 'y ~ x') +
    scale_color_manual(values=c("#009392"))+
    scale_fill_manual(values=c("#009392"), guide="none")+
    guides(color=guide_legend(override.aes=list(fill=c("#009392")))) + 
    theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) 
  
  
  
  #==================================================
  # split by distance from frontline
  #==================================================
  
  library(lubridate)
  
  data_long <- read.csv(file.path(data_dir_inter, paste0("coeffs_", "distance_sum_", ".csv"))) %>% setDT()
  
  data_long$month  <- substr(data_long$id, 6, 7)
  
  
  cols <- colnames(data_long)
  cols <- cols[grepl("TYPE", cols)]
  basecols <- cols[!grepl("distance", cols)]
  
  sd <- 276.5610572377576 # this is the standard deviation of distance_to_front and comes from the frontline file
  mn <- 269.8046590875043 # I have pyspark set up but not sparklyr hence why I read it there
  
  newcols <- c()
  for (col in basecols){
    distcol <- paste0("X.", col, "...distance_to_front.")
    upcol <- paste0(col, "upper")
    lowcol <- paste0(col, "lower")
    data_long[, (upcol) := get(col) + get(distcol) * (sd)]
    data_long[, (lowcol) := get(col)]
    
    newcols <- c(newcols, upcol, lowcol)
  }
  
  data_long <- gather(data_long[,- ..cols], period, movement,all_of(newcols),  factor_key=TRUE) %>% setDT()
  

  
  data_long[,period := as.character(period)]
  data_long$conf <- substr(data_long$period, nchar(data_long$period)-4, nchar(data_long$period))
  data_long$period <- gsub('e_TYPE_|upper|lower', '', data_long$period)
  data_long$period <- gsub('[.]', '-', data_long$period)
  
  data_long$period_num <- as.numeric(data_long$period)
  data_long$period_num <- (data_long$period_num)/60
  
  data_long$time_period = 0
  data_long$time_period[data_long$month=="03"|data_long$month=="04"] = 1
  data_long$time_period[data_long$month=="05"|data_long$month=="06"] = 2
  data_long$time_period[data_long$month=="07"|data_long$month=="08"|data_long$month=="09"] = 3
  
  # trim outliers in each period
  for (prd in unique(data_long$time_period)) {
    seven <- quantile(data_long[time_period == prd]$movement, 0.75, na.rm = TRUE)
    twenty <- quantile(data_long[time_period == prd]$movement, 0.25)
    iqr <- seven - twenty
    data_long <- data_long[!((time_period == prd) & ((movement >= (seven  + 10*iqr)) |
                                                          (movement <  (twenty- 10 * iqr))))]
  }
  
  plot_data <- data_long %>% 
    group_by(period_num, time_period, conf) %>%
    summarise(x_mean = mean(movement), x_sd = sd(movement))
  
  group <- "conf"
  y_var <- ", distance sum (m)"
  facet_labels <- c("March/April", "May/June", "July/August/September")
  labels = c("0 km", "275 km (1 sd)")
  legend <- "Distance from front"
  
  plot_data$group <- plot_data[[group]]
  plot_data <- plot_data[!is.na(plot_data$group),]
  
  p <- ggplot(plot_data, aes(x=period_num, y=x_mean, group=group, colour = factor(group, labels = labels))) + 
    geom_point(alpha = 0.5) +
    theme_minimal() +theme(text=element_text(size=21)) +
    ylab(paste0("Movement ", y_var)) + 
    xlab("Minutes after alert") + 
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(color = legend) + 
    stat_smooth(aes(fill=factor(group)), data = plot_data[plot_data$period_num<0,], method = "loess", formula = 'y ~ x')  +
    stat_smooth(aes(fill=factor(group)), data = plot_data[plot_data$period_num>=0,], method = "loess", formula = 'y ~ x') +
    scale_color_manual(values=c("#009392", "#d0587e"))+
    scale_fill_manual(values=c("#009392", "#d0587e"), guide="none")+
    guides(color=guide_legend(override.aes=list(fill=c("#009392",  "#d0587e")))) + 
    theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
    facet_wrap(~ factor(time_period, labels = facet_labels))
  
  
  ggsave(file.path(out_dir_inter, 'distance_sum_inter.pdf'), width = 15, height = 5)
  
  
  
  
  
  #==================================================
  # app downloads
  #==================================================
  
  app <- readxl::read_excel(file.path(final_dir, 'similarweb_downloads.xlsx'), sheet = "Повітряна_тривога") %>% setDT()
  app[, `Unique Installs` := as.numeric(`Unique Installs`) / 1000000]
  app <- app[!is.na(`Unique Installs`)]
  app[, Date := as.Date(Date)]
  app <- app[Date < "2022-10-01"]
  
  ggplot(app, aes(x = Date, y = `Unique Installs`)) + geom_line(color ="#009392", size = 2) + 
    ylim(0,3.5) + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                                           panel.grid.minor.y = element_blank()) + 
    ylab("Monthly unique installs, millions") + xlab("") + scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")
  
  ggsave(file.path("output", "figs", 'app_downloads.pdf'), width = 10, height = 5)
  
   
  
  
  
  
  
  
  
  
  
  
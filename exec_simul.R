hest_by_cutoff <- function(time, status, types, counts) {
	time_levels <- sort(unique(time))
	num_time_levels <- length(time_levels)
	types_levels <- unique(types)
	num_types <- length(unique(types))
	sur <- list()
	haz <- list()
	for(i in 1:num_types) {
		sur[[types_levels[i]]] <- matrix(NA, num_time_levels, num_time_levels)
		haz[[types_levels[i]]] <- matrix(NA, num_time_levels, num_time_levels)
	}
	
	for(i in 1:num_time_levels) {
		print(paste("i = ", i))
		selection <- (time <= time_levels[i])
		survs <- Surv(rep(time[selection], times = counts[selection]),
				rep(status[selection], times = counts[selection]))
		selected_types <- rep(types[selection], times = counts[selection])
		sf <- survfit(survs ~ selected_types)
		local_strata_nm <- sub("selected_types=", "", names(sf$strata), perl = T)
		local_strata <- rep(local_strata_nm, times = sf$strata)
		for(j in 1:length(sf$surv)) {
			print(paste("j = ", j))
			sur[[local_strata[j]]][sf$time[j], i] <- sf$surv[j]
			haz[[local_strata[j]]][sf$time[j], i] <- sf$n.event[j] / sf$n.risk[j]
		}
	}
	return(list(sur=sur, haz=haz))
}

ests <- with(inf_table_all, hest_by_cutoff(time, infected, vaccine, counts))
surv_ests <- ests$sur
haz_ests <- ests$haz

#How estimated survival changes over time
plotCols <- rainbow(34)
plot(1:34, surv_ests$Unvaccinated[1, ], type = "l", xlab = "Cutoff", ylab = "Estimated survival", xlim = c(0, 35), ylim = c(0, 1), main = "Estimated Survival vs Time (Unvaccinated)",
	col = plotCols[1])


for(i in 2:34) {
	points(1:34, surv_ests$Unvaccinated[i, ], type = "l", col = plotCols[i])
}

#How estimated hazard rates changes over time
plotCols <- rainbow(34)
plot(1:34, haz_ests$Unvaccinated[1, ], type = "l", xlab = "Cutoff", ylab = "Estimated hazard", xlim = c(0, 35), ylim = c(0, 0.52), main = "Estimated Hazard vs Time (Unvaccinated)",
	col = plotCols[1])


for(i in 2:34) {
	points(1:34, haz_ests$Unvaccinated[i, ], type = "l", col = plotCols[i])
}

#How estimated hazard ratio changes over time
par(mfrow = c(1, 3))
hr_ests <- list()
for(v in c("Janssen", "Moderna", "Pfizer")){
	hr_ests[[v]] <- haz_ests[[v]] / haz_ests$Unvaccinated
}

#JJ
plotCols <- rainbow(34)
plot(1:34, hr_ests$Janssen[1, ], type = "l", xlab = "Cutoff", ylab = "Estimated hazard ratio", xlim = c(0, 35), ylim = c(0, 1), main = "Estimated Hazard Ratio vs Time (Unvaccinated)",
	col = plotCols[1])

for(i in 2:34) {
	points(1:34, hr_ests$Janssen[i, ], type = "l", col = plotCols[i])
}

#Moderna
plotCols <- rainbow(34)
plot(1:34, hr_ests$Moderna[1, ], type = "l", xlab = "Cutoff", ylab = "Estimated hazard ratio", xlim = c(0, 35), ylim = c(0, 1), main = "Estimated Hazard Ratio vs Time (Unvaccinated)",
	col = plotCols[1])

for(i in 2:34) {
	points(1:34, hr_ests$Moderna[i, ], type = "l", col = plotCols[i])
}

#Pfizer
plotCols <- rainbow(34)
plot(1:34, hr_ests$Pfizer[1, ], type = "l", xlab = "Cutoff", ylab = "Estimated hazard ratio", xlim = c(0, 35), ylim = c(0, 1), main = "Estimated Hazard Ratio vs Time (Unvaccinated)",
	col = plotCols[1])


for(i in 2:34) {
	points(1:34, hr_ests$Pfizer[i, ], type = "l", col = plotCols[i])
}


	
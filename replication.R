library(survival)
library(KMsurv)
library(GGally)

Surv_wrapper <- function(time, status, counts) {
	return(Surv(rep(time, times = counts), rep(status, times = counts)))
}

##To simulate censoring half of the censored individuals before infections
Surv_wrapper2 <- function(time, status, counts) {
	n <- length(time)
	if (length(status) != n || length(counts) != n) {
		print("arguments for Surv_wrapper need to be of identical length")
		return(FALSE)
	}
	
	adjust_factor <- 0.1
	adjust_mat <- matrix(0, sum(status == 0) * 2 + sum(status == 1), 2)
	row_num <- 1
	for (i in 1:n) {
		count <- counts[i]
		if (status[i] == 0) {
			if (count %% 2 == 1) {
				rolloff <- round(runif(1))
				count1 <- floor(count / 2) + rolloff
				count2 <- floor(count / 2) + (1 - rolloff)
			} else {
				count1 <- floor(count / 2)
				count2 <- count1
			}
			adjust_mat[row_num, ] <- c(-1, count1)
			row_num <- row_num + 1
			adjust_mat[row_num, ] <- c(0, count2)
		} else {
			adjust_mat[row_num, ] <- c(0, count)
		}
		row_num <- row_num + 1
	}
	t <- rep(time, times = counts)
	s <- rep(status, times = counts)
	t <- t + rep(adjust_mat[, 1] * adjust_factor, times = adjust_mat[, 2])
	# return(list(t = t, s = s, adjust_mat = adjust_mat))
	return(Surv(t, s))
}

compute_haz <- function(sf) {
	vac_names <- names(sf$strata)
	vac_names <- gsub("vactypes=(\\w+)", "\\1", vac_names, perl = T)
	vacInds <- cumsum(sf$strata)
	vacInds_start <- c(1, vacInds[-length(vacInds)] + 1)
	res <- list()
	for (i in 1:length(vac_names)) {
		start <- vacInds_start[i]
		finish <- vacInds[i]
		numObs <- finish - start + 1
		vname <- vac_names[i]
		print(vname)
		tDiffs <- sf$time[start:finish] - c(0, sf$time[start:(finish - 1)])
				
		hazs <- matrix(0, 2, numObs)
		hazs[1, ] <- sf$time[start:finish]
		hazs[2, ] <- (sf$n.event[start:finish]/sf$n.risk[start:finish])/tDiffs
		res[[vname]] <- hazs
	}
	return(res)
}

compute_hr <- function(haz) {
	unvac_times <- haz$Unvaccinated[1, ]
	unvac_hazs <- haz$Unvaccinated[2, ]
	len <- length(unvac_times)
	res <- list(t = unvac_times)
	for (treat in names(haz)) {
		if (treat != "Unvaccinated") {
			hr <- rep(NA, len)
			treat_haz <- haz[[treat]]
			for (i in 1:length(unvac_times)) {
				unvac_time <- unvac_times[i]
				td_min <- min(abs(treat_haz[1, ] - unvac_time))
				treat_ind <- which(abs(treat_haz[1, ] - unvac_time) == td_min)
				hr[i] <- treat_haz[2, treat_ind]/unvac_hazs[i]
			}
			res[[treat]] <- hr
		}
	}
	return(data.frame(res))
}

		
		
	

par(mfrow = c(1,1))

inf_table_all_unvac <- inf_table_all[inf_table_all$vaccine == "Unvaccinated",]
inf_surv_all_unvac <- with(inf_table_all_unvac, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_unvac)

#lifetab from KMsurv provides an alternative way to do actuarial counts
inftab_all_unvac <- lifetab(0:34, sum(inf_table_all_unvac$counts),
					inf_table_all_unvac$counts[!inf_table_all_unvac$infected],
					inf_table_all_unvac$counts[inf_table_all_unvac$infected])
inf_sf_all_unvac <- survfit(inf_surv_all_unvac ~ 0)

inf_table_all_jj <- inf_table_all[inf_table_all$vaccine == "Janssen",]
inf_surv_all_jj <- with(inf_table_all_jj, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_jj)

inf_table_all_moderna <- inf_table_all[inf_table_all$vaccine == "Moderna",]
inf_surv_all_moderna <- with(inf_table_all_moderna, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_moderna)

inf_table_all_pb <- inf_table_all[inf_table_all$vaccine == "Pfizer",]
inf_surv_all_pb <- with(inf_table_all_pb, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_pb)

##altogether
inf_surv_all <- with(inf_table_all, Surv_wrapper(time, infected, counts))
vactypes <- rep(inf_table_all$vaccine, times = inf_table_all$counts)
inf_sf_all <- survfit(inf_surv_all ~ vactypes)
ggsurv(inf_sf_all)

inf_haz_all <- compute_haz(inf_sf_all)
inf_hr_all <- compute_hr(inf_haz_all)

plot(inf_hr_all$t, inf_hr_all$Janssen, col = "orange", type = "l")
points(inf_hr_all$t, inf_hr_all$Moderna, col = "blue", type = "l")
points(inf_hr_all$t, inf_hr_all$Pfizer, col = "green", type = "l")



##simulated counts with Surv_wrapper2
##corresponds roughly to lifetable method in SAS's PROC LIFETEST
inf_surv_all2 <- with(inf_table_all, Surv_wrapper2(time, infected, counts))
inf_sf_all2 <- survfit(inf_surv_all2 ~ vactypes)
ggsurv(inf_sf_all2)


##age <50
inf_surv_under50 <- with(inf_table_under50, Surv_wrapper(time, infected, counts))
vactypes <- with(inf_table_under50, rep(vaccine, times = counts))
inf_sf_under50 <- survfit(inf_surv_under50 ~ vactypes, conf.type = "none")
ggsurv(inf_sf_under50, CI = FALSE)

inf_haz_under50 <- compute_haz(inf_sf_under50)
inf_hr_under50 <- compute_hr(inf_haz_under50)

plot(inf_hr_under50$t, inf_hr_under50$Janssen, col = "orange", type = "l")
points(inf_hr_under50$t, inf_hr_under50$Moderna, col = "blue", type = "l")
points(inf_hr_under50$t, inf_hr_under50$Pfizer, col = "green", type = "l")


##age >= 65
inf_surv_50_64 <- with(inf_table_50_64, Surv_wrapper(time, infected, counts))
vactypes <- with(inf_table_50_64, rep(vaccine, times = counts))
inf_sf_50_64 <- survfit(inf_surv_50_64 ~ vactypes, conf.type = "none")
ggsurv(inf_sf_50_64, CI = FALSE)

inf_haz_50_64 <- compute_haz(inf_sf_50_64)
inf_hr_50_64 <- compute_hr(inf_haz_50_64)

plot(inf_hr_50_64$t, inf_hr_50_64$Janssen, col = "orange")
points(inf_hr_50_64$t, inf_hr_50_64$Moderna, col = "blue")
points(inf_hr_50_64$t, inf_hr_50_64$Pfizer, col = "green")

##age 65+
inf_surv_over65 <- with(inf_table_over65, Surv_wrapper(time, infected, counts))
vactypes <- with(inf_table_over65, rep(vaccine, times = counts))
inf_sf_over65 <- survfit(inf_surv_over65 ~ vactypes, conf.type = "none")
ggsurv(inf_sf_over65, CI = FALSE)

inf_haz_over65 <- compute_haz(inf_sf_over65)
inf_hr_over65 <- compute_hr(inf_haz_over65)

plot(inf_hr_over65$t, inf_hr_over65$Janssen, col = "orange", type = "l", main = "Over 65 Group", ylab = "Hazard", xlab = "Week")
points(inf_hr_over65$t, inf_hr_over65$Moderna, col = "blue", type = "l")
points(inf_hr_over65$t, inf_hr_over65$Pfizer, col = "green", type = "l")


##cox ph
inf_table_combined <- rbind(inf_table_under50, inf_table_50_64, inf_table_over65)
inf_table_combined$age_group <- as.factor(rep(c("under50", "50_64", "over65"),
						times = c(dim(inf_table_under50)[1], dim(inf_table_50_64)[1], dim(inf_table_over65)[1])))
inf_surv_combined <- with(inf_table_combined, Surv_wrapper(time, infected, counts))
vactypes <- with(inf_table_combined, rep(vaccine, times = counts))


for (vac in unique(inf_table_combined$vaccine)) {
	inf_table_combined[[paste("is_", vac, sep = "")]] <- inf_table_combined$vaccine == vac
}
inf_survcovars_combined <- list()
for(field in names(inf_table_combined)) {
	if (field != "counts") {
		inf_survcovars_combined[[field]] <- rep(inf_table_combined[[field]], times = inf_table_combined$counts)
	}
}
inf_survcovars_combined <- data.frame(inf_survcovars_combined)

dim(inf_surv_combined)
dim(inf_survcovars_combined)
		



##test against known result
test_sf <- survfit(inf_surv_combined ~ vactypes)
test_sf$cumhaz - inf_sf_all$cumhaz #difference expected because of inconsistency in total
age_group <- with(inf_table_combined, rep(age_group, times = counts))
inf_sf_combined <- survfit(inf_surv_combined ~ vactypes + age_group)
##

inf_sf_all

inf_combined <- make_dummy(inf_table_combined)
inf_combined <- inf_combined[(inf_combined$counts != 0),]
rownames(inf_combined) <- 1:nrow(inf_combined)

inf_surv_combined <- Surv(inf_combined$time, inf_combined$infected)

inf_cox0 <- coxph(Surv(time, infected) ~ age_group + is_Janssen + is_Moderna + is_Pfizer, weights = counts, data = inf_combined, iter.max = 1000)
inf_cox

#is the Surv object reweighed using weights provided?
S <- Surv(inf_combined$time, inf_combined$infected)
inf_cox <- coxph(S ~ age_group + is_Janssen + is_Moderna + is_Pfizer, weights = counts, data = inf_combined, iter.max = 1000)
inf_cox
#yes

inf_cox <- coxph(inf_surv_combined ~ age_group + is_Janssen + is_Moderna + is_Pfizer, weights = counts, data = inf_combined, iter.max = 1000)
inf_cox

inf_cox <- coxph(inf_surv_combined ~ age_group + is_Janssen + is_Moderna + is_Pfizer + time + time_Janssen + time_Moderna + time_Pfizer, weight = counts, data = inf_combined, ties="efron", x = T, y = T, iter.max = 1000)
inf_cox


##Remark: convergence issue with time covariate, which is exceptionally large in magnitude (19)
inf_noj <- inf_combined[(inf_combined$vaccine != "Janssen"),]
rownames(inf_noj) <- 1:nrow(inf_noj)
inf_surv_noj <- Surv(inf_noj$time, inf_noj$infected)
inf_cox_noj <- coxph(inf_surv_noj ~ age_group + is_Moderna + is_Pfizer + time + time_Moderna + time_Pfizer, weights = counts, data = inf_noj, iter.max = 1000)
coxph(inf_surv_noj ~ age_group + is_Moderna + is_Pfizer + time_Moderna + time_Pfizer, weights = counts, data = inf_noj, iter.max = 1000)
#Janssen is not the cause

bh <- basehaz(inf_cox)
plot(inf_combined$time, inf_cox$residuals)
plot(bh$time, bh$hazard)
plot(bh$time, log(bh$hazard))

logbh <- log(bh$hazard)
lm_bh <- lm(logbh ~ bh$time)
lm_bh
anova(lm_bh)
##interpretation: baseline hazard is exponential and mostly cancels out effect of the time covar

##try using tt instead
trans <- function (x, t, ...) {
	return(x*t)
}

inf_cox_pp <- coxph(inf_surv_combined ~ age_group + is_Janssen + is_Moderna + is_Pfizer + tt(is_Janssen) + tt(is_Moderna) + tt(is_Pfizer), data = inf_combined, tt = trans, weight = counts, ties = "efron")
inf_cox_pp2 <- coxph(inf_surv_combined ~ age_group + is_Janssen + is_Moderna + is_Pfizer + tt(is_Janssen) + tt(is_Moderna) + tt(is_Pfizer), data = inf_combined, tt = trans, weight = counts, ties = "efron")
##this works





inf_zph2 <- cox.zph(inf_cox, transform = function(time){ log(time + 20) }, terms = F)
inf_cox2 <- coxph(inf_surv_combined ~ tt(vaccine) + age_group,
			data = inf_survcovars_combined,
			tt = function(x, t, ...){
				vac <- (x == c("Janssen", "Moderna", "Pfizer")) * t
				cbind(Janssen = (x == "Janssen"), Moderna = (x == "Moderna"), Pfizer = (x == "Pfizer"))
			})
save(inf_cox2, file="inf_cox2", compress = F)

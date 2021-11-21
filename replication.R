library(survival)
library(KMsurv)

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

##simulated counts with Surv_wrapper2
inf_surv_all2 <- with(inf_table_all, Surv_wrapper2(time, infected, counts))
inf_sf_all2 <- survfit(inf_surv_all2 ~ vactypes)
ggsurv(inf_sf_all2)



inf_surv_under50 <- with(inf_table_under50, Surv_wrapper(time, infected, counts))
vactypes <- with(inf_table_under50, rep(vaccine, times = counts))
inf_sf_under50 <- survfit(inf_surv_under50 ~ vactypes, conf.type = "none")
ggsurv(inf_sf_under50, CI = FALSE)

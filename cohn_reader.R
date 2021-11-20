library(readxl)

cohn_all <- read_xlsx("table_S3.xlsx", sheet = "A")
cohn_under50 <- read_xlsx("table_S3.xlsx", sheet = "B")
cohn_50_64 <- read_xlsx("table_S3.xlsx", sheet = "C")
cohn_over65 <- read_xlsx("table_S3.xlsx", sheet = "D")

mort_under65 <- read_xlsx("table_S3.xlsx", sheet = "E")
mort_over65 <- read_xlsx("table_S3.xlsx", sheet = "F")
mort_low_comorb <- read_xlsx("table_S3.xlsx", sheet = "G")
mort_high_comorb <- read_xlsx("table_S3.xlsx", sheet = "H")

summary_infs <- function(df_inf) {
	num_rows <- dim(df_inf)[1]
	odds_unvac <- df_inf[[2]] / df_inf[[3]]
	odds_jj <- df_inf[[4]] / df_inf[[5]]
	odds_moderna <- df_inf[[6]] / df_inf[[7]]
	odds_pb <- df_inf[[8]] / df_inf[[9]]
	
	weeks <- 1:num_rows
	
	plot(weeks, odds_unvac, pch = 1, type = "b", col = "red", lty = 1, xlab = "Week", ylab = "Odds of infection",
		ylim = c(0, 1.2))

	points(weeks, odds_jj, pch = 2, col = "orange", type = "b", lty = 1)
	points(weeks, odds_moderna, pch = 3, col = "green", type = "b", lty = 1)
	points(weeks, odds_pb, pch = 4, col = "blue", type = "b", lty = 1)
	
	or_jj <- odds_jj / odds_unvac
	or_moderna <- odds_moderna / odds_unvac
	or_pb <- odds_pb / odds_unvac
	
	plot(weeks, or_jj, pch = 2, col = "orange", type = "b", lty = 1, xlab = "Week", ylab = "Odds ratio",
		ylim = c(0, 1.5))
	points(weeks, or_moderna, pch = 3, col = "green", type = "b", lty = 1)
	points(weeks, or_pb, pch = 4, col = "blue", type = "b", lty = 1)
	
	odds <- list(unvac = odds_unvac, jj = odds_jj, moderna = odds_moderna, pb = odds_pb)
	or <- list(jj = or_jj, moderna = or_moderna, pb = or_pb)
	res <- list(odds = odds, or = or)
	return(res)
}

make_datatable_inf <- function(df) {
	num_rows <- dim(df)[1]
	num_cols <- dim(df)[2]
	names_vec <- names(df)[2:num_cols]
	time <- rep(df[[1]], num_cols - 1)
	vaccine <-  rep(gsub("(\\w+) (?:\\+|\\-)", "\\1", names_vec, perl = T), each = num_rows) #vector of vaccinate names
	infected <- rep(grepl("\\+", names_vec), each = num_rows) #T for +, F for -
	counts <- rep(0, num_rows * (num_cols - 1))
	for (i in 2:num_cols) {
		counts[((i - 2) * num_rows + 1):(num_rows * (i - 1))] <- df[[i]]
	}

	res <- data.frame(time, vaccine, infected, counts)
	return(res)
}

make_datatable_mort <- function(df) {
	num_rows <- dim(df)[1]
	num_cols <- dim(df)[2]
	names_vec <- names(df)[2:num_cols]
	time <- rep(df[[1]], num_cols - 1)
	vaccinated <- rep(grepl("Vac", names_vec), each = num_rows)
	infected <- rep(grepl("PCR\\+", names_vec), each = num_rows)
	dead <- rep(!grepl("Non-deaths", names_vec), each = num_rows)
	counts <- rep(0, num_rows * (num_cols - 1))
	for (i in 2:num_cols) {
		counts[((i - 2) * num_rows + 1):(num_rows * (i - 1))] <- df[[i]]
	}

	res <- data.frame(time, vaccinated, infected, dead, counts)
	return(res)
}
	

par(mfcol = c(2, 4))

summ_all <- summary_infs(cohn_all)
summ_under50 <- summary_infs(cohn_under50)
summ_50_64 <- summary_infs(cohn_50_64)
summ_over65 <- summary_infs(cohn_over65)

inf_table_all <- make_datatable_inf(cohn_all)
inf_table_under50 <- make_datatable_inf(cohn_under50)
inf_table_50_64 <- make_datatable_inf(cohn_50_64)
inf_table_over65 <- make_datatable_inf(cohn_over65)

mort_table_under65 <- make_datatable_mort(mort_under65)
mort_table_over65 <- make_datatable_mort(mort_over65)
mort_table_low_comorb <- make_datatable_mort(mort_low_comorb)
mort_table_high_comorb <- make_datatable_mort(mort_high_comorb)


## Data QA check

# S3A should be the sum of S3B-D. Verify
sum_matches <- rep(NA, dim(cohn_all)[2] - 1)
sum_computed <- Reduce(f = "+", x = list(cohn_under50, cohn_50_64, cohn_over65))
for (i in 2:dim(cohn_all)[2]) {
	sum_matches[i - 1] <- prod(cohn_all[[i]] == sum_computed[[i]])
}

which(!sum_matches)

which(sum_computed[[6 + 1]] != cohn_all[[6 + 1]])
sum_computed[[7]][7] - cohn_all[[7]][7]
# mismatched cell is Moderna/PCR-/Week 7
# mismatch confirmed to exist in source tables (S3A-D)


# S3A column sums should be equal or greater (in case of multiple assays/patient)
# compared to table S1 vaccination status counts

inf_freq_mat <- matrix(NA, 4, 2)
vacs <- unique(inf_table_all$vaccine)
inf_status <- c(F, T) # PCR- first for compatibility with Table S1
for (i in 1:length(vacs)) {
	for (j in 1:length(inf_status)) {
		inf_freq_mat[i, j] <- sum(inf_table_all$counts[(inf_table_all$vaccine == vacs[i]) & (inf_table_all$infected == inf_status[j])])
	}
}

S1_vac <- read_xlsx("table_S1.xlsx", sheet = "Vaccine")
S1_vac_count <- as.matrix(S1_vac[,2:3])
sum(S1_vac_count) #780225, matches the reported total

vac_count_deviation <- S1_vac_count - inf_freq_mat
vac_count_reldev <- vac_count_deviation / inf_freq_mat #relative deviation
vac_count_deviation
vac_count_reldev
sum(vac_count_deviation)
## Table S1 counts are smaller than computed counts for unvaccinated, larger for the vaccinated categories
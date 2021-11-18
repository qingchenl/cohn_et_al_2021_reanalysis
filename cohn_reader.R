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

par(mfcol = c(2, 4))

summ_all <- summary_infs(cohn_all)
summ_under50 <- summary_infs(cohn_under50)
summ_50_64 <- summary_infs(cohn_50_64)
summ_over65 <- summary_infs(cohn_over65)
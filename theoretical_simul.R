# Demonstrate how unrepresentative sample can affect apparent hazard ratio
# Lifetime randomly generated using exponential distribution, so Markovian

N <- 10^5

set.seed(1234)

time_steps <- 9
x <- 1:time_steps

times_a <- rexp(N, 0.02)
times_b <- rexp(N, 0.01)

bin_counts <- function(times, cutoffs) {
	N <- length(times)
	cutoffs <- sort(cutoffs) # make sure the cutoffs are ordered
	numBins <- length(cutoffs) + 1
	bins <- rep(0, numBins)
	bins[1] <- sum(times <= cutoffs[1])
	print(bins)
	for (i in 2:(numBins - 1)) {
		bins[i] = sum((times > cutoffs[i - 1]) & (times <= cutoffs[i]))
	}
	bins[numBins] <- sum(times > cutoffs[length(cutoffs)])
	bin_names <- c(cutoffs, paste(">", as.character(cutoffs[length(cutoffs)]), sep = ""))
	#print("naming")
	#print(cat(">", as.character(cutoffs[length(cutoffs)])))
	#print(bin_names)
	names(bins) <- bin_names
	return(bins)
}

deaths_a <- bin_counts(times_a, x)
deaths_a

deaths_b <- bin_counts(times_b, x)
deaths_b

haz_calc <- function(N, deaths, losses) {
	num_bins <- length(deaths)
	if (missing(losses)) {
		losses <- rep(0, num_bins)
	}
	N_eff <- N
	haz <- rep(0, num_bins)
	for (i in 1:num_bins) {
		N_eff <- N_eff - losses[i]
		haz[i] <- deaths[i] / N_eff
		N_eff <- N_eff - deaths[i]
	}
	return(haz)
}
	
haz_a <- haz_calc(N, deaths_a)
haz_a

haz_b <- haz_calc(N, deaths_b)
haz_b

haz_ratio <- haz_b / haz_a
haz_ratio

# Fix positivity rate
tested_a <- ceiling(sum(deaths_a[1:(length(deaths_a) - 1)]) * 3)
tested_b <- ceiling(sum(deaths_b[1:(length(deaths_b) - 1)]) * 3)

haz_a_app <- haz_calc(tested_a, deaths_a)
haz_b_app <- haz_calc(tested_b, deaths_b)
haz_ratio_app <- haz_b_app / haz_a_app

haz_a_app
haz_b_app
haz_ratio_app

plot(1:9, haz_ratio_app[1:9])

# Uneven positivity rate
tested_a <- ceiling(sum(deaths_a[1:(length(deaths_a) - 1)]) * 10)
tested_b <- ceiling(sum(deaths_b[1:(length(deaths_b) - 1)]) * 3)

haz_a_app <- haz_calc(tested_a, deaths_a)
haz_b_app <- haz_calc(tested_b, deaths_b)
haz_ratio_app <- haz_b_app / haz_a_app

haz_a_app
haz_b_app
haz_ratio_app

plot(1:9, haz_ratio_app[1:9])
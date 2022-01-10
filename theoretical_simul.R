# Demonstrate how unrepresentative sample can affect apparent hazard ratio
# Lifetime randomly generated using exponential distribution, so Markovian

library(survival)

N <- 10^5

set.seed(1234)

time_steps <- 9
x <- 1:time_steps

times_a <- rexp(N, 0.02)
times_b <- rexp(N, 0.01)

# utility function to produce all possible n-ary Cartesian product terms
list_combs <- function(levels) {
	combs <- matrix(0, prod(levels), length(levels))
	#level_mat <- matrix(levels, prod(levels), levels, byrow = T)
	inds <- 0:(prod(levels) - 1)
	
	for (i in length(levels):1) {
		if (i == length(levels)) {
			combs[, i] <- inds %% levels[i]
		} else {
			combs[, i] <- floor(inds / prod(levels[length(levels):(i+1)])) %% levels[i]
		}
	}
	return(combs + 1)
}

list_combs(c(2,3,5))

#returns number of observations between cutoffs
bin_counts <- function(times, cutoffs) {
	N <- length(times)
	cutoffs <- sort(cutoffs) # make sure the cutoffs are ordered
	numBins <- length(cutoffs) + 1
	bins <- rep(0, numBins)
	bins[1] <- sum(times <= cutoffs[1])
	# print(bins)
	for (i in 2:(numBins - 1)) {
		bins[i] = sum((times > cutoffs[i - 1]) & (times <= cutoffs[i]))
	}
	bins[numBins] <- sum(times > cutoffs[length(cutoffs)])
	bin_names <- c(cutoffs, paste(">", as.character(cutoffs[length(cutoffs)]), sep = ""))
	names(bins) <- bin_names
	return(bins)
}

deaths_a <- bin_counts(times_a, x)
deaths_a

deaths_b <- bin_counts(times_b, x)
deaths_b

wrexp_rp <- function(n, rate, cutoff, prop) { #weighted sampling from exponential distribution, random proportion
	if (missing(n)) {
		n <- 1
	}
	if (prop <= 0) {
		print("wrexp_rp: prop argument needs to be positive")
		return(NULL)
	}
	p <- runif(n)
	p_cutoff <- pexp(cutoff, rate)
	p_exp <- ((p <= prop) * (p / prop) * p_cutoff) + (p > prop) * (p_cutoff + ((p - prop) / (1 - prop)) * (1 - p_cutoff))
	return(qexp(p_exp, rate))
}

wrexp_cs <- function(n, rate, cutoff1, cutoff2) { #sampling from exponential distribution, conditional on section
	if (missing(n)) {
		n <- 1
	}
	tol <- 10^-10
	if (cutoff2 - cutoff1 == 0) {
		print("wrexp_cs: cutoffs can't be identical")
		return(NULL)
	} else if (abs(cutoff2 - cutoff1) < tol) {
		warning("cutoffs closer than tolerance")
	}
	if (cutoff2 < cutoff1) {
		lims <- c(cutoff2, cutoff1)
	} else {
		lims <- c(cutoff1, cutoff2)
	}
	
	p <- runif(n)
	plims <- pexp(lims, rate)
	p_exp <- plims[1] + p * (plims[2] - plims[1])
	return(qexp(p_exp, rate))
}
	
	

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

times_to_obs <- function(lifetimes, cutoff) {
	time <- lifetimes * (lifetimes <= cutoff) + cutoff * (lifetimes > cutoff) #use lifetime if not greater than cutoff, cutoff time (when they are censored) otherwise
	infected <- lifetimes <= cutoff #status TRUE (1) for event if it happens before cutoff, status FALSE (0) if it happens after cutoff
	res <- data.frame(time, infected)
	attr(res, "coxable") <- TRUE
	attr(res, "cutoff") <- cutoff
	return(res)
}

draw_sample <- function(coxable, times, samp_sizes) {
	num_samples <- length(times)
	total_records <- nrow(coxable)
	
	#check for obvious issues in input
	if (num_samples != length(samp_sizes)) {
		print("draw_sample: times and samp_sizes of different lengths were provided")
		return(NULL)
	}
	if (attr(coxable, "coxable") != TRUE) {
		warning("coxable argument for draw_sample may not be in the right format")
	}
	if (sum(sample_sizes) > total_records) {
		print("draw_sample: coxable doesn't have enough records to draw the specified samples from")
		return(NULL)
	}
	
	#make sure times are sorted; samp_sizes will take same ordering as times
	indx <- order(samp_sizes)
	samp_sizes <- samp_sizes[indx]
	names(samp_sizes) <- 1:num_samples
	times <- times[indx]
	names(times) <- 1:num_samples
	
	if (!is.null(attr(coxable, "cutoff")) && attr(coxable, "cutoff") < times[num_samples]) {
		warning("some of times may exceed original followup time")
	}
	
	infected <- rep(NA, num_samples * 2) #logical vector for infection status
	counts <- rep(NA, num_samples * 2) #count of individuals in this group
	remaining <- 1:num_samples
	time <- rep(times, each = 2)
	
	for (i in 1:num_samples) {
		#check that the sample can be drawn
		if (sample_sizes[i] > length(remaining)) {
			print(paste("draw_sample: sample of size ", sample_sizes[i], " required at time ", i, ", but only ", length(remaining), " records remain", sep = ""))
			return(NULL)
		}
		
		#take sample and count infected/censored
		samp <- sample(remaining, sample_sizes[i], replace = F)
		num_infected <- sum(coxable$time[samp] <= times[i] & coxable$infected[samp])
		num_censored <- sum(coxable$time[samp] > times[i] | !coxable$infected[samp])
		
		#fill respective vector entries with data
		index1 <- i * 2 - 1
		index2 <- index1 + 1
		infected[index1] <- TRUE
		infected[index2] <- FALSE
		counts[index1] <- num_infected
		counts[index2] <- num_censored
		
		allInfected <- which(coxable$time <= times[i])
		#update list of remaining entries
		remaining <- remaining[!(remaining %in% union(samp, allInfected))]
	}
	res <- data.frame(time, infected, counts)
	attr(res, "redrawn") <- TRUE
	attr(res, "coxable") <- TRUE
	attr(res, "counts") <- TRUE
	if (!is.null(attr(coxable, "cutoff"))) {
		attr(res, "cutoff") <- min(c(attr(coxable, "cutoff"), max(times))) #new cutoff is time of last observation or cutoff of the original, whichever is earlier
	}
	return(res)
}
	
haz_a <- haz_calc(N, deaths_a)
haz_a

haz_b <- haz_calc(N, deaths_b)
haz_b

haz_ratio <- haz_b / haz_a
haz_ratio

#using coxph:



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


## Take a set of rates and generate dataset
gen_data <- function(rates, N, groups) {
	if (missing(rates)|missing(N)) {
		print("gen_data: rates and N are required arguments")
		return(NULL)
	} else if (length(N) > 1 & length(N) != length(rates)) {
		print("gen_data: either provide one N (identical sample sizes) or one N for each rate")
		return(NULL)
	} else if (length(groups) > 1 & length(groups) != length(rates)) {
		print("gen_data: either provide one group name or one group name for each rate")
		return(NULL)
	}
	if (missing(groups)) {
		groups = "group"
	}
	if (length(rates) > 1) {
		if (length(N) == 1) {
			N <- rep(N, length(rates))
		}
		if (length(groups) == 1) {
			group_name_template <- groups
			for (i in 1:length(rates)) {
				groups[i] <- paste(group_name_template, i, sep = "")
			}
		}
	}
	
	res <- data.frame(times = numeric(0), group = character(0))
	for (i in 1:length(rates)) {
		observations <- data.frame(time = rexp(N[i], rates[i]), group = rep(groups[i], N[i]))
		res <- rbind(res, observations)
	}
	attr(res, "surv_time") <- TRUE
	return(res)
}

set_status <- function(surv_times, cutoff, status_name = "infected") {
	if(!is.data.frame(surv_times)) {
		data.frame(time = surv_times)
	}

	if (missing(cutoff)) {
		surv_times[[status_name]] <- rep(T, nrow(surv_times))
	} else {
		surv_times[[status_name]] <- surv_times$time <= cutoff
		attr(surv_times, "cutoff") <- cutoff
	}
	attr(surv_times, "coxable") <- TRUE
	return(surv_times)
}

gd <- gen_data(c(0.1, 0.5, 0.01), 10^3, c("group1", "group2", "group3"))
gd <- set_status(gd, status_name = "status")

coxph(Surv(time, status) ~ group, data = gd) #hazard ratios of 5 and 0.1, as expected

gd2 <- set_status(gd, 2, status_name = "status") #censoring data after cutoff
coxph(Surv(time, status) ~ group, data = gd2) #same result, as expected

gd2$is_group2 <- gd2$group == "group2"
gd2$is_group3 <- gd2$group == "group3"
coxph(Surv(time, status) ~ is_group2 + is_group3, data = gd2)

trans <- function (x, t, ...) {
	return(x*t)
}

coxph(Surv(time, status) ~ group + tt(is_group2) + tt(is_group3), data = gd2, tt = trans)


cuts <- seq(0.5, 20, by = 0.2)

check_match <- function(df, conditions) {
	num_rows <- nrow(df)
	#res <- rep(0, num_rows)
	cond_names <- names(conditions)
	checks <- matrix(0, num_rows, length(conditions))
	for (i in 1:length(conditions)) {
		checks[, i] <- df[[cond_names[i]]] == conditions[[i]]
	}
	return(as.logical(apply(checks, 1, prod)))
}

make_counts <- function(rawdata, cuts, factors, status_name = "status", time_name = "time") {
	template <- rawdata[F,] #dataframe with appropriate columns and 0 rows
	template$counts <- numeric(0)
	
	if (missing(factors)) {
		factors <- character(0)
		for (nm in names(rawdata)) {
			if (nm != status_name & nm != time_name) {
				factors <- c(factors, nm) #inefficient, but if there are too many factors we have bigger problems to worry about
			}
		}
	}
	factor_levels <- rep(1, length(factors))
	factor_values <- list()
	for (i in 1:length(factors)) {
		factor_values[[i]] <- unique(rawdata[[factors[i]]])
		#print(factor_values)
		factor_levels[i] <- length(factor_values[[i]])
	}
	
	comb_mat <- list_combs(factor_levels)
	#print(comb_mat)
	curr_facs <- list()
	#selections <- list()
	for (i in 1:nrow(comb_mat)) {
		for (j in 1:length(factors)) {
			curr_facs[[factors[j]]] <- factor_values[[j]][comb_mat[i, j]]
		}
		#return(curr_facs)
		raw_selection <- check_match(rawdata, curr_facs)
		#selections[[i]] <- raw_selection
		#print(head(raw_selection))
		if (sum(raw_selection) > 0) { #at least 1 row that matches criteria
			temp_data <- bin_counts(rawdata[[time_name]][raw_selection], cuts)
			tm_point <- names(temp_data)
			for (tm in tm_point) {
				if(temp_data[[tm]] > 0) {
					rw <- curr_facs
					rw$counts <- temp_data[[tm]]
					if(is.na(suppressWarnings(as.numeric(tm)))) {
						rw[[time_name]] <- max(cuts)
						rw[[status_name]] <- FALSE
					} else {
						rw[[time_name]] <- as.numeric(tm)
						rw[[status_name]] <- TRUE
					}
					#print(template)
					#print(rw)
					template <- rbind(template, rw)
				}
			}
		}
	}
	return(template)
}


gd_discr <- make_counts(gd)
coxph(Surv(time, status) ~ group, data = gd_discr, weight = counts)


gen_data2 <- function(rates, N, groups, distr, prop, cutoff1, cutoff2) {
	if (missing(rates)||missing(N)) {
		print("gen_data: rates and N are required arguments")
		return(NULL)
	} else if (length(N) > 1 && length(N) != length(rates)) {
		print("gen_data: either provide one N (identical sample sizes) or one N for each rate")
		return(NULL)
	} else if (length(groups) > 1 && length(groups) != length(rates)) {
		print("gen_data: either provide one group name or one group name for each rate")
		return(NULL)
	}
	if (missing(groups)) {
		groups = "group"
	}
	if (missing(distr)) {
		distr = "exp"
	}
	
	print(distr)
	print("check")
	if (length(rates) > 1) {
		if (length(N) == 1) {
			N <- rep(N, length(rates))
		}
		if (length(groups) == 1) {
			group_name_template <- groups
			for (i in 1:length(rates)) {
				groups[i] <- paste(group_name_template, i, sep = "")
			}
		}
		print(missing(prop))
		if (!missing(prop) && length(prop) == 1) {
			prop <- rep(prop, length(rates))
		}
		
		if (!missing(cutoff1) && length(cutoff1) == 1) {
			cutoff1 <- rep(cutoff1, length(rates))
		}
		if (!missing(cutoff2)) {
			if (length(cutoff2) == 1) {
				cutoff2 <- rep(cutoff2, length(rates))
			}
			if (sum(cutoff2 <= cutoff1) > 1) {
				print("gen_data2: cutoff2 needs to be greater than cutoff1")
				return(NULL)
			}
		}
	}
	print("check1")
	if (distr == "exp") {
		rdistr <- rexp
	} else if (distr == "wrexp_rp") {
		if (missing(prop) || missing(cutoff1)) {
			print("gen_data2: wrexp_rp distribution requires prop and cutoff1 arguments")
			return(NULL)
		} else {
			rdistr <- wrexp_rp
		}
	} else if (distr == "wrexp_cs") {
		if (missing(cutoff1) || missing(cutoff2)) {
			print("gen_data2: wrexp_cs distribution requires cutoff1 and cutoff2 arguments")
		} else {
			rdistr <- wrexp_cs
		}
	} else {
		print(paste("gen_data2: ", distr, " is not a recognized distribution"))
		return(NULL)
	}
	
	res <- data.frame(times = numeric(0), group = character(0))
	for (i in 1:length(rates)) {
		if (distr == "exp") {
			arg_list <- list(n = N[i], rate = rates[i])
		} else if (distr == "wrexp_rp") {
			arg_list <- list(n = N[i], rate = rates[i], cutoff = cutoff1[i], prop = prop[i])
		} else if (distr == "wrexp_cs") {
			arg_list <- list(n = N[i], rate = rates[i], cutoff1 = cutoff1[i], cutoff2 = cutoff2[i])
		}
		times <- do.call(rdistr, arg_list)
		observations <- data.frame(time = times, group = rep(groups[i], N[i]))
		res <- rbind(res, observations)
	}
	attr(res, "surv_time") <- TRUE
	return(res)
}

#tests
df <- gen_data2(c(0.1, 0.5), 1000, c("group1", "group2"))
coxph(Surv(time, rep(1, nrow(df))) ~ group, data = df)

df <- gen_data2(c(0.1, 0.5), 1000, c("group1", "group2"), "wrexp_rp", c(pexp(0.1, 10), pexp(0.5, 10)), 10) #expect HR of ~5
coxph(Surv(time, rep(1, nrow(df))) ~ group, data = df)

df <- gen_data2(c(0.1, 0.5), 1000, c("group1", "group2"), "wrexp_rp", c(0.1, 0.1), 10)
coxph(Surv(time, rep(1, nrow(df))) ~ group, data = df)

df <- gen_data2(c(0.1, 0.5), 1000, c("group1", "group2"), "wrexp_cs", cutoff1 = 0, cutoff2 = 1000) #expect HR of ~5
coxph(Surv(time, rep(1, nrow(df))) ~ group, data = df)

df <- gen_data2(c(0.1, 0.5), 1000, c("group1", "group2"), "wrexp_cs", cutoff1 = 1, cutoff2 = 5)
coxph(Surv(time, rep(1, nrow(df))) ~ group, data = df)
#





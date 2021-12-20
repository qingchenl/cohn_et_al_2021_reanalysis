# figuring out what coxph is doing


#terms.inner is not exported, manually copying code from github
terms.inner <- function(x) {
    if (inherits(x, "formula")) {
        if (length(x) ==3) c(terms.inner(x[[2]]), terms.inner(x[[3]]))
        else terms.inner(x[[2]])
    }
    else if (inherits(x, "call") && 
             (x[[1]] != as.name("$") && x[[1]] != as.name("["))) {
        if (x[[1]] == '+' || x[[1]]== '*' || x[[1]] == '-' || x[[1]] ==':') {
            # terms in a model equation, unary minus only has one argument
            if (length(x)==3) c(terms.inner(x[[2]]), terms.inner(x[[3]]))
            else terms.inner(x[[2]])
        }
        else if (x[[1]] == as.name("Surv"))
                 unlist(lapply(x[-1], terms.inner))
        else terms.inner(x[[2]])
    }
    else(deparse(x))
}

#generate dummy variable columns for vaccination and interaction terms
make_dummy <- function(df) {
	res <- df
	vac_types <- c("Janssen", "Moderna", "Pfizer")
	for(vac in vac_types) {
		temp_name <- paste("is_", vac, sep = "")
		res[[temp_name]] <- (df$vaccine == vac)
		res[[paste("time_", vac, sep = "")]] <- res[[temp_name]] * res$time
	}
	return(res)
}

#convert counts to individuals
expand_counts <- function(df) {
	fields <- names(df)
	res <- list()
	for (fld in fields) {
		if (fld != "counts") {
			res[[fld]] <- rep(df[[fld]], times = df$counts)
		}
	}
	return(data.frame(res))
}


#probe function
coxph_part <- function(formula, data, weights, subset, na.action,
        init, control, ties= c("efron", "breslow", "exact"),
        singular.ok =TRUE,  robust,
        model=FALSE, x=FALSE, y=TRUE,  tt, method=ties, 
        id, cluster, istate, statedata, nocenter=c(-1, 0, 1), ...) {

	ties <- match.arg(ties)
	Call <- match.call()
	#print(ties)
	extraArgs <- list(...)
	#print(extraArgs)
	## if (length(extraArgs))... is omitted because we aren't using any

	if (missing(control)) control <- coxph.control(...) 
	#print(control)
	ss <- "cluster"
	if (is.list(formula))
		Terms <- if (missing(data)) terms(formula[[1]], specials=ss) else
                 		terms(formula[[1]], specials=ss, data=data)
    	else Terms <- if (missing(data)) terms(formula, specials=ss) else
                 terms(formula, specials=ss, data=data)

	tcl <- attr(Terms, 'specials')$cluster
	#print(tcl)
	#omit everything tcl related
    	indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster", "id", "istate"),
                  names(Call), nomatch=0)
	#print(indx)
	if (indx[1] ==0) stop("A formula argument is required")
    	tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    	#print(tform)
	tform[[1L]] <- quote(stats::model.frame)  # change the function called
	#print(tform)
	#take the !is.list(formula) branch in if statement starting on line 69
	multiform <- F
	covlist <- NULL
	dformula <- formula
	
	special <- c("strata", "tt", "frailty", "ridge", "pspline")
    	tform$formula <- if(missing(data)) terms(formula, special) else
                                      terms(formula, special, data=data)

	if (!is.null(attr(tform$formula, "specials")$tt)) {
        coxenv <- new.env(parent= environment(formula))
        assign("tt", function(x) x, envir=coxenv)
        environment(tform$formula) <- coxenv
    	}
	mf <- eval(tform, parent.frame())
	Terms <- terms(mf)
	n <- nrow(mf)
    	Y <- model.response(mf) #807953 rows as expected
	isSurv2 <- inherits(Y, "Surv2")
	id <- model.extract(mf, "id")
       istate <- model.extract(mf, "istate")
	# print(isSurv2) #as expected, FALSE
	
	type <- attr(Y, "type") #right
	data.n <- nrow(Y)
	multi <- FALSE

	if (control$timefix) Y <- aeqSurv(Y)
	if (length(attr(Terms, 'variables')) > 2) { # a ~1 formula has length 2
        ytemp <- terms.inner(formula[1:2])
		#print(ytemp)
        suppressWarnings(z <- as.numeric(ytemp)) # are any of the elements numeric?
        ytemp <- ytemp[is.na(z)]  # toss numerics, e.g. Surv(t, 1-s)
		#print(z)
        xtemp <- terms.inner(formula[-2])
		#print(xtemp)
        if (any(!is.na(match(xtemp, ytemp))))
            warning("a variable appears on both the left and right sides of the formula")
    	}
	strats <- attr(Terms, "specials")$strata
	hasinteractions <- FALSE
   	dropterms <- NULL

	istrat <- NULL

	timetrans <- attr(Terms, "specials")$tt
	#return(Y)
	if (length(timetrans)) {
		timetrans <- untangle.specials(Terms, 'tt')
		ntrans <- length(timetrans$terms)
		tt <- list(tt)
		
		#is.list(tt) checks are non-problematic
		if (ncol(Y) == 2) {
			if(length(strats) == 0) {
				sorted <- order(-Y[,1], Y[,2])
				newstrat <- rep.int(0L, nrow(Y))
				newstrat[1] <- 1L
			} 
		}
	}
		
	
	return(newstrat)

}

#inf_survcovars_combined$age_group2 <- as.character(inf_survcovars_combined$age_group)

#benchmarking
ro_inds <- order(-inf_table_combined$time, inf_table_combined$infected)
inf_table_ro <- inf_table_combined[ro_inds,]

rownames(inf_table_ro) <- 1:nrow(inf_table_ro)
inf_table_ro <- make_dummy(inf_table_ro)
inf_table_uo <- make_dummy(inf_table_combined)

#sample_size <- 10^5

sample_size <- sum(inf_table_combined$count) #full sample

samp_inds <- sample(1:sum(inf_table_combined$count), size = sample_size)
samp_inds <- sort(samp_inds)

sample_ro <- expand_counts(inf_table_ro)[samp_inds, ]
sample_uo <- expand_counts(inf_table_uo)[samp_inds, ]

rownames(sample_ro) <- 1:sample_size
rownames(sample_uo) <- 1:sample_size


transform <- function(x, t, ...){
				#vac <- (x == c("Janssen", "Moderna", "Pfizer")) * t
				cbind(Janssen = (x == "Janssen"), Moderna = (x == "Moderna"), Pfizer = (x == "Pfizer"))
			}
cox_ro <- coxph(Surv(time, infected) ~ age_group + vaccine, data = sample_ro)

cox_ro_alt <- coxph(Surv(time, infected) ~ age_group + is_Janssen + is_Moderna + is_Pfizer, data = sample_ro)

#spot check the results match
cox_ro$coef
cox_ro_alt$coef

-cox_ro$coef[5] + cox_ro$coef[3] #for Moderna
cox_ro_alt$coef[4]
#checks out

begin_ro <- Sys.time()
cox_ro <- coxph(Surv(time, infected) ~ age_group + is_Janssen + is_Moderna + is_Pfizer + time_Janssen + time_Moderna + time_Pfizer, data = sample_ro)
end_ro <- Sys.time()

begin_uo <- Sys.time()
cox_uo <- coxph(Surv(time, infected) ~ age_group + is_Janssen + is_Moderna + is_Pfizer + time_Janssen + time_Moderna + time_Pfizer, data = sample_uo, iter.max = 10000)
end_uo <- Sys.time()

print(end_ro - begin_ro)
print(end_uo - begin_uo)

#pre-sorting does not make a difference

#full model
cox_uo <- coxph(Surv(time, infected) ~ age_group + is_Janssen + is_Moderna + is_Pfizer + time_Janssen + time_Moderna + time_Pfizer, data = sample_uo)

#reduced model (intercept coerced)
cox_uo_alt <- coxph(Surv(time, infected) ~ age_group + time_Janssen + time_Moderna + time_Pfizer, data = sample_uo)

x <- seq(1, 34, by = 0.001)
y_Janssen <- exp(cox_uo$coef[3] + cox_uo$coef[6]*x)
y_Moderna <- exp(cox_uo$coef[4] + cox_uo$coef[7]*x)
y_Pfizer <- exp(cox_uo$coef[5] + cox_uo$coef[8]*x)

plot(x, y_Janssen, type = "l", col = "orange")
points(x, y_Moderna, type = "l", col = "blue")
points(x, y_Janssen, type = "l", col = "green")
points(inf_hr_all$t, inf_hr_all$Janssen, col = "orange")
points(inf_hr_all$t, inf_hr_all$Moderna, col = "blue")
points(inf_hr_all$t, inf_hr_all$Pfizer, col = "green")


cox_trace <- coxph_part(inf_surv_combined ~ tt(vaccine) + age_group2,
			data = inf_survcovars_combined,
			tt = function(x, t, ...){
				vac <- (x == c("Janssen", "Moderna", "Pfizer")) * t
				cbind(Janssen = (x == "Janssen"), Moderna = (x == "Moderna"), Pfizer = (x == "Pfizer"))
			})

test_data <- expand_counts(make_dummy(inf_table_combined[(inf_table_combined$time > 5), ]))

cox_uo <- coxph(Surv(time, infected) ~ age_group + is_Janssen + is_Moderna + is_Pfizer + time_Janssen + time_Moderna + time_Pfizer, data = test_data)

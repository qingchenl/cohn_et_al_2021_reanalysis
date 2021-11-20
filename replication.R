library(survival)
Surv_wrapper <- function(time, status, counts) {
	return(Surv(rep(time, times = counts), rep(status, times = counts)))
}

inf_table_all_unvac <- inf_table_all[inf_table_all$vaccine == "Unvaccinated",]
inf_surv_all_unvac <- with(inf_table_all_unvac, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_unvac)

inf_table_all_jj <- inf_table_all[inf_table_all$vaccine == "Janssen",]
inf_surv_all_jj <- with(inf_table_all_jj, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_jj)

inf_table_all_moderna <- inf_table_all[inf_table_all$vaccine == "Moderna",]
inf_surv_all_moderna <- with(inf_table_all_moderna, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_moderna)

inf_table_all_pb <- inf_table_all[inf_table_all$vaccine == "Pfizer",]
inf_surv_all_pb <- with(inf_table_all_pb, Surv_wrapper(time, infected, counts))
plot(inf_surv_all_pb)
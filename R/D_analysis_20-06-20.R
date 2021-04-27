#Running analysis


####calculating arrival age####
###calculating with new function 
system.time(db_measures_Tyrannidae<- DivB_metrics(W = W, tree = tree, ancestral.area = ancestral.area, biogeo = biogeo, 
                                      diversification = c("jetz", "freck"), PD = TRUE, PE = TRUE, 
                                      age.arrival = TRUE, age.no.ancestor = NA, 
                                      dispersal.from = TRUE, ED.type = "equal.splits"))


max(do.call(data.frame,lapply(db_measures_Tyrannidae$Db_comm_metrics, function(x) replace(x, is.infinite(x),NA)))[,2], na.rm = T)
max(db_measures_Tyrannidae$Db_comm_metrics[, 1])
age_arrival<- diversification.assembly(W= W, tree= tree, ancestral.area= ancestral.area, biogeo= binogeo)
ages<- age_arrival$age_arrival
saveRDS(ages, here::here("R", "agesResult.rds"))

#mean age for each assemblage
mean_age<- apply(ages, 1, function(x) mean(x[which(x != 0)])) #mean arrival age for each assemblage

######Calating NRI #######
dis <- cophenetic(tree) #cophenetic distance matrix
org<- SYNCSA::organize.syncsa(comm=W, phylodist=dis) #organizing matrices
nri <- ses.mpd(org$community,org$phylodist,null.model="taxa.labels",runs = 999) #nri calculation
saveRDS(nri, here::here("R", "nriRes.rds")) #saving nri results
ses_mpd_res<- nri$mpd.obs.z
ses.mpd_res_noNA<- ifelse(!is.na(ses_mpd_res), ses_mpd_res, 0)

#binding Age and NRI
anova_data<- data.frame(mean_age= mean_age, NRI= nri_res_noNA, local= temp_trop)
anova_data_ses<- data.frame(mean_age= mean_age, NRI= ses.mpd_res_noNA, local= temp_trop)
###Anova####
res_anovaAge<- anova.1way(mean_age~temp_trop, data = anova_data_ses, nperm=10000) #anova with age
res_anovaNRI<- anova.1way(NRI~temp_trop, data = anova_data_ses, nperm=10000) #anova with NRI values

db_measures_Tyrannidae$Db_comm_metrics

#Running analysis


####calculating arrival age####
age_arrival<- diversification.assembly(W= W, tree= tree, ancestral.area= ancestral.area, biogeo= biogeo)
ages<- age_arrival$age_arrival

#mean age for each assemblage
mean_age<- apply(ages, 1, function(x) mean(x[which(x != 0)])) #mean arrival age for each assemblage

######Calating NRI #######
dis <- cophenetic(tree) #cophenetic distance matrix
org<- SYNCSA::organize.syncsa(comm=W, phylodist=dis) #organizing matrices
nri <- ses.mpd(org$community,org$phylodist,null.model="taxa.labels",runs = 999) #nri calculation
saveRDS(nri, here::here("R", "nriRes.rds")) #saving nri results
nri_res <- -1 * nri$mpd.obs.z #nri values
nri_res_noNA<- ifelse(!is.na(nri_res), nri_res, 0)

#binding Age and NRI
anova_data<- data.frame(mean_age= mean_age, NRI= nri_res_noNA, local= temp_trop)

###Anova####
res_anovaAge<- anova.1way(mean_age~temp_trop, data = anova_data, nperm=10000) #anova with age
res_anovaNRI<- anova.1way(NRI~temp_trop, data = anova_data, nperm=10000) #anova with NRI values


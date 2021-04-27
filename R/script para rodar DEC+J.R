##############################################################
#SCRIPT PARA RODAR DEC+J EM 3 CLADOS DE CANIDAE
##############################################################
library(optimx)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)
library(ape)
library(phytools)




#install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")

#Now let's load some of the source files necessary to run BioGeoBears:
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)





#And finally, let's get the directory of the example data that installed with BioGeoBears:
#extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))

#####################################################################################
###########****************************************************###########################
################################IMPORTANTE######################################
#####################################################################################


#diretório PC do vale
extdata_dir= setwd("C:\\Users\\window\\Desktop\\Doutorado\\CAP1- Reconstrução de areas ancestrais e modelagem do passado\\ANALISE\\Reconstrução de areas ancestrais")

#diretório PC de casa
extdata_dir= setwd("C:\\Users\\lucas\\OneDrive\\Documentos\\MEGA\\UFRGS\\DOUTORADO\\CAP1- Reconstrução de areas ancestrais e modelagem do passado\\ANALISE\\Reconstrução de areas ancestrais")


getwd()
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################



#tem que chamar os pacotes de novo 


#load the tree
tree_file_name = np(paste(addslash(extdata_dir), "canis.newick", sep=""))
tr = read.tree(tree_file_name)

#Let's plot the tree:
plot(tr)
axisPhylo()




#Now we need to load data on the geographic range of each extant species in our phylogey.
geo_file_name = np(paste(addslash(extdata_dir), "canis_classificacao_3.data", sep=""))
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)

#Let's take a look at the geographic range data
tipranges


# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 4



####################################################
####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=F)
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=F)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=F)

# Large numbers of areas have problems:
numstates_from_numareas(numareas=8, maxareas=8, include_null_range=F)

# ...unless you limit the max_range_size:
numstates_from_numareas(numareas=8, maxareas=4, include_null_range=F)
####################################################
####################################################



#################################################################
###########DEC

#First we will run the analysis using the DEC model. The DEC is the default model in
#BioGeoBears, so it is relatively straightforward to setup.
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

#Give BioGeoBEARS the location of the example input files:
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name


#And let's configure our analysis. If you are running this on your own dataset be sure to
#adjust the maximum range size (that is the maximum number of areas a lineage can inhabit
#at an given point - usually this is just the total number of areas in your dataset).
BioGeoBEARS_run_object$max_range_size = 4
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = F
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE


#Now we are ready to run the analysis:
results_DEC = bears_optim_run(BioGeoBEARS_run_object)



#################################################################
###########DEC+J

#Now we will run the analysis using the DEC+J model. This include the "jump" parameter
#for long distance dispersal / founder speciation events.
BioGeoBEARS_run_object = define_BioGeoBEARS_run()


#Give BioGeoBEARS the location of the example input files:
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name

trfn=tree_file_name
geogfn=geo_file_name

#If you are running this on your own dataset be
#sure to adjust the maximum range size. These settings are all the same as they were for the
#DEC model above.
BioGeoBEARS_run_object$max_range_size = 4
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = F
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE


#And now finally, set up DEC+J model. We'll use the maximum likelihood parameter value
#estimates from the DEC analysis we already ran to get good starting points for the hillclimbing
#heuristic.
dstart = results_DEC$outputs@params_table["d","est"]
estart = results_DEC$outputs@params_table["e","est"]

#Set the starting values in our analysis object:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

#The DEC is a 2-parameter model nested within the 3-parameter DEC+J, so we need to add
#J as a new free parameter to estimate. We also need to assign it an initial value.
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

#Now we are ready to run the analysis:
results_DECJ = bears_optim_run(BioGeoBEARS_run_object)




###########################PLOTS  DEC and DEC+J ################################################


######Plot the DEC ancestral states. These plots will be written to a PDF file.
analysis_titletxt = "DEC on Caninae"
results_object = results_DEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=F, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                         cornercoords_loc=scriptdir, include_null_range=F, tr=tr, tipranges=tipranges)



#######Plot the DEC+J ancestral states to the PDF file.
analysis_titletxt ="DEC+J on Caninae"
results_object = results_DECJ
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=F, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.5, tipcex=0.7, statecex=0.3, splitcex=0.3, titlecex=0.8, plotsplits=TRUE,
                         cornercoords_loc=scriptdir, include_null_range=F, tr=tr, tipranges=tipranges)





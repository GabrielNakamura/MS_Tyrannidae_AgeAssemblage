#Tudo que ? igual eu tiro e so carrego
#O que for especifico da analise vai dentro do for

#igual

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4

BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)


# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

#aqui diferente
objeto_resp<- vector(mode = "list",length= numerodearvores )

for(i in 1:numerodearvores){
  #DEC model
  runslow = TRUE
  resfn = "Psychotria_DEC_M0_unconstrained_v1.Rdata"
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    save(res, file=resfn)
    resDEC = res
  } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
  }
  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  resfn = "Psychotria_DEC+J_M0_unconstrained_v1.Rdata"
  runslow = TRUE
  #DEC+J model
  if (runslow)
  {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
    
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resDECj = res
  } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
  }
  
  #DIVALIKE model parameters
  # Set up DIVALIKE model
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  runslow = TRUE
  resfn = "Psychotria_DIVALIKE_M0_unconstrained_v1.Rdata"
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resDIVALIKE = res
  } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKE = res
  }
  testeDEC <- res$ML_marginal_prob_each_state_at_branch_top_AT_node
  
  # Set up DIVALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
  resfn = "Psychotria_DIVALIKE+J_M0_unconstrained_v1.Rdata"
  runslow = TRUE
  if (runslow)
  {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
    
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resDIVALIKEj = res
  } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKEj = res
  }
  
  #set up parameters for BAYAREA model
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # Check the inputs
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  runslow = TRUE
  resfn = "Psychotria_BAYAREALIKE_M0_unconstrained_v1.Rdata"
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resBAYAREALIKE = res
  } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKE = res
  }
  
  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  resfn = "Psychotria_BAYAREALIKE+J_M0_unconstrained_v1.Rdata"
  runslow = TRUE
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resBAYAREALIKEj = res
  } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKEj = res
  }
  teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
  teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
  row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
  restable = put_jcol_after_ecol(restable)
  # With AICs:
  AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
  restable = cbind(restable, AICtable)
  restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
  restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
  restable_AIC_rellike
  model<- min(restable_AIC_rellike)
  if(model="nome da coluna"){
    objeto_resp[[i]]<- resDEC
  }
  if(model= "outro nome"){
    objeto_resp[[i]]<- DECJ
  }
  
}

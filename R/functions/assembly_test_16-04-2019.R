####migration assembly function with diversification####
#####input arguments#####
# W = composition matrix, in lines the sample units and columns species
# tree = phylogenetic tree from class phylo;
# ancestral.area = data.frame. Lines with node names from tree and one column containing the ancestral area of each node;
# biogeo = data.frame. Lines with names of sample units, must be the same names that in W and one column containing the Ecoregion that each sample unit belongs.

######outputs#######
# Ecoregion.per.node = a matrix object. lines containig the nodes and the colums the species. Each cell present the Ecorregion that each ancestral of species i was found;
# Diversification.Assembly_Jetz = numeric. Biome Diversification values calculated according to a modifyied version of Jetz's metric of diversification;
# Diversification.Assembly_Freck =  numeric. Biome Diversification values calculated according to a modifyied version of Freckleton's metric of diversification;
# Diversification_total = matrix. Diversification values calculate for species according to Jetz and Freckleton metrics.

library(geiger)
library(ape)
library(ade4)
library(phylobase)
library(picante)
library(sjmisc)
migration.assembly<-function(W, tree, ancestral.area, biogeo){
  s <- length(tree$tip.label)
  n <- tree$Nnode
  names_spComm<- colnames(W)
  PDt<- picante::pd(samp = W, tree = tree, include.root = T)$PD
  EDtotal<- evol.distinct(tree = tree, type = "equal.splits")
  Jetz_total<- 1/EDtotal$w
  Freck_total<- numeric(length = length(tree$tip.label))
  names(Freck_total)<- tree$tip.label
  type<- "equal.splits" #argument to compute equal splits metrics in ER calculation - internal function
  T_freck<- max(cophenetic(tree))/2
  for(l in 1:length(tree$tip.label)){
    nodes_Freck<- .get.nodes(tree, tree$tip.label[l]) 
    #nodes_Freck<- nodes_Freck[1:(length(nodes_Freck) - 1)] #extract the internal nodes that make the path to species j, from the tip to root
    #T_freck<- sum(tree$edge.length[which(tree$edge[, 
    #                                               2] %in% nodes_Freck)], tree$edge.length[which.edge(tree, 
    #                                                                                                  paste("s",l,sep=""))])
    Freck_total[l]<- length(nodes_Freck)/T_freck
  }
  names(Jetz_total)<- tree$tip.label
  spxnode <- matrix(0,s,n)
  ages<-abs(node.depth.edgelength(phy = tree)
            -max(node.depth.edgelength(tree)))[-c(1:length(tree$tip.label))]
  ages<- data.frame(age= ages)
  rownames(ages)<- paste("N",(s+1):(s+(s-1)),sep="")
  tree_test1 <- suppressWarnings(phylo4(tree)) # precisa converter para esse formato para a funcao ancestors() funcionar
  for (i in 1:s){
    #i=1
    nos_ancestrais <- ancestors(tree_test1,i)
    nos_ancestrais <- nos_ancestrais - s  
    spxnode[i,nos_ancestrais] <- 1
  }
  rownames(spxnode)<-tree$tip.label #names for columns 
  colnames(spxnode)<- paste("N",(s+1):(s+(s-1)),sep="") #names for Nodes
  spxnode<- t(spxnode) #mesmo que Node
  AS<-spxnode #Ancestral State matrix 
  for(i in 1:nrow(spxnode)){
    pres<-which(spxnode[i,]==1)
    AS[i,pres]<-as.character(ancestral.area[i,1]) #matriz - a informacao em cada celula representa o estado do ancestral de cada especie
  }
  for(i in 1:nrow(AS)){
    zero<-which(AS[i,]==0)
    AS[i,zero]<-NA #NAs indicam nos que nao contem a especie 
  }
  matrix_XJetz<-W #matrix to receive the results of local Jetz metric
  matrix_XFreck<- W #matrix to receive the results of local Freckleton metric
  age_arrival<-  W
  node_local<- vector(mode = "list") 
  PD_local<- matrix(NA, nrow= nrow(W), ncol= 1, dimnames= list(rownames(W), "PD_local"))
  #i=3
  for(i in 1:nrow(W)){
    pres<- which(W[i,]==1)
    pres<- names_spComm[pres]
    nodes_species<- vector(mode= "list")
    tip_edgeLength<- numeric()
    #j=2
    for(j in 1:length(pres)){
      nodes_sp<- AS[,pres[j]][!is.na(AS[,pres[j]])]
      nodes_sp<- nodes_sp[length(nodes_sp):1] #nodes for species j in community i
      #if(biogeo[i,1]!=nodes_sp[1]) #arrumar aqui para abranger mais de uma letra
      if(str_contains(nodes_sp[1],biogeo[i,1])!=TRUE)
      { #check if the most recent ancestor was in the same bioregion that the observed species
        #matrix_XJetz[i, pres[j]]<- ages[names(nodes_sp[1]),1]/2 #if TRUE put as the arrival age at the half of the value of the most recent ancestor and zero 
        matrix_XJetz[i,pres[j]]<- 0.00001 
        matrix_XFreck[i,pres[j]]<- 0.00001
      } else{
        nodes_Freck<- .get.nodes(tree, pres[j]) #get the internal nodes from tip to root for species j
        nodes_Freck_internal<- nodes_Freck[1:(length(nodes_Freck))] #organize the internal nodes for species j
        Div_Freck<- length(nodes_Freck)/T_freck #modified equation to calculate local diversity based in Equation 4 Freckleton et al paper  Am.Nat 2008 
        nodes_all<- numeric(length = length(nodes_sp)) #test if all ancestors of species j are in the same ecoregion that local i
        for(m in 1:length(nodes_sp)){
          nodes_all[m]<- str_contains(nodes_sp[m], biogeo[i,1])
        }
        if(all(nodes_all==1)){ #if all ancestors of species j are in the same ecoregion of local 1 this will be TRUE
          x<- names(nodes_sp[length(nodes_sp)]) # if TRUE, take node 8 as the reference node for calculation of local diversification
        }else
          {
          if(length(nodes_sp)==1){
            x<- names(nodes_sp) # for the case in which the specie do not present internal branches
          }
          #k=1
          for(k in 1:length(nodes_sp)){ #find for the most ancient ancestor of species j that was present in the local i
            if(str_contains(nodes_sp[k],biogeo[i,1])!=TRUE){
              x<- names(nodes_sp)[k-1] # take the name of the node (ancestor) of species j that was present in the local i
              break
            } 
            #if(biogeo[i,1]!=nodes_sp[k]){
            #  x<- names(nodes_sp)[k-1] # case in which species present internal branches
            #  break
            #}
          } #extract the most ancient ancestral that was presented at Ecoregion of species j in local i
        }
        nodes_div<- sort(c(nodepath(tree, from = as.numeric(strsplit(x, split = "N")[[1]][2]), 
                                    to = which(tree$tip.label==pres[j])))[-length(nodepath(tree, from = as.numeric(strsplit(x, split = "N")[[1]][2]), to = which(tree$tip.label==pres[j])))],decreasing = TRUE) #organize the path of the most ancient ancestral that was presented at ecoregion of local i 
        nodes_species[[j]]<- nodes_div #nodes for species j, this was done to compute PD local
        if(length(nodes_div)==1){
          internal.brlen_div <- tree$edge.length[which(tree$edge[, 
                                                                 2] %in% nodes_div)] #branch lenghts (in times) for internal branch lengts of the most ancient ancestral 
          #ED_div <- internal.brlen_div
        } else{
          nodes_div <- nodes_div[1:(length(nodes_div) - 1)] #internal nodes that form the path from most ancient ancestral that was presented at Ecoregion of local i to species j
          internal.brlen_div <- tree$edge.length[which(tree$edge[, 
                                                                 2] %in% nodes_div)] #ages for internal nodes of nodes_div object
        }
        Div_Freck_local<- (length(nodes_div)+1)/T_freck #modifyed equation 4 from Freckleton et al (2008) considering only the nodes that diversified in ecoregion of local i, plus one is only a correction of the previous step
        if (length(internal.brlen_div) != 0) { #starts the calculation for ED measure from Redding et al.
          internal.brlen_div <- internal.brlen_div * switch(type, equal.splits = sort(rep(0.5, 
                                                                                          length(internal.brlen_div))^c(1:length(internal.brlen_div))), 
                                                            fair.proportion = {
                                                              for (j in 1:length(nodes_div)) {
                                                                sons <- .node.desc(tree, nodes_div[j])
                                                                n.descendents <- length(sons$tips)
                                                                if (j == 1) portion <- n.descendents else portion <- c(n.descendents, 
                                                                                                                       portion)
                                                              }
                                                              1/portion
                                                            })
          
        }
        tip_edgeLength[j]<- tree$edge.length[which.edge(tree, pres[j])]
        
        ED_div <- sum(internal.brlen_div, tree$edge.length[which.edge(tree, 
                                                                      pres[j])]) #Local ED - modifyed ED considering only the edges of ancestros inside the biogeo of local i for species j
        EDtotal_spp<- EDtotal$w[which(EDtotal$Species==pres[j])]
        #diff_total<- EDtotal$w[which(EDtotal$Species==pres[j])]-ED_div
        JetzTot<- 1/EDtotal$w[which(EDtotal$Species==pres[j])] #Diversificatio calculated according to Jetz for species j
        JetzLocal<- (JetzTot*(ED_div/EDtotal_spp)) #Modified Jetz to calculate only for local diversification
        matrix_XJetz[i, pres[j]]<- JetzLocal #Jetz local diversification
        matrix_XFreck[i, pres[j]]<- Div_Freck_local #Freckleton local diversification
        #matrix_X[i,pres[j]]<-ages[x,1] 
        age_arrival[i,pres[j]]<- ages[x,] #recebe a idade de chegada
      }
    }
    unique_internal_blength<- tree$edge.length[which(tree$edge[, 
                                                               2] %in% unique(nodes_species))]
    PD_local[i,]<- sum(unique_internal_blength, sum(tip_edgeLength))
  }
  Diversification_table<- data.frame(Jetz_total= Jetz_total, Freckleton_total= Freck_total) # original Jetz and Freckleton measures of diversification
  return(list(Ecoregion.per.node=AS,Diversification.Assembly_Jetz=matrix_XJetz, Diversification.Assembly_Freck= matrix_XFreck, Diversification_total= Diversification_table, 
              PD_local= PD_local, PD_total= PDt, age_arrival= age_arrival))
}
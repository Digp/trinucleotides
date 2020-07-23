library(veriNA3d)
library(RNeo4j)
library(minipack)

#source("validation_function-working.R")
#source("Extract_ntinfo_function.R")

graph = startGraph("http://localhost:7474/db/data/", username="neo4j", password="tabletennis1999")


# Retreive all PDB IDs from the Protein Dat Bank and take only those "nuc" and "prot-nuc"

#alllist <- queryEntryList(justIDs = FALSE)
#nuclist <-  alllist$pdbID[alllist$type %in% c("nuc", "prot-nuc")]


#Working with small PDBID dataset
nuclist <- c("1bau", "100d")
 

#Loop over all the PDBID
for (i in 1:length(nuclist)){
  
  # Get the PDB ID
  pdbID <- nuclist[i]
  
  #Perform validation analysis on pdbID
  data <- validation(pdbID)
  
  #Extract ntinfo function
  ntinfo <- extract_ntinfo(pdbID)
  localenv <- ntinfo$localenv
  localenv <- as.data.frame(localenv)
  ntinfo$eta[is.na(ntinfo$eta)] <- NaN
  ntinfo$theta[is.na(ntinfo$theta)] <- NaN
  ntinfo$kappa[is.na(ntinfo$kappa)] <- NaN
  ntinfo$alpha[is.na(ntinfo$alpha)] <- NaN
  ntinfo$beta[is.na(ntinfo$beta)] <- NaN
  ntinfo$gamma[is.na(ntinfo$gamma)] <- NaN
  ntinfo$delta[is.na(ntinfo$delta)] <- NaN
  ntinfo$epsilon[is.na(ntinfo$epsilon)] <- NaN
  ntinfo$zeta[is.na(ntinfo$zeta)] <- NaN
  ntinfo$nu0[is.na(ntinfo$nu0)] <- NaN
  ntinfo$nu1[is.na(ntinfo$nu1)] <- NaN
  ntinfo$nu2[is.na(ntinfo$nu2)] <- NaN
  ntinfo$nu3[is.na(ntinfo$nu3)] <- NaN
  ntinfo$nu4[is.na(ntinfo$nu4)] <- NaN
  ntinfo$chi[is.na(ntinfo$chi)] <- NaN
  ntinfo$pu_amp[is.na(ntinfo$pu_amp)] <- NaN
  ntinfo$pu_phase[is.na(ntinfo$pu_phase)] <- NaN
  ntinfo$Dp[is.na(ntinfo$Dp)] <- NaN
  
  data <- as.data.frame(cbind(data, ntinfo))
  

  # Create nodes and releationships directly from a dataframe

  x <- createNode(graph, "PDBID", name = pdbID)
  
  for (j in 1:length(data$seq)){
    
    p <- createNode(graph, "NUCLEOTIDE", name = data[j, ]$id_dssr)
    
    createRel(x, "CONTAINS", p)
    #createRel(e, "IDENTIFIER", p)
    #createRel(e, "TRINUCLEOTIDE", l)
    p <- updateProp(p, localenv = data[j, ]$localenv)
    p <- updateProp(p, seq = data[j, ]$seq)
    p <- updateProp(p, clash = data[j, ]$clashes)
    p <- updateProp(p, pucker = data[j, ]$pucker_outlier)
    p <- updateProp(p, bond_lengths = data[j, ]$bond_lengths)
    p <- updateProp(p, suite = data[j, ]$suite_outlier)
    p <- updateProp(p, bond_angles = data[j, ]$bond_angles)
    p <- updateProp(p, chirals = data[j, ]$chirals)
    p <- updateProp(p, planes = data[j, ]$planes)
    p <- updateProp(p, rsrz = data[j, ]$rsrz)
    p <- updateProp(p, theta = data[j, ]$theta)
    p <- updateProp(p, eta = data[j, ]$eta)
    p <- updateProp(p, kappa = data[j, ]$kappa)
    p <- updateProp(p, alpha = data[j, ]$alpha)
    p <- updateProp(p, beta = data[j, ]$beta)
    p <- updateProp(p, gamma = data[j, ]$gamma)
    p <- updateProp(p, delta = data[j, ]$delta)
    p <- updateProp(p, epsilon = data[j, ]$epsilon)
    p <- updateProp(p, zeta = data[j, ]$zeta)
    p <- updateProp(p, pu_amp = data[j, ]$pu_amp)
    p <- updateProp(p, pu_phase = data[j, ]$pu_phase)
    p <- updateProp(p, nu0 = data[j, ]$nu0)
    p <- updateProp(p, nu1 = data[j, ]$nu1)
    p <- updateProp(p, nu2 = data[j, ]$nu2)
    p <- updateProp(p, nu3 = data[j, ]$nu3)
    p <- updateProp(p, nu4 = data[j, ]$nu4)
  }
}

for (i in 1:length(names(data)))
#Some queries
# MATCH (p{name:"1bau"}) -[:CONTAINS]->(n{clash:true}) RETURN p,n
# MATCH (p) -->() RETURN p,n
# MATCH (n) DETACH DELETE n



# UKCI2022_sepis.R

##### SEPSIS SEQUENCE ANOMALIES ------------------

data(sepsis)  # from seqDetect package but use my RData file
load("sepsis.seq.RData")  # The data in sepsis took a lot of preprocessing to get in a fit state for PST !!
set.seed(101) # Set Seed so that same sample can be reproduced in future
# Now randomly select 75% of data as sample from total 'n' rows of the data  
seqsample <- sample.int(n = nrow(sepsis.seq), size = floor(.75*nrow(sepsis.seq)), replace = FALSE)
sepsis.seq.train <- sepsis.seq[seqsample, ]
sepsis.seq.test  <- sepsis.seq[-seqsample, ]

bob <- data.frame(lapply(sepsis.seq, as.character), stringsAsFactors=FALSE)
sepSTR <- NULL
allSTR <- ""
tempSTR <- ""

for(i in 1:200){
tempSTR <- ""
bob[i,] <- str_replace(bob[i,], "ERRegistration", "A")
bob[i,] <- str_replace(bob[i,], "Leucocytes", "B")
bob[i,] <- str_replace(bob[i,], "CRP", "C")
bob[i,] <- str_replace(bob[i,], "LacticAcid", "D")
bob[i,] <- str_replace(bob[i,], "ERTriage", "E")
bob[i,] <- str_replace(bob[i,], "ERSepsisTriage", "F")
bob[i,] <- str_replace(bob[i,], "IVLiquid", "G")
bob[i,] <- str_replace(bob[i,], "AdmissionNC", "H")
bob[i,] <- str_replace(bob[i,], "IVAntibiotics", "I")
bob[i,] <- str_replace(bob[i,], "AdmissionIC", "X")
bob[i,] <- str_replace(bob[i,], "ReleaseA", "1")
bob[i,] <- str_replace(bob[i,], "ReleaseB", "2")
bob[i,] <- str_replace(bob[i,], "ReleaseC", "3")
bob[i,] <- str_replace(bob[i,], "ReleaseD", "4")
bob[i,] <- str_replace(bob[i,], "ReleaseE", "5")
bob[i,] <- str_replace(bob[i,], "ReturnER", "Z")

bob2 <- paste(bob[i,], sep="", collapse=NULL)
bob2 <- lapply(bob2, function(z){ z[!is.na(z) & z != ""]})
bob2 <- paste(bob2, sep="-", collapse="")
bob2 <- gsub('[%]', '', bob2)
bob2 <- gsub('[0]', '', bob2)
bob2 <- gsub("\\0","",bob2)

tempSTR <- bob2
allSTR <- strcat(allSTR,bob2)  # allSTR format is needed for motif searching

tempSTR <- strsplit(tempSTR, split = "")
tempSTR <- unlist(tempSTR)
tempSTR <- paste(tempSTR, collapse="-")
sepSTR[i] <- tempSTR# sepSTR format is needed for HMM    
}

motifs <- find_motifs(allSTR)
motifs_df <- as.data.frame(unlist(motifs[[2]]))
xtable(motifs_df)

# converts the sequence motif into a series of identifiers prefixed by "M" e.g. "M7" etc
# motifs[[2]] contains the list of motif names and order in which they appear
#motiforder <- motif_sequence(sepSTR,motifs[[2]])  
sepsis.motif.train <- TraMineR::seqdef(sepSTR[1:100])
sepsis.motif.test  <- TraMineR::seqdef(sepSTR[101:200])

# PHASE 1: Build a HMM using seqHMM on training data
hmm_sepsis <- build_hmm(observations = sepsis.motif.train, n_states = 5)
plot(hmm_sepsis,ncol.legend = 6)

hmm_fitted <- fit_model(hmm_sepsis)
# save important parameters
tr <- hmm_fitted$model$transition_probs
emiss <- hmm_fitted$model$emission_probs
init <- hmm_fitted$model$initial_probs

# PHASE 2: pass new sequence through the model
hmm_new_sepsis <- build_hmm(observations = sepsis.motif.test ,transition_probs = tr,
                     emission_probs = emiss, initial_probs = init)
plot(hmm_new)

# join the motifs in one string but separate them by "-" as required by seqdef() function
# e.g. "M2-M3-M5-M5-M6-M6-M6-M7-M5-M1"
motiforder <- paste(motiforder, collapse="-")

require("igraph")
set.seed(1234)
plot(hmm_sepsis,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))

plot(hmm_new,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))

# THE LAND OF WEIRD AND WONDERFUL PLOTS
subm <- seqsubm(sepsis.motif.train, method = "CONSTANT", with.miss = TRUE)
srhdist.om <- seqdist(sepsis.motif.train, method = "OM", sm = subm, with.miss = TRUE)
clusterward <- agnes(srhdist.om, diss = TRUE, method = "ward")
cluster3 <- cutree(clusterward, k = 3)
cluster3 <- factor(cluster3, labels = c("Type 1", "Type 2", "Type 3"))
#seqfplot(seq.SRH.test, group = cluster3, pbarw = T)
plot(clusterward,which.plots=TRUE)

# freq event subsequences
bf.seqe <- seqecreate(sepsis.motif.train)
fsubseq <- seqefsub(bf.seqe, min.support = 5)
head(fsubseq)

condition <- seqecontain(fsubseq, event.list = c("M25>M41"))
fsubseq[condition]
couts <- seqsubm(sepsis.motif.test, method = "TRATE")


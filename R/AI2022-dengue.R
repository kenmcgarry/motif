## dengue virus DNA sequence, this virus is composed of 10,735 DNA symbols, coding for 10 proteins
load("dengvirus.RData")

seq <- toupper(seq)
denseq <- paste(seq,collapse ="")
#head(denseq)

# for DNA important motifs are : start codon:"ATG"; stop codons: are "TGA", "TAA", and "TAG"
motifs <- find_motifs(substr(denseq,1,500))
motifs_df <- as.data.frame(unlist(motifs[[2]]))
#xtable(motifs_df)

# converts the sequence motif into a series of identifiers prefixed by "M" e.g. "M7" etc
# motifs[[2]] contains the list of motif names and order in which they appear
motiforder <- motif_sequence(denseq,motifs[[2]])  
# Now gather motifs into batches of 100, then join the motifs in one string but separate them by "-" 
motifbatch <- split(motiforder, ceiling(seq_along(motiforder)/20))
length(motifbatch)

for(i in 1:length(motifbatch)){
motifbatch[i] <- paste(unlist(motifbatch[i]), collapse="-")
}

seq.deng <- TraMineR::seqdef(unlist(motifbatch))
seq.deng.train <- TraMineR::seqdef(seq.deng[1:100,])
seq.deng.test  <- TraMineR::seqdef(seq.deng[101:161,])

# PHASE 1: Build a HMM using seqHMM on training data
hmm_deng <- build_hmm(observations = seq.deng.train, n_states = 4)
#plot(hmm_deng,ncol.legend = 6)

hmm_fitted <- fit_model(hmm_deng)
# save important parameters
tr <- hmm_fitted$model$transition_probs
emiss <- hmm_fitted$model$emission_probs
init <- hmm_fitted$model$initial_probs

# PHASE 2: pass new sequence through the model
hmm_new_deng <- build_hmm(observations = seq.deng.test ,transition_probs = tr,
                    emission_probs = emiss, initial_probs = init)
#plot(hmm_new)

require("igraph")
set.seed(1234)
plot(hmm_deng,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))

plot(hmm_new_deng,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))

# THE LAND OF WEIRD AND WONDERFUL PLOTS
subm <- seqsubm(seq.deng.train, method = "CONSTANT", with.miss = TRUE)
srhdist.om <- seqdist(seq.deng.train, method = "OM", sm = subm, with.miss = TRUE)
clusterward <- agnes(srhdist.om, diss = TRUE, method = "ward")
cluster3 <- cutree(clusterward, k = 3)
cluster3 <- factor(cluster3, labels = c("Type 1", "Type 2", "Type 3"))
#seqfplot(seq.SRH.test, group = cluster3, pbarw = T)
plot(clusterward,which.plots=TRUE)

# freq event subsequences
bf.seqe <- seqecreate(seq.deng.train)
fsubseq <- seqefsub(bf.seqe, min.support = 5)
head(fsubseq)

condition <- seqecontain(fsubseq, event.list = c("M25>M41"))
fsubseq[condition]
couts <- seqsubm(seq.deng.test, method = "TRATE")




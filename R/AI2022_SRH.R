# UKCI2022_SRH.R

## EXPERIMENT 3 - The self-rated health (SRH) data set contains sequences for 2612 respondents of a survey conducted 
## by the Swiss Household Panel (SHP), aged between 20 and 80 years at the start of the survey. The data is 
## organized into 11 variables of 2621 records (one reading for each person over the 11 years of the survey 1999-2009).
## This survey has missing data. The categories are: 
## G1	(very well)
## G2	(well)
## M	(so, so (average))
## B2	(not very well)
## B1	(not well at all)
## * (missing)

# SRH is a large dataset so we build the PST on 75% of training data
data("SRH", package = "PST")

set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now randomly select 75% of data as sample from total 'n' rows of the data  
#seqsample <- sample.int(n = nrow(SRH.seq), size = floor(.75*nrow(SRH.seq)), replace = FALSE)
#SRH.seq.train <- SRH.seq[seqsample, ]
#SRH.seq.test  <- SRH.seq[-seqsample, ]

bob <- data.frame(lapply(SRH.seq, as.character), stringsAsFactors=FALSE)
bob2 <- paste(unlist(bob), sep="-", collapse="")
bob2 <- gsub('[*]', 'M', bob2)  # replace "*" with "M"

motifs <- find_motifs(substr(bob2,1,1000))
motifs_df <- as.data.frame(unlist(motifs[[2]]))
#xtable(motifs_df)

for(i in 1:length(motifs_df)){
  motifs_df[i] <- paste(unlist(motifs_df[i]), collapse="-")
}

seq.SRH <- TraMineR::seqdef(unlist(motifs_df))
seq.SRH.train <- TraMineR::seqdef(seq.deng[1:20,])
seq.SRH.test  <- TraMineR::seqdef(seq.deng[21:40,])

# PHASE 1: Build a HMM using seqHMM on training data
hmm_srh <- build_hmm(observations = seq.SRH.train, n_states = 5)
#plot(hmm_srh,ncol.legend = 6)

hmm_fitted <- fit_model(hmm_srh)
# save important parameters
tr <- hmm_fitted$model$transition_probs
emiss <- hmm_fitted$model$emission_probs
init <- hmm_fitted$model$initial_probs

# PHASE 2: pass new sequence through the model
hmm_new_srh <- build_hmm(observations = seq.SRH.test,transition_probs = tr,
                     emission_probs = emiss, initial_probs = init)
#plot(hmm_new)

# igraph plot
require("igraph")
set.seed(1234)
plot(hmm_srh,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))

plot(hmm_new_srh,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))


## get data ready for motif discovery
SRH.str <- seqformat(SRH.seq,from = "STS", to = "STS")
SRH.str <- na.omit(SRH.str)
SRH.str <- paste(SRH.str, sep="", collapse=NULL) 
SRH.str <- paste(SRH.str,collapse = "")

motifs <- find_motifs(substr(SRH.str,start=1,stop=600))

# converts the sequence motif into a series of identifiers prefixed by "M" e.g. "M7" etc
# motifs[[2]] contains the list of motif names and order in which they appear
motiforder <- motif_sequence(seq1,motifs[[2]]) 

# THE LAND OF WEIRD AND WONDERFUL PLOTS
subm <- seqsubm(seq.SRH.train, method = "CONSTANT", with.miss = TRUE)
srhdist.om <- seqdist(seq.SRH.train, method = "OM", sm = subm, with.miss = TRUE)
clusterward <- agnes(srhdist.om, diss = TRUE, method = "ward")
cluster3 <- cutree(clusterward, k = 3)
cluster3 <- factor(cluster3, labels = c("Type 1", "Type 2", "Type 3"))
#seqfplot(seq.SRH.test, group = cluster3, pbarw = T)
plot(clusterward,which.plots=TRUE)

# freq event subsequences
bf.seqe <- seqecreate(seq.SRH.train)
fsubseq <- seqefsub(bf.seqe, min.support = 5)
head(fsubseq)

condition <- seqecontain(fsubseq, event.list = c("M25>M41"))
fsubseq[condition]
couts <- seqsubm(seq.SRH.test, method = "TRATE")

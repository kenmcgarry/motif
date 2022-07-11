# UKCI2022_markovchains.R
# https://stats.stackexchange.com/questions/35308/hidden-markov-models-and-anomaly-detection
# Email comms with Satu Helske (seqHMM package):-
# For indices of (dis)similarity between different sequences or sequence subsamples, the sequence 
# analysis literature has different approaches, for example this might be suitable: 
# Liao, T. F., & Fasang, A. E. (2021). Comparing Groups of Life-Course Sequences Using the Bayesian 
# Information Criterion and the Likelihood-Ratio Test. Sociological Methodology, 51(1), 44-85.
#
# Regarding your second question, do you wish to calculate the most probable path of hidden states 
# for the new sequence? For that you can use build_hmm and hidden_paths: pass your existing model 
# parameters and the sequences to the build_hmm function and then use the hidden_paths function for 
# finding the most probable sequences of hidden state for each unit (do not re-estimate the parameters 
# with the fit_model function in between).
#---------------------------------------------------------------

library(seqHMM)
# PHASE 1: Build a HMM using seqHMM on training data
hmm_model <- build_hmm(observations = seq.motifs.train1, n_states = 3)
plot(hmm_model,ncol.legend = 6)

hmm_fitted <- fit_model(hmm_model)
# save important parameters
tr <- hmm_fitted$model$transition_probs
emiss <- hmm_fitted$model$emission_probs
init <- hmm_fitted$model$initial_probs

# PHASE 2: pass new sequence through the model
hmm_new <- build_hmm(observations = seq.motifs.test1,transition_probs = tr,
                     emission_probs = emiss, initial_probs = init)
plot(hmm_new)

# use the hidden_paths function for finding the most probable sequences of hidden state for each unit
mpp <- hidden_paths(hmm_new)


# PHASE 3: What are the differences between original HMM model and new sequence data?
results <- compare_hmm(hmm_model,hmm_new)


# Simulating sequences of observed and hidden states given parameters of a hidden Markov model.
sim <- simulate_hmm(n_sequences = 20, 
                    initial_probs = init,
                    transition_probs = tr,
                    emission_probs = emiss,
                    sequence_length = 10)

sim$observations


# igraph plot
require("igraph")
set.seed(1234)
plot(hmm_biofam,
     layout = layout_nicely, pie = FALSE,
     vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
     edge.curved = FALSE, edge.width = 1,
     loops = TRUE, edge.loop.angle = -pi/8, edge.arrow.size = 0.5,
     trim = 0.01, label.signif = 3,
     xlim = c(-1, 1.3))

plot(hmm_biofam,
     layout = matrix(c(1, 2, 3, 4, 2,  1, 1, 1, 1, 0), ncol = 2),
     xlim = c(0.5, 4.5), ylim = c(-0.5, 1.5), rescale = FALSE,
     edge.curved = c(0, -0.8, 0.6, 0, 0, -0.8, 0),
     cex.edge.width = 0.8, edge.arrow.size = 0.5,
     legend.prop = 0.3, ncol.legend = 3,
     vertex.label.dist = 1.1, combine.slices = 0.02,
     combined.slice.label = "others (emission prob. < 0.02)")

# State inference for doing useful things: hidden_paths(), posterior_probs(), and forward_backward()
# The forward_backward function computes scaled forward and backward probabilities of a HMM.
# Compute forward and backward probabilities
#fb <- forward_backward(hmm_new)
# The most probable hidden state at time t
# given the observations up to time t for the first subject:
#apply(fb$forward_probs[, , 1], 2, which.max)

# Compute the most probable hidden state paths given the data and the model
mpp <- hidden_paths(hmm_new)

# Compute posterior probabilities
pb <- posterior_probs(hmm_new)
# Locally most probable states for the first subject:
pb[, , 1]


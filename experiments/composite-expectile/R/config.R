# Shared configuration for the experiment.
P0 <- 0.95      # weight switches on where the contamination is ~1% of the tail
P1 <- 0.98      # and is fully on where it is ~0.01%
N_PANEL <- 24   # composite Gauss-Legendre: 24 panels x 8 nodes = 192 levels,
N_GL    <- 8    # which holds the quadrature error below 0.15% of a fitted
                # return level (see scripts/98-validate.R)
W_FUN <- make_weight(P0, P1)
GRID  <- make_level_grid(P0, N_PANEL, N_GL)

N_OBS  <- 100   # sample size: a long annual-maximum record
N_REP  <- 2000  # Monte Carlo replicates
RETURN_PERIODS <- exp(seq(log(2), log(1000), length.out = 40))

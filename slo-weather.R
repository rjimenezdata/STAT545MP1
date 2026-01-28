## DATA ------------------------------------------------------------------------
library(tidyverse)

# read in noaa data for slo: daily precipitation
raw <- read_csv('https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/US1CASL0018.csv')

# 12ish years of data
raw$DATE |> range()

# daily precipitation
prcp <- raw$PRCP

# define "heavy": 75th percentile of positive precipitation
cutoff <- prcp %>%
  .[!is.na(.) & . > 0] %>%
  quantile(probs = 0.75, na.rm = TRUE) %>%
  as.numeric()

# precipitation state
prcp_state <- tibble(prcp = prcp) %>%
  mutate(
    state = case_when(
      is.na(prcp)        ~ NA_character_,
      prcp <= 0          ~ "none",
      prcp <= cutoff     ~ "light",
      TRUE               ~ "heavy"
    ),
    state = fct_relevel(state, "none", "light", "heavy")
  ) %>%
  pull(state) |>
  as.character()

# quick check
prcp_state

# share of days at each level -- note missing data
table(prcp_state, useNA = 'always') |> prop.table()

## CONDITIONAL LIKELIHOOD ESTIMATION (by hand) ---------------------------------

# use complete transitions only (remove transitions to/from NA)
origin <- prcp_state[-length(prcp_state)]
destination <- prcp_state[-1]
complete <- !is.na(origin) & !is.na(destination)
N <- table(origin[complete], destination[complete])
p.hat <- prop.table(N, margin = 1)
p.hat

# long-run behavior
v <- eigen(t(p.hat))$vectors[,1]
pi <- v/sum(v)
pi

## CONDITIONAL LIKELIHOOD ESTIMATION (package) ---------------------------------
library(markovchain)

# quick and dirty -- ignore NAs and stitch across gaps 
# (NOT CORRECT, but works okay if NA rate is low)
fit_approx <- prcp_state[!is.na(prcp_state)] |>
  markovchainFit(method = 'mle')
fit_approx$estimate

# split into contiguous blocks and remove NAs, then fit
# (CORRECT)
x <- prcp_state
is_miss <- is.na(x)
block_id <- cumsum(is_miss & !dplyr::lag(is_miss, default = FALSE))
blocks <- split(x, block_id)
blocks <- lapply(blocks, \(b) b[!is.na(b)])
blocks <- Filter(\(b) length(b) >= 2, blocks)
fit <- markovchainFit(data = blocks, method = 'mle')

# transition probabilities
fit$estimate
p.hat # for comparison

# stationary distribution
steadyStates(fit$estimate)
pi # for comparison

# k-step prediction
k <- 4
conditionalDistribution(fit$estimate^k, state = 'none')

# expected return times
1/pi

## DIAGNOSTICS -----------------------------------------------------------------

# data adequacy check
# (looking for sufficient effective sample size)
transition_counts <- fct_na_value_to_level(prcp_state, level = 'unobserved') |>
  createSequenceMatrix()
rowSums(transition_counts) # this is the main check
transition_counts # can also look for row sparsity

# compare stationary distribution to empirical marginals
# (mismatch indicates inhomogeneity or low mixing)
table(prcp_state) |> prop.table()
steadyStates(fit$estimate)

# simulation to check run lengths for dry spells
# (looking to see if the model reproduces time behavior of data)
set.seed(12526)
get_spells <- function(seq, target = "none") {
  r <- rle(seq)                  # run-length encoding
  r$lengths[r$values == target]  # lengths of runs where the value == target
}
obs_spells <- unlist(lapply(blocks, get_spells))
sim_spells <- unlist(lapply(blocks, function(b) {
  sim_seq <- markovchainSequence(
    n = length(b),
    markovchain = fit$estimate,
    t0 = b[1],
    include = TRUE
  )
  get_spells(sim_seq, target = "none")
}))
spell_summary <- function(spells, cutoff = 7) {
  spells <- spells[!is.na(spells)]
  if (length(spells) == 0) {
    return(c(n_spells = 0, mean = NA, median = NA, p90 = NA, p95 = NA, ge_cutoff = NA))
  }
  c(
    n_spells  = length(spells),
    mean      = mean(spells),
    median    = median(spells),
    p90       = unname(quantile(spells, 0.90)),
    p95       = unname(quantile(spells, 0.95)),
    ge_cutoff = mean(spells >= cutoff)
  )
}
rbind(
  observed  = spell_summary(obs_spells),
  simulated = spell_summary(sim_spells)
)

## FIT BY SEASON ---------------------------------------------------------------
# here we're checking for time homogeneity, which probably doesn't hold

# assign each day to a season
d <- as.Date(raw$DATE)
m <- as.integer(format(d, "%m"))
season <- ifelse(m %in% c(10, 11, 12, 1, 2, 3), "Oct-Mar", "Apr-Sep")

# create blocks within a season by breaking on missing data and nonconsecutive dates
make_blocks <- function(x_sub, d_sub) {
  n <- length(x_sub)
  if (n < 2) return(list())
  
  # TRUE where a NEW block should start (at position t)
  new_block <- logical(n)
  new_block[1] <- TRUE
  
  # break if current day is NA
  new_block[is.na(x_sub)] <- TRUE
  
  # break if previous day or current day is NA (avoid NA-adjacent transitions)
  new_block[c(FALSE, is.na(x_sub[-1]) | is.na(x_sub[-n]))] <- TRUE
  
  # break if dates aren't consecutive
  new_block[c(FALSE, diff(d_sub) != 1)] <- TRUE
  
  block_id <- cumsum(new_block)
  
  blks <- split(x_sub, block_id)
  blks <- lapply(blks, function(b) b[!is.na(b)])
  Filter(function(b) length(b) >= 2, blks)
}
idx_by_season <- split(seq_along(prcp_state), season)
blocks_by_season <- lapply(idx_by_season, function(idx) {
  make_blocks(x[idx], d[idx])
})
str(blocks_by_season)

# your turn: conditional mles for each season


# PoPLR

An implementation of PoPLR[1] that allows for missing 
values in first visits.

> [1] *Visual Field Progression in Glaucoma: Estimating the Overall Significance of Deterioration with Permutation Analyses of Pointwise Linear Regression (PoPLR)*
> Neil O'Leary; Balwantray C. Chauhan; Paul H. Artes
> *Investigative Ophthalmology & Visual Science*
> October 2012, Vol.53, 6776-6784. https://doi.org/10.1167/iovs.12-10049


# cpp version

Compiled on mac with 
    g++ -Wc++11-extensions -std=c++14 -I /opt/homebrew/include/ poplr.cpp

Compiled on ubuntu with 
    g++ -O6 -std=c++17 -I /home/aht/.cmdstan/cmdstan-2.35.0/stan/lib/stan_math/lib/boost_1.84.0 -fopenmp  poplr.cpp

# Input pr_tt_given_mt files

These 3 I used for the ARREST 10-2 paper and IPS 2026, so they are derived from the
ZEST procedure used in that paper
> Pr (tt in -10, -9, ..., 40 | mt in -1, 0, 1, .., 35) = pr_tt_given_mt[mt + 1][tt + 11]
> Normalised so each row sums to 1 (except mt == 36 which should not exist for this ZEST)

  1. arrest10_2_sd1.5_fp0_fn0_pr_tt_given_mt.csv
  2. arrest10_2_sd2_fp15_fn03_pr_tt_given_mt.csv
  3. arrest10_2_sd1.5_fp15_fn03_pr_tt_given_mt.csv

This is derived specifically for the ZEST used in hIPPOS (arrest11 with prior8).
See `~/src/hIPPOS/analysis/probs/make_probs.r`

  1. arrest11_fp15_fn03_pr_tt_given_mt.csv

# Input for $P(mt < 17 | tt)$

These 3 I used for the ARREST 10-2 paper and IPS 2026, so they are derived from the
ZEST procedure used in that paper
> Pr (tt in -10, -9, ..., 40 | mt in -1, 0, 1, .., 35) = pr_tt_given_mt[mt + 1][tt + 11]
> Normalised so each row sums to 1 (except mt == 36 which should not exist for this ZEST)

  1. arrest10_2_sd1.5_fp15_fn03_pr_mtlt17_given_tt.csv
  2. arrest10_2_sd2.0_fp15_fn03_pr_mtlt17_given_tt.csv

This is derived specifically for the ZEST used in hIPPOS (arrest11 with prior8).
See `~/src/hIPPOS/analysis/probs/make_probs.r`

  1. arrest11_fp15_fn03_pr_mtlt17_given_tt.csv

# Horizontal vs Vertical paritioning

It seems that using Horizontal partitioning might be too sensitive (low specificity) as you
are throwing away locations for each block. 
For example, in ip027, block 1 has the 52 24-2 locations, and block 2 only 3 locations.
The 3-location block has a low p-value, but if you add in just one of the 52 locations
with first visit removed, p jumps to 0.06.

But vertical paritioning throws away ability to re-order columns in time. 

Why not just permute the lot and use PLR with missing values?
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
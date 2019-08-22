-- the following function computes Buchsbaum-Eisenbud
-- multipliers of a free resolution based on B-E's paper
-- "Some structure theorems for finite free resolutions"
-- the only difference is that in order to make each multiplier
-- homogeneous we give its domain the appropriate degree
BEmult = F -> (
    -- ranks of the differentials
    n := length F;
    r := apply(1..n,i->rank F.dd_i);
    
    )
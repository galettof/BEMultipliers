-- the following function computes Buchsbaum-Eisenbud
-- multipliers of a free resolution based on B-E's paper
-- "Some structure theorems for finite free resolutions"
-- the only difference is that in order to make each multiplier
-- homogeneous we give its domain the appropriate degree
BEmult = F -> (
    -- length of the resolution
    n := length F;
    -- ranks of the differentials
    r := apply(1..n,i->rank F.dd_i);
    -- format of the resolution
    f := apply(n+1,i->rank F_i);
    -- the a will house the next multiplier
    -- the last multiplier is simply an exterior power
    a := exteriorPower(r_(n-1),F.dd_n);
    -- mults is the list of multipliers that will be returned
    mults := {a};
    -- i is a running index
    i := n-1;
    while (i > 0) do (
	-- the G's are rank one modules generated in the
	-- appropriate degree to make everything homogeneous.
	-- in fact this G is the domain of the next multiplier a
	G := exteriorPower(f_i,F_i);
	-- this code doesn't use the diagram of B-E's paper
	-- instead it uses its dual and accounts for degrees.
	-- e's are the exterior powers of the differentials
	e := (dual exteriorPower(r_(i-1),F.dd_i))**G;
	-- we need the iso of free modules Wedge^i F->Wedge^j F^*
	-- from the pairing Wedge^i F ** Wedge^j F->Wedge^(i+j) F
	-- where F is free has rank i+j.
	-- because of the way M2 computes exterior powers
	-- this iso is simply the matrix with 1's on the
	-- antidiagonal and 0's everywhere else.
	-- we obtain this as a permutation matrix
	I := (id_(target e))^(reverse toList(0..rank(target e)-1));
	-- get next multiplier by factoring as in dual diagram
	b := e // (I*a);
	-- now dualize back and fix the degrees
	a = dual (b ** (dual G));
	mults = prepend(a,mults);
	i = i-1;
	);
    return mults;
    )
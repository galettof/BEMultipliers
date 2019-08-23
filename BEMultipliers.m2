newPackage(
     "BEMultipliers",
     Version => "0.1",
     Date => "August 23, 2019",
     AuxiliaryFiles => false,
     Authors => {{Name => "Federico Galetto",
     	       Email => "galetto.federico@gmail.com",
	       HomePage => "https://math.galetto.org"}},
     Headline => "Buchsbaum-Eisenbud multipliers of free resolutions",
     PackageImports => {"SimpleDoc"}
     )


export {
    "buchsbaumEisenbudMultipliers",--method
    "bem",--shortcut
    "BEmults",--CacheTable key
    "exteriorDuality",--method
    "dualMultiplier"--method
    }


-------------------------------------------------------------------
-------------------------------------------------------------------
-- Exported methods
-------------------------------------------------------------------
-------------------------------------------------------------------


-- the following method computes Buchsbaum-Eisenbud
-- multipliers of a free resolution based on B-E's paper
-- "Some structure theorems for finite free resolutions"
-- The only difference is that in order to make each multiplier
-- homogeneous we give its domain the appropriate degree
buchsbaumEisenbudMultipliers = method()
bem = buchsbaumEisenbudMultipliers


-- WARNING: currently no safety checks are implemented!

-- this computes up to the k-th multiplier and returns it
buchsbaumEisenbudMultipliers(ChainComplex,ZZ) := Matrix => (F,k) -> (
    -- check if stored and return
    if F.cache#?BEmults then (
	if F.cache#BEmults#?k then return F.cache#BEmults#k;
	);
    -- otherwise create the hash tables to store results in cache
    F.cache#BEmults = new MutableHashTable;
    -- now we start the actual computation
    -- length of the resolution
    n := length F;
    -- ranks of the differentials
    r := apply(1..n,i->rank F.dd_i);
    -- format of the resolution
    f := apply(n+1,i->rank F_i);
    -- the a will house the next multiplier
    -- the last multiplier is simply an exterior power
    a := exteriorPower(r_(n-1),F.dd_n);
    F.cache#BEmults#n = a;
    -- i is a running index
    i := n-1;
    while (i >= k) do (
	-- the G's are rank one modules generated in the
	-- appropriate degree to make everything homogeneous.
	-- in fact this G is the domain of the next multiplier a
	G := exteriorPower(f_i,F_i);
	-- this code doesn't use the diagram of B-E's paper
	-- instead it uses its dual and accounts for degrees.
	-- e's are the exterior powers of the differentials
	e := (dual exteriorPower(r_(i-1),F.dd_i))**G;
	-- next: change of basis using exterior duality
	w := exteriorDuality(F,i);
	c := (dual w)**G;
	-- get next multiplier by factoring as in dual diagram
	b := e // (c*a);
	-- now dualize back and fix the degrees
	a = dual (b ** (dual G));
	F.cache#BEmults#i = a;
	i = i-1;
	);
    return F.cache#BEmults#k;
    )

-- this returns all multipliers in a list
buchsbaumEisenbudMultipliers(ChainComplex) := List => F -> (
    n := length F;
    -- we create a zero multiplier and prepend it to the list
    -- of multipliers created with the earlier function
    -- this is purely for convenience of indices
    R := ring F;
    {map(R^0,R^0,0)} | toList apply(1..n,k->bem(F,k))
    )


-- we need the iso of free modules Wedge^i F->Wedge^j F^*
-- from the pairing Wedge^i F ** Wedge^j F->Wedge^(i+j) F
-- where F is free has rank i+j.
-- Because of the way M2 lists subsets, this matrix
-- has +/-1 on the antidiagonal, 0's elsewhere
-- where the sign depends on the number of inversions
-- of the sets indexing row and column
exteriorDuality = method(TypicalValue => Matrix)

-- this computes the duality as a matrix of free abelian groups
exteriorDuality(ZZ,ZZ) := (r,n) -> (
    U := subsets(toList(0..n-1),r);
    V := subsets(toList(0..n-1),n-r);
    b := binomial(n,r);
    M := mutableMatrix map(ZZ^b,ZZ^b,0);
    for i to b-1 do (
	M_(b-i-1,i) = (-1)^(numberOfInversions(U_i,V_(b-i-1)));
	);
    return matrix M;
    )


-- this gives the duality for the k-th module in a free resolution
exteriorDuality(ChainComplex,ZZ) := (F,k) -> (
    r := rank F.dd_k;
    n := rank F_k;
    -- need a rank one free module to approriately twist degrees
    G := exteriorPower(n,F_k);
    domain := exteriorPower(r,F_k);
    codomain := exteriorPower(n-r,dual F_k) ** G;
    M := promote(exteriorDuality(r,n),ring F);
    map(codomain,domain,M)
    )


dualMultiplier = method(TypicalValue => Matrix)

-- returns the dual of a BE multiplier with the appropriate
-- degrees of domain and codomain (via a twist)
dualMultiplier(ChainComplex,ZZ) := (F,k) -> (
    G := exteriorPower(rank F_(k-1),F_(k-1));
    (dual bem(F,k)) ** G
    )


-------------------------------------------------------------------
-------------------------------------------------------------------
-- Unexported functions
-------------------------------------------------------------------
-------------------------------------------------------------------


-- u and v are complementary sublists of 0..n-1
-- this function counts the number of inversion to unshuffle u|v
numberOfInversions = (u,v) -> (
    c := 0;
    for j in v do (
	for i in u do (
	    if j-i<0 then c=c+1;
	    );
	);
    return c;
    )

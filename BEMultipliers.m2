newPackage(
     "BEMultipliers",
     Version => "0.2",
     Date => "August 24, 2019",
     AuxiliaryFiles => false,
     Authors => {{Name => "Federico Galetto",
     	       Email => "galetto.federico@gmail.com",
	       HomePage => "https://math.galetto.org"}},
     HomePage => "https://github.com/galettof/BEMultipliers",
     Headline => "Buchsbaum-Eisenbud multipliers of free resolutions",
     PackageImports => {"SimpleDoc"}
     )


export {
    "buchsbaumEisenbudMultipliers",--method
    "bem",--shortcut
    "BEmults",--CacheTable key
    --"exteriorDuality",--method
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
bem = buchsbaumEisenbudMultipliers = method()
--bem = buchsbaumEisenbudMultipliers


-- WARNING: currently no safety checks are implemented!

-- this computes up to the k-th multiplier and returns it
buchsbaumEisenbudMultipliers(ZZ,ChainComplex) := Matrix => (k,F) -> (
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
	w := exteriorDuality(i,F);
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
    {map(R^0,R^0,0)} | toList apply(1..n,k->bem(k,F))
    )


dualMultiplier = method(TypicalValue => Matrix)

-- returns the dual of a BE multiplier with the appropriate
-- degrees of domain and codomain (via twist by a rank one module)
dualMultiplier(ZZ,ChainComplex) := (k,F) -> (
    G := exteriorPower(rank F_(k-1),F_(k-1));
    ((dual bem(k,F)) ** G) * exteriorDuality(k-1,F)
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
exteriorDuality(ZZ,ChainComplex) := (k,F) -> (
    r := rank F.dd_k;
    n := rank F_k;
    -- need a rank one free module to approriately twist degrees
    G := exteriorPower(n,F_k);
    domain := exteriorPower(r,F_k);
    codomain := exteriorPower(n-r,dual F_k) ** G;
    M := promote(exteriorDuality(r,n),ring F);
    map(codomain,domain,M)
    )


-------------------------------------------------------------------
-------------------------------------------------------------------
-- Documentation
-------------------------------------------------------------------
-------------------------------------------------------------------


beginDocumentation()
doc ///
    Key
    	BEMultipliers
    Headline
    	Buchsbaum-Eisenbud multipliers of free resolutions
    Description
    	Text
	    @EM "BEMultipliers"@ is a Macaulay2 package that
	    can be used to compute Buchsbaum-Eisenbud multipliers
	    of minimal free resolutions of modules over polynomial
	    rings as introduced in
	    @HREF("https://doi.org/10.1016/S0001-8708(74)80019-8",
	    "Buchsbaum, Eisenbud - Some structure theorems for
	    finite free resolutions")@.
///

doc ///
    Key
    	BEmults
    Headline
    	cache key for Buchsbaum-Eisenbud multipliers
    Description
    	Text
	    When @TO "buchsbaumEisenbudMultipliers"@ is called to
	    compute Buchsbaum-Eisenbud multipliers of a resolution,
	    it stores the results of the computation in a
	    @TO "HashTable"@ in the cache of the resolution.
	    The key to access this information is @TT "BEmults"@.
///

doc ///
     Key
     	  buchsbaumEisenbudMultipliers
     Headline
     	  compute Buchsbaum-Eisenbud multipliers of a resolution
     Description
     	  Text
	       Use this method to compute Buchsbaum-Eisenbud
	       multipliers of a free resolution over a polynomial
	       ring.
///

doc ///
     Key
     	  "bem"
     Headline
     	  compute Buchsbaum-Eisenbud multipliers of a resolution
     Description
     	  Text
	       Use @TT "bem"@ as a synonym for @TO "buchsbaumEisenbudMultipliers"@.
///

doc ///
     Key
     	  (buchsbaumEisenbudMultipliers,ChainComplex)
     Headline
     	  compute Buchsbaum-Eisenbud multipliers of a resolution
     Usage
     	  a = buchsbaumEisenbudMultipliers(F)
     Inputs
     	  F:ChainComplex
     Outputs
     	  a:List
     Description
     	  Text
	       Use this method to compute all Buchsbaum-Eisenbud
	       multipliers of a free resolution @TT "F"@ over a
	       polynomial ring. By default, Macaulay2 does not
	       check that @TT "F"@ is actually a resolution.
	       
	       The output is a list @TT "a"@ containing all the
	       multipliers in increasing order. In order to match
	       the indexing of the paper by Buchsbaum and Eisenbud,
	       a zero map is inserted at the beginning of the list.
	       In other words, @TT "a_i"@ is the i-th multiplier.
	       
	       The results of the computation are stored in the
	       cache of @TT "F"@ with the key @TO "BEmults"@.
     	  Example
	       R=QQ[x,y,z]
	       K=koszul vars R
	       a=buchsbaumEisenbudMultipliers(K)
	       peek K.cache#BEmults
     SeeAlso
     	  BEmults
///

doc ///
     Key
     	  (buchsbaumEisenbudMultipliers,ZZ,ChainComplex)
     Headline
     	  return the k-th Buchsbaum-Eisenbud multiplier of a free resolution
     Usage
     	  a = buchsbaumEisenbudMultipliers(k,F)
     Inputs
     	  k:ZZ
     	  F:ChainComplex
     Outputs
     	  M:Matrix
     Description
     	  Text
	       Use this method to return a single Buchsbaum-Eisenbud
	       multiplier of a free resolution @TT "F"@ over a
	       polynomial ring. By default, Macaulay2 does not
	       check that @TT "F"@ is actually a resolution.
	       Note that the definition of the multipliers is
	       recursive, so all the previous ones are computed as
	       well (and stored).
	       
	       The output is the matrix of the k-th multiplier.
     	  Example
	       R=QQ[x,y,z]
	       K=koszul vars R
	       buchsbaumEisenbudMultipliers(2,K)
	       peek K.cache#BEmults
     SeeAlso
     	  BEmults
///

doc ///
     Key
     	  dualMultiplier
	  (dualMultiplier,ZZ,ChainComplex)
     Headline
     	  compute the dual of a Buchsbaum-Eisenbud multiplier
     Usage
     	  M = dualMultiplier(k,F)
     Inputs
     	  k:ZZ
     	  F:ChainComplex
     Outputs
     	  M:Matrix
     Description
     	  Text
	       Use this method to compute the dual of a
	       Buchsbaum-Eisenbud multipliers of a free
	       resolution over a polynomial ring. The dual of a
	       multiplier is the dual as a module homomorphism but
	       shifted in the appropriate degree. This is necessary
	       in order to make the commutative diagram of the
	       First Structure Theorem of Buchsbaum and Eisenbud
	       a diagram of graded modules.
     	  Example
	       R=QQ[x,y,z]
	       K=koszul vars R
	       a=buchsbaumEisenbudMultipliers(K)
     	  Text
	       We can use the dual of a multiplier to check the
	       First Structure Theorem holds.
     	  Example
	       exteriorPower(rank K.dd_2,K.dd_2) == bem(2,K) * dualMultiplier(3,K)
	       exteriorPower(rank K.dd_1,K.dd_1) == bem(1,K) * dualMultiplier(2,K)
	       
///

end


uninstallPackage "BEMultipliers"
restart
installPackage "BEMultipliers"
viewHelp

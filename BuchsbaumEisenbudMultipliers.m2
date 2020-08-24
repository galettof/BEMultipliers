newPackage(
     "BuchsbaumEisenbudMultipliers",
     Version => "0.2",
     Date => "August 24, 2019",
     AuxiliaryFiles => false,
     Authors => {{Name => "Federico Galetto",
     	       Email => "galetto.federico@gmail.com",
	       HomePage => "https://math.galetto.org"},
	   {Name => "Keller VandeBogert",
     	       Email => "kellerlv@math.sc.edu",
	       HomePage => "https://sites.google.com/view/kellervandebogert/home?authuser=0"}},
     HomePage => "https://github.com/galettof/BEMultipliers",
     Headline => "Buchsbaum-Eisenbud multipliers of free resolutions",
     PackageImports => {"SimpleDoc"}
     )


export {
    "aMultiplier",--method
    "cMultiplier",--method
    "ComputeRanks",--option
    "exteriorDuality"--method
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
aMultiplier = method(Options => {ComputeRanks => false})


-- WARNING: currently no safety checks are implemented!

-- this computes up to the k-th multiplier and returns it
aMultiplier(ZZ,ChainComplex) := Matrix => op -> (k,F) -> (
    -- get the index of the last nonzero module of the resolution
    n := max F;
    while F_n == 0 do n = n-1;
    if (k > n or k <= min F) then
    	error ("the first argument should be the homological
	    degree of the domain of a nonzero differential");
    -- check if stored and return
    if F.cache#?aMultiplier then (
	if F.cache#aMultiplier#?k then return F.cache#aMultiplier#k;
	)
    -- otherwise create the hash tables to store results in cache
    else (
	F.cache#aMultiplier = new MutableHashTable;
	);
    -- now we start the actual computation
    -- the a will house the next multiplier
    -- the last multiplier is simply an exterior power
    if op.ComputeRanks then (
	r := rank F.dd_n;
	) else (
	r = rank(n,F);
	);
    a := exteriorPower(r,F.dd_n);
    F.cache#aMultiplier#n = a;
    -- i is a running index
    i := n-1;
    while (i >= k) do (
    	if op.ComputeRanks then (
	    r = rank F.dd_i;
	    ) else (
	    r = rank(i,F);
	    );
	-- this code doesn't use the diagram of B-E's paper
	-- instead it uses its dual and accounts for degrees.
	-- e's are the exterior powers of the differentials
	e := dual exteriorPower(r,F.dd_i);
	-- G is a rank one module generated in the
	-- appropriate degree to make everything homogeneous.
	w := dual exteriorDuality(r,F_i);
	G := dual exteriorPower(rank F_i,F_i);
	-- get next multiplier by factoring as in dual diagram
	b := e // (w*(a**G));
	a = dual b;
	F.cache#aMultiplier#i = a;
	i = i-1;
	);
    return F.cache#aMultiplier#k;
    )

-- this computes the map a_j^k for lower order multipliers
aMultiplier(ZZ,ZZ,ChainComplex) := Matrix => op -> (j,k,F) -> (
    -- check if stored and return
    if F.cache#?aMultiplier then (
	if F.cache#aMultiplier#?(j,k) then return F.cache#aMultiplier#(j,k);
	);
    -- if F.cache#aMultiplier does not exist, it is created
    -- when aMultiplier is called
    -- EXT MULTIPLICATION IS MISSING
    m := wedgeProduct(rank(k,F),j,F_(k-1)) * (aMultiplier(k,F) ** id_(exteriorPower(j,F_(k-1))));
    return F.cache#aMultiplier#(j,k) = m;
    )

-- Weyman's lower order multipliers
cMultiplier = method(Options => {ComputeRanks => false})

cMultiplier(ZZ,ZZ,ChainComplex) := Matrix => op -> (j,k,F) -> (
    if not (j-1)*(length F)<=j*(k-1)-2 then (
	error"Not defined for chosen parameters";
	)
    else (
    	-- check if stored and return
    	if F.cache#?cMultiplier then (
	    if F.cache#cMultiplier#?(j,k) then return F.cache#cMultiplier#(j,k);
	    )
    	else (
	    F.cache#cMultiplier = new MutableHashTable;
	    );
    	if op.ComputeRanks then (
	    r := rank F.dd_(k-1);
	    ) else (
	    r = rank(k-1,F);
	    );
	-- see aMultiplier for comments on this code
	e := dual exteriorPower(r-j,F.dd_(k-1));
	w := dual exteriorDuality(r-j,F_(k-1));
	G := dual exteriorPower(rank F_(k-1),F_(k-1));
	b := e // (w*(aMultiplier(j,k,F)**G));
	return F.cache#cMultiplier#(j,k) = dual b;
	);
    )


-- Below is the iso of free modules Wedge^i F->Wedge^j F^* ** Wedge^(i+j) F
-- from the pairing Wedge^i F ** Wedge^j F->Wedge^(i+j) F
-- where F is free of rank i+j.
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

-- this is the duality Wedge^r F ->(Wedge^(n-r) F)^* ** Wedge^n F
-- for a free R-module F of rank n
exteriorDuality(ZZ,Module) := (r,F) -> (
    if isFreeModule(F) then (
    	n := rank F;
    	-- need a rank one free module to approriately twist degrees
    	G := exteriorPower(n,F);
    	domain := exteriorPower(r,F);
    	codomain := exteriorPower(n-r,dual F) ** G;
    	M := promote(exteriorDuality(r,n),ring F);
    	return map(codomain,domain,M);
    	) else (
	    error "Only implemented for free modules";
    	);
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


-- compute the ranks of the differentials in a free resolution
-- using the Buchsbaum-Eisenbud exactness criterion
-- F is the resolution, k is for the k-th differential
rank(ZZ,ChainComplex) := (k,F) -> (
    n := max F;
    r := rank F_n;
    while n>k do (
	n = n-1;
	r = rank F_n - r;
	);
    return r;
    )


-------------------------------------------------------------------
-------------------------------------------------------------------
-- Documentation
-------------------------------------------------------------------
-------------------------------------------------------------------


beginDocumentation()
doc ///
    Key
    	BuchsbaumEisenbudMultipliers
    Headline
    	Buchsbaum-Eisenbud multipliers of free resolutions
    Description
    	Text
	    @EM "BuchsbaumEisenbudMultipliers"@ is a Macaulay2 package that
	    can be used to compute Buchsbaum-Eisenbud multipliers
	    of minimal free resolutions of modules over polynomial
	    rings as introduced in
	    @HREF("https://doi.org/10.1016/S0001-8708(74)80019-8",
	    "Buchsbaum, Eisenbud - Some structure theorems for
	    finite free resolutions")@.
///

doc ///
     Key
     	  aMultiplier
     	  (aMultiplier,ZZ,ChainComplex)
     Headline
     	  compute a Buchsbaum-Eisenbud multiplier
     Usage
     	  aMultiplier(k,F)
     Inputs
     	  k:ZZ
     	  F:ChainComplex
     Outputs
     	  :Matrix
     Description
     	  Text
	       Use this method to compute a Buchsbaum-Eisenbud
	       multiplier of a free resolution @TT "F"@ over a
	       polynomial ring. By default, Macaulay2 does not
	       check that @TT "F"@ is actually a resolution.
	       The output is the matrix of the k-th multiplier.
     	  Example
	       R=QQ[x,y,z]
	       K=koszul vars R
	       aMultiplier(2,K)
     	  Text
	       We can check this multiplier satisfies the First
	       Structure Theorem of Buchsbaum and Eisenbud.
	       Note the exterior duality isomorphism must be made
	       explicit in order for the equality to hold.
	       The module @TT "E"@ is a rank one graded free
	       module that twists the map into the right degree.
     	  Example
	       E=exteriorPower(rank K_2,K_2)
	       exteriorPower(rank K.dd_2,K.dd_2) == aMultiplier(2,K) * ((dual aMultiplier(3,K))**E) * exteriorDuality(rank K.dd_2,K_2)
     	  Text
	       Note that the definition of the multipliers is
	       recursive, so all the previous ones are computed as
	       well (and stored).
     	  Example
	       peek K.cache#aMultiplier
///

doc ///
     Key
     	  (aMultiplier,ZZ,ZZ,ChainComplex)
     Headline
     	  compute a Buchsbaum-Eisenbud multiplier
     Usage
     	  aMultiplier(j,k,F)
     Inputs
     	  j:ZZ
     	  k:ZZ
     	  F:ChainComplex
     Outputs
     	  :Matrix
     Description
     	  Text
	       Use this method to compute a Buchsbaum-Eisenbud
	       multiplier of a free resolution @TT "F"@ over a
	       polynomial ring. By default, Macaulay2 does not
	       check that @TT "F"@ is actually a resolution.
	       The output is the matrix of the k-th multiplier.
     	  Example
	       R=QQ[x_1..x_4]
	       I=ideal apply(subsets(gens R,2),product)
	       RI=res I
	       aMultiplier(2,RI)
     	  Text
	       We can check this multiplier satisfies the First
	       Structure Theorem of Buchsbaum and Eisenbud.
	       Note the exterior duality isomorphism must be made
	       explicit in order for the equality to hold.
	       The module @TT "E"@ is a rank one graded free
	       module that twists the map into the right degree.
     	  Example
	       E=exteriorPower(rank RI_2,RI_2)
	       exteriorPower(rank RI.dd_2,RI.dd_2) == aMultiplier(2,RI) * ((dual aMultiplier(3,RI))**E) * exteriorDuality(rank RI.dd_2,RI_2)
     	  Text
	       Note that the definition of the multipliers is
	       recursive, so all the previous ones are computed as
	       well (and stored).
     	  Example
	       peek RI.cache#aMultiplier
///

doc ///
    Key
    	ComputeRanks
	[aMultiplier,ComputeRanks]
    Headline
    	Explicitly compute ranks of differentials
    Description
    	Text
	    The methods @TO "aMultiplier"@ and @TT "cMultiplier"@
	    compute the ranks of the differentials in a resolution
	    using the Buchsbaum-Eisenbud exactness criterion.
	    Setting this optional argument to @TO "true"@ computes
	    ranks directly using the function @TO "rank"@.
///

doc ///
     Key
     	  exteriorDuality
     	  (exteriorDuality,ZZ,ZZ)
     	  (exteriorDuality,ZZ,Module)
     Headline
     	  exterior duality isomorphism of free modules
     Usage
     	  exteriorDuality(k,n)
     	  exteriorDuality(k,F)
     Inputs
     	  k:ZZ
     	  F:Module
     Outputs
     	  :Matrix
     Description
     	  Text
	       If @TT "F"@ is a free module of rank @TT "n"@,
	       this module constructs the isomorphism
	       $$\wedge^k F \to (\wedge^{n-k} F)^* \otimes \wedge^n F$$
	       induced by the perfect pairing
	       $$\wedge^k F \otimes \wedge^{n-k} F\to \wedge^n F.$$
	       The form @TT "exteriorDuality(k,n)"@ produces a
	       matrix over the integers, while the form
	       @TT "exteriorDuality(k,F)"@ produces a matrix over
	       the ring of the module @TT "F"@.
     	  Example
	       R=QQ[x,y,z]
	       F=R^5
	       exteriorDuality(2,5)
	       exteriorDuality(2,F)
///

end


uninstallPackage "BuchsbaumEisenbudMultipliers"
restart
installPackage "BuchsbaumEisenbudMultipliers"
installPackage("BuchsbaumEisenbudMultipliers",RemakeAllDocumentation=>true)

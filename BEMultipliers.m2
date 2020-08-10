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
    "aMultiplier",--method
    "exteriorDuality",--method
    "dualMultiplier",--method
    "lowerBEM"
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
aMultiplier = aMultiplier = method()


-- WARNING: currently no safety checks are implemented!

-- this computes up to the k-th multiplier and returns it
aMultiplier(ZZ,ChainComplex) := Matrix => (k,F) -> (
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
    a := exteriorPower(rank F.dd_n,F.dd_n);
    F.cache#aMultiplier#n = a;
    -- i is a running index
    i := n-1;
    while (i >= k) do (
	-- this code doesn't use the diagram of B-E's paper
	-- instead it uses its dual and accounts for degrees.
	-- e's are the exterior powers of the differentials
	e := dual exteriorPower(rank F.dd_i,F.dd_i);
	-- G is a rank one module generated in the
	-- appropriate degree to make everything homogeneous.
	w := dual exteriorDuality(rank F.dd_i,F_i);
	G := dual exteriorPower(rank F_i,F_i);
	-- get next multiplier by factoring as in dual diagram
	b := e // (w*(a**G));
	a = dual b;
	F.cache#aMultiplier#i = a;
	i = i-1;
	);
    return F.cache#aMultiplier#k;
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



-------------------------------------------------------------------
-------------------------------------------------------------------
-- Code below provided by Keller VandeBogert on 8/4/2020
-------------------------------------------------------------------
-------------------------------------------------------------------

lowerBEM = (j,k,F) -> (
    if not (j-1)*(length F)<=j*(k-1)-2 then (
	error"Not defined for chosen parameters";
	)
    else (
    -- if (j-1)*(length F)<=j*(k-1)-2 then (
	fRank:=rank F_(k-1);
	dRank1:=rank(F.dd_(k));
	dRank2:=rank(F.dd_(k-1));
	inda:=wedgeProduct(dRank1,j,F_(k-1))*((aMultiplier(k,F))**id_(exteriorPower(j,F_(k-1))));
	--wedgeD:= dual exteriorPower(dRank2-j,F.dd_(k-1));
	w:=exteriorPower(dRank2-j,F.dd_(k-1));
	extD:=exteriorDuality(dRank2-j,k-1,F);
	G := exteriorPower(rank F_(k-1),F_(k-1));
	return (((dual (w*extD))**G)//inda);
	--return dual((extD*(matrix entries wedgeD))//(matrix entries inda));
	);
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

end


uninstallPackage "BEMultipliers"
restart
installPackage "BEMultipliers"
viewHelp

doc ///
Key
  CanonicalComplexes
Headline
  Construct standard classical complexes with correct choice of bases and compute Buchsbaum-Eisenbud multipliers on complexes.
Description
  Text
   This package contains functions for computing certain well known complexes. Many of the ideals or modules
   resolved by these complexes can be easily resolved much more quickly by the "res" command, but
   Macaulay2 often does not choose bases in a fashion that one might prefer.
   
   By using complexes with coherent choice of bases, one can run tests and deduce natural choices for operations
   on chain complexes involved in more theoretical computations. Common examples may be for the computation
   of comparison maps/morphisms of complexes or DG-algebra structures.
   
   Buchsbaum-Eisenbud multipliers of all orders have also been implemented as in the structure theorem
   proved by Buchsbaum and Eisenbud and later generalized by Weyman for lower order minors.
   
   Examples need to be added here.
SeeAlso
  "ChainComplexExtras"
  "ChainComplexOperations"
  "DGAlgebras"
  "SchurComplexes"
///

doc ///
Key
  genList
  (genList, ZZ, ZZ, Module)
Headline
    Outputs the basis order used for the Schur modules appearing in the corresponding L-complexes.
Usage
    G = genList(i,j,F)
Inputs
    i,j: ZZ
    	The indices on the corresponding tensor product of symmetric and exterior powers that
	the Schur module injects into.
    F: Module
    	A free module.
Outputs
    G: MatrixExpression
    	A matrix expression listing the basis order that is used for the Schur modules appearing in 
	the corresponding L-complexes.
Description
    Text
    	Let $F = \mathbb{Q}^3$ with standard basis denoted $f_1, \ f_2, \ f_3$. We may compute the basis 
	of $L_2^1 (\mathbb{Q}^3)$ and $L_2^2 (\mathbb{Q}^3)$ with the ordering that Macaulay2 is using as follows:
    Example
    	genList(1,2,QQ^3)
	genList(2,2,QQ^3)
    Text
    	The above tells us that the first generator of $L_2^1 (\mathbb{Q}^3)$ is given by
	$$f_2 \otimes f_1^2 - f_1 \otimes f_1 f_2$$
	and the third generator of $L_2^2 (\mathbb{Q}^3)$ is given by
	$$f_2 \wedge f_3 \otimes f_1 f_3 - f_1 \wedge f_3 \otimes f_2 f_3 + f_1 \wedge f_2 \otimes f_3^2.$$
SeeAlso
    Lcomplex
///


doc ///
Key
  Lcomplex
  (Lcomplex, ZZ, Matrix)
  (Lcomplex, ZZ, Ideal)
  (Lcomplex, ZZ, Ring)
Headline
    Outputs the L-complex corresponding to a map $\phi : F \to R$, where $F$ is a free module.
Usage
    F=Lcomplex(b,phi)
Inputs
    b: ZZ
    	This means that all terms appearing in the complex are of the form $L_b^i$ for $i$ ranging
	over the appropriate parameters.
    phi: Matrix
    	A matrix $\phi$ whose target is the rank one free module $R$, where $R$ is the ambient ring.
Outputs
    F: ChainComplex
    	The L-complexes of Buchsbaum and Eisenbud. The $0$-th homology will be $R/ \im (\phi)^b$, and 
	if $\im (\phi)$ has large enough grade, this will be a minimal free resolution.
Description
    Text
    	Let $R = k[x_1 , x_2 , x_3]$ and $\phi : Re_1 \oplus Re_2 \opls Re_3 \to R$ be the map induced
	by sending $e_i \mapsto x_i$. The image has grade $3$, so the L-complex will be acyclic for all
	values of $b$.
    Example
    	R=QQ[x_1..x_3];
	F=Lcomplex(3,vars R)
	prune HH F
    Text
	As expected, the complex is acyclic with $H_0 (F) = R/(x_1,x_2,x_3)^3$. The command genList may
	be used to see the basis being used in each homological degree of F.
SeeAlso
    genList
    Kcomplex
///

doc ///
Key
  Kcomplex
  (Kcomplex, ZZ, Matrix)
  (Kcomplex, ZZ, Ideal)
  (Kcomplex, ZZ, Ring)
Headline
    Outputs the K-complex corresponding to a map $\phi : F \to R$, where $F$ is a free module.
Usage
    F=Kcomplex(b,phi)
Inputs
    b: ZZ
    	This means that all terms appearing in the complex are of the form $K_b^i$ for $i$ ranging
	over the appropriate parameters.
    phi: Matrix
    	A matrix $\phi$ whose target is the rank one free module $R$, where $R$ is the ambient ring.
Outputs
    F: ChainComplex
    	The K-complexes of Buchsbaum and Eisenbud. The $0$-th homology will be isomorphic to the canonical module
	of $R/ \im(\phi)^{b+1}$ (up to a grading shift) if $\im (\phi)$ has large enough grade, and will be a minimal free resolution.
Description
    Text
    	Let $R = k[x_1 , x_2 , x_3]$ and $\phi : Re_1 \oplus Re_2 \opls Re_3 \to R$ be the map induced
	by sending $e_i \mapsto x_i$. The image has grade $3$, so the K-complex will be acyclic for all
	values of $b$.
    Example
    	R=QQ[x_1..x_3];
	F=Kcomplex(2,vars R)
	prune HH F
    Text
	As expected, the complex is acyclic. It remains to see that the $0$th homology is correct. We can easily
	compute the canonical module of $R/(x_1,x_2,x_3)^3$:
    Example
    	prune Ext^3 (R^1/(ideal vars R)^3,R^1)
    Text
    	It is easy to see that this is isomorphic to $H_0 (F)$ up to a shift. It is also well known that 
	the L and K-complexes are dual to each other. We verify that this is true in the following example:
    Example
    	G=Lcomplex(3,vars R);
	E=extend(G,(dual F)[-3],id_(G_0))
SeeAlso
    genList
    Lcomplex
///

doc ///
Key
    EN
    (EN, ZZ, Matrix)
Headline
    The Eagon-Northcott family of complexes associated to a matrix $M$.
Usage
    F=EN(i,M)
Inputs
    i: ZZ
    	An integer specifying the parameter arising in the associated Eagon-Northcott
	family of complexes.
    M: Matrix
    	A matrix whose entries the differentials of the resulting complex will depend on.
Outputs
    F: ChainComplex
    	The associated complex in the Eagon-Northcott family of complexes.
Description
    Text
    The Eagon-Northcott family of complexes is described explicitly in Eisenbud's book 
    "Commutative Algebra with a View Toward Algebraic Geometry" on page 595. Let $M$ denote an
    $n \times m$ matrix with $n \leq m$ and assume that $I_n (M)$ (the ideal of $n \times n$ minors of $M$)
    has grade $m-n+1$ (which is the case for $M$ generic). Then, the complexes arising in the
    Eagon-Northcott family will be acyclic for $-1 \leq i$. The special cases $\texttt{EN} (0,M)$ and 
    $\texttt{EN}(1,M)$ are the classical Eagon-Northcott and Buchsbaum-Rim complexes, respectively.
    Example
    	M=genMat(2,4)
	EN(0,M)
	prune HH oo
	EN(1,M)
	prune HH_0 oo == prune coker M
    Text
    	The Eagon-Northcott family of complexes is "closed under taking duals". More precisely, it is always true
	that $\texttt{EN} (i,M)^* = \texttt{EN} (m-n-i,M)$. We can verify this directly for a single
	example:
    Example
    	extend(EN(0,M),(dual EN(2,M))[-3],id_(Q^1))
    Text
        This also allows us to produce families of self-dual complexes when $m-n$ is an even integer:
    Example
    	EN(2,genMat(2,6))
	EN(3,genMat(3,9))
	EN(1,genMat(3,5))
	oo.dd
SeeAlso
    genEN
    genENBlist
    Ccomplex
///

doc ///
Key
    genEN
    (genEN, ZZ, List)
Headline
    The generalized Eagon-Northcott complex associated to a determinantal facet ideal.
Usage
    F=genEN(i,P)
Inputs 
    i: ZZ
        An integer, specifying the size of the minors you wish to take.
    P: List
    	A list of lists, where entry lists the corresponding maximal cliques of the simplicial
	complex.
Outputs
    F: ChainComplex
    	The generalized Eagon-Northcott complex associated to the determinantal facet ideal.
Description
    Text
    	The generalized Eagon-Northcott complex was introduced by Herzog, Kiani, and Madani. It arises
	as the linear strand of certain classes of determinantal facet ideals associated to clutters. 
	It can be described in the following way: let $M$ denote a generic $n\times m$ matrix with $n\leq m$
	and $\Delta$ an $(n-1)$-pure simplicial complex. Then the generalized Eagon-Northcott complex
	is the subcomplex of the classical Eagon-Northott complex formed by restricting to all basis elements
	of the form $f_\sigma \otimes g^{(\alpha)} \in \bigwedge^{n+\ell} F \otimes D_\ell (G^*)$, 
	where $\sigma \in Delta$ and $M$ is viewed as a morphism of free modules $M : F \to G$.
    Example
    	F=genEN(2,{{1,2,3,4},{3,4,5}})
    	F.dd
	prune HH F
    Text
    	This can also recover the classical Eagon-Northcott complex on a generic matrix as so:
    Example
    	genEN(2,{{1,2,3,4}})
    Text
    	The generalized Eagon-Northcott complex can be used to detect the existence of minimal nonfaces
	in a simplicial complex. More precisely, it is known that the generalized Eagon-Northcott complex
	associated to a simplicial complex is acyclic in linear degrees if and only if the associated
	clique complex of $\Delta$ has no minimal nonfaces.
SeeAlso
    genENBlist
    EN
    Ccomplex
///


doc ///
Key
    genENBlist
    (genENBlist, ZZ, ZZ, Matrix)
Headline
    Display the basis list used for the generalized Eagon-Northcott complex.
Usage
    L=genENBlist(i,j,P,G)
Inputs
    i,j: ZZ
    	Integers specifying that the basis list will be viewed as a submodule of $\bigwedge^i F \otimes D_b (G^*)$
	if $M : F \to G$. 
    P: List
    	A list giving the maximal cliques of the simplicial complex associated to the determinantal
	facet ideal.
    G: Module
        The target free module, whose rank is the size of the minors of the matrix.
Outputs
    L: List
    	The list giving the exact order of basis elements used in the genEN command. 
Description
    Text
    	Consider a $1$-pure simplicial complex with associated clique complex $\{ \{ 1,2,3 \}, \{2,3,4 \} \}$. 
	Then the basis elements/order used in the associated generalized Eagon-Northcott complex
	can be computed as follows:
    Example
    	F=genEN(2,{{1,2,3},{2,3,4}})
    	for i from 1 to 2 do print genENBlist(i+1,i-1,{{1,2,3},{2,3,4}},QQ^2)
    Text
    	This means that in homological degree $2$, the first basis element of $F$ is
	$f_1 \wedge f_2 \wedge f_3 \otimes g_1^*$. 
SeeAlso
    genEN
    EN
    Ccomplex
///

doc ///
Key
    Ccomplex
    (Ccomplex, ZZ, ZZ, Matrix)
Headline
    The family of complexes constructed by Kustin in the paper "Canonical Complexes Associated to a Matrix".
Usage
    F=Ccomplex(i,a,phi)
Inputs
    i,a: ZZ
    	Integers specifying the parameters used in the family of complexes.
    phi: Matrix
    	An $n \times m$ matrix with $m \leq n$.
Description
    Text
    	The C-complexes were constructed by Kustin and can be viewed as a generalization of
	the Eagon-Northcott complexes. In particular, the command $\texttt{Ccomplex} (i,1,M)$ will
	be isomorphic to $\texttt{EN} (i,M^*)$. We can see this in the following example:
    Example
    	F1=Ccomplex(0,1,M=genMat(4,2))
	F2=EN(0,transpose M)
	extend(F1,F2,id_(F1_0))
    Text
    	The C-complexes are more general than the Eagon-Northcott family, however. The following
	is example 7.7 in the paper by Kustin:
    Example
    	C=Ccomplex(0,2,genMat(4,3))
	C.dd
    Text
    	Since both differentials are quadratic, this complex does not arise from any complex in
	the Eagon-Northcott family.
SeeAlso
    EN
    genEN
    genENBlist
///

doc ///
Key
    tateLike
    (tateLike, ChainComplex,List,ZZ)
    (tateLike, ChainComplex)
Headline
    Construct a Tate-like complex associated to a complex that is a cyclic module in homological
    degree 0.
Usage
    F=tateLike(C,L,i)
Inputs
    C: ChainComplex
    	A chain complex with $C_0 = R$ a rank $1$ free module.
    L: List
    	A list of length equal to the length of C, specifying the topmost "exponent vector"
        to appear in the Tate-like complex.
    i: ZZ
    	An integer specifying which homological degree you would like the topmost term to appear
        in. If i is not provided, this term is placed in the "natural" homological degree, determined
	by L.
Outputs
    F: ChainComplex
    	A Tate-like complex associated to the input data.
Description
    Text
    	Let $C : 0 \to C_n \to  \cdots \to C_1 \to C_0 = R$ denote a complex of length $n$. Then, command
	$\texttt{tateLike} (F,\{ a_1 , \dotsc , a_n \} , i)$ creates a complex whose term in homological
	degree $i$ is $\bigwedge^{a_1} F_1 \otimes D_{a_2} (F_2) \otimes \cdots$, where exterior or
	divided powers are taken depending on whether the index is odd or even, respecively. The differential
	in this complex is induced by the formal product rule, sometimes referred to as the "Tate differential".
	For example, if $C = C_3 \to C_2 \to C_1 \to R$ is a complex, then $\texttt{tateLike} (C, \{1,2,0 \} , 3)$
	will output the complex
	$$C_1 \otimes D_2 (C_2) \to D_2(C_2) \oplus \bigwedge^2 C_1 \otimes C_2 \to C_1 \otimes C_2 \oplus \bigwedge^3 C_1.$$
	If $i$ had not been specified, then the topmost term would have been placed in homological
	degree $1 \cdot 1 + 2 \cdot 2 = 5$. In general, if a homological degree is not specified,
	the topmost term will be placed in homological degree $a_1 + 2a_2 + \cdots + n a_n$.
    Example
       Q=QQ[x_1..x_3]
       K=koszul vars Q;
       tateLike(K,{1,2,0},2)
       prune HH oo
   Text
       If $T$ denotes the Tate-like complex above, then one can verify that $H_1 (T)$ is isomorphic
       to the canonical module of $Q/(x_1,x_2,x_3)^3$. Tate-like complexes can also be used to compute
       DG algebra structures on complexes. If the length of the complex is more than $3$, then the 
       resulting product is not necessarily associative. However, in the length $\leq 3$ case, the
       product is guaranteed to be associative.
       
       In the following example, we compute an associative DG algebra structure on the L-complex
       $\texttt{Lcomplex} (2,Q)$.
    Example
    	Q=QQ[x_1..x_3]
	F=Lcomplex(2,Q);
	TL1 = tateLike(F,{2,0,0})
	(extend(F,TL1,id_(Q^1)))_2 --this gives the product between deg 1 elements
	genList(0,2,QQ^3)
	genList(1,2,QQ^3)
	genList(2,2,QQ^3)
    Text
    Identifying $L_b^a (Q^3)$ with the cokernel of the natural map $\bigwedge^{a+1} F \otimes S_{b-1} (F) \to \bigwedge^a F \otimes S_b (F)$,
    we find, for instance,
    $$(f_1 \otimes f_1) \cdot (f_1 \otimes f_2) = x_1 f_1 \wedge f_2 \otimes f_1,$$
    $$(f_1 \otimes f_1) \cdot (f_2 \otimes f_2) = x_2 f_1 \wedge f_2 \otimes f_1 + x_1 f_1 \wedge f_2 \otimes f_2.$$
    Likewise, for the product between elements of homological degree $1$ and $2$:
    Example
	TL2 = tateLike(F,{1,1,0},3)
	(extend(F,TL2,id_(Q^1)))_3
    Text
    	The above shows that
	$$(f_1 \otimes f_1) \cdot (f_1 \wedge f_2 \otimes f_3) = -x_1 f_1 \wedge f_2 \wedge f_3 \otimes f_1.$$
SeeAlso
    EN
    Ccomplex
///

doc ///
Key
    gulliksenNegard
    (gulliksenNegard, Matrix)
Headline
    Compute the classical Gulliksen-Negard complex.
Usage
    F=gulliksenNegard(M)
Inputs
    M: Matrix
    	M is an $n \times n$ matrix.
Outputs
    F: ChainComplex
    	The classical Gulliksen-Negard complex.
Description
    Text
    	Given an $n \times n$ matrix $M$ with $I_{n-1} (M)$ of grade $4$ (for instance, $M$ is generic), the Gulliksen-Negard
	complex yields a minimal free resolution of the quotient ring defined by $I_{n-1} (M)$.
	In particular, one finds that $I_{n-1} (M)$ is always a grade $4$ Gorenstein ideal. The
	construction of the Gulliksen-Negard complex implemented here comes from Bruns-Vetter's book
	"Determinantal Rings".
    Example
    	M=genMat(4,4)
	F=gulliksenNegard M
    	betti F
	betti gulliksenNegard genMat(5,5)
    Text
    	We can verify directly that this complex is self-dual:
    Example
    	F=gulliksenNegard genMat(3,3);
	F.dd
    	extend(F,(dual F)[-4],id_(F_0))	
SeeAlso
    EN
    Ccomplex
    tateLike
///

doc ///
Key
    extComult
    (extComult, ZZ, ZZ, Module)
Headline
    Compute comultiplication in the exterior algebra.
Usage
    M = extComult(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of comultiplication.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the comultiplication map $\bigwedge^{p+q} F \to \bigwedge^p F \otimes \bigwedge^q F$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extContract
    divComult
    divProduct
    symComult
    symProduct
    divContract
    revdivContract
///
	
doc ///
Key
    extContract
    (extContract, ZZ, ZZ, Module)
Headline
    Compute the contraction map in the exterior algebra.
Usage
    M = extContract(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of contraction.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the contraction map $\bigwedge^p F^* \otimes \bigwedge^q F \to \bigwedge^{q-p} F$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    divComult
    divProduct
    symComult
    symProduct
    divContract
    revdivContract
///

doc ///
Key
    divComult
    (divComult, ZZ, ZZ, Module)
Headline
    Compute the comultiplication map in the divided power algebra.
Usage
    M = divComult(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of comultiplication.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the comultiplication map $D_{p+q} (F) \to D_p (F) \otimes D_q (F)$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    extContract
    divProduct
    symComult
    symProduct
    divContract
    revdivContract
///

doc ///
Key
    divProduct
    (divProduct, ZZ, ZZ, Module)
Headline
    Compute the multiplication map in the divided power algebra.
Usage
    M = divProduct(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of multiplication.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the multiplication map $ D_p (F) \otimes D_q (F) \to D_{p+q} (F)$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    extContract
    divComult
    symComult
    symProduct
    divContract
    revdivContract
///


doc ///
Key
    symComult
    (symComult, ZZ, ZZ, Module)
Headline
    Compute the comultiplication map in the symmetric algebra.
Usage
    M = symComult(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of comultiplication.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the multiplication map $ S_{p+q} (F)\to S_p (F) \otimes S_q (F)$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    extContract
    divComult
    divProduct
    symProduct
    divContract
    revdivContract
///


doc ///
Key
    symProduct
    (symProduct, ZZ, ZZ, Module)
Headline
    Compute the multiplication map in the symmetric algebra.
Usage
    M = symProduct(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of multiplication.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the multiplication map $S_p (F) \otimes S_q (F)\to S_{p+q} (F)$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    extContract
    divComult
    divProduct
    symComult
    divContract
    revdivContract
///

doc ///
Key
    divContract
    (divContract, ZZ, ZZ, Module)
Headline
    Compute the contraction map on divided powers.
Usage
    M = divContract(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of the contraction map.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the contraction map $S_p (F^*) \otimes D_q (F) \to D_{q-p} (F)$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    extContract
    divComult
    divProduct
    symComult
    symProduct
    revdivContract
///

doc ///
Key
    revdivContract
    (revdivContract, ZZ, ZZ, Module)
Headline
    Compute the contraction map on divided powers.
Usage
    M = revdivContract(p,q,F)
Inputs
    p,q: ZZ
    	Integers specifying the parameters of the contraction map.
    F: Module
    	A free module.
Outputs
    M: Matrix
    	The matrix representation of the contraction map $D_q (F) \otimes S_p (F^*) \to D_{q-p} (F)$
	with respect to the default Macaulay2 basis order.
SeeAlso
    extComult
    extContract
    divComult
    divProduct
    symComult
    symProduct
    divContract
///


doc ///
Key
    genMat
    (genMat, ZZ, ZZ)
    (genMat, ZZ, ZZ, Ring)
Headline
    Defines the coordinate ring of a generic matrix and outputs the matrix.
Usage
    M=genMat(n,m,R)
Inputs
    n,m: ZZ
    	Integers specifying the size of the desired matrix.
    R: Ring
    	A ring specifying which coordinate ring the matrix will live in. If no ring $R$ is given, 
	the default ring will be $\mathbb{Q}$.
Outputs
    M: Matrix
    	A generic $n \times m$ matrix over the ring $Q = R [ x_{i,j} \ | \ 1 \leq i \leq n, \ 1 \leq j \leq m]$.
Description
    Example
    	M=genMat(2,4)
	EN(2,M)
	genMat(4,8)
SeeAlso
    EN
    Ccomplex
///
	

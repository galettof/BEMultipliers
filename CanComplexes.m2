needsPackage "KustinMiller"
needsPackage "ChainComplexOperations"
needsPackage "ChainComplexExtras"

--=======================================================================
--Multilinear algebra codes

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


--computes sign of permutation
sgn = v -> (c:=0;
for i in v do (
	for j in v_{0..position(v,l->l==i)} do (
		if i-j < 0 then c = c+1;
	);
);
return (-1)^c;
)



--this code computes the comultiplication map
--wedge^(p+q) F -> wedge^p F ** wedge^q F 
--with respect to the standard basis orders
extComult = (p,q,F) -> (L:={};
      n:=rank F;
      b1:=binomial(n,p)-1;
      b2:=binomial(n,q)-1;
      b3:=binomial(n,p+q)-1;
      L1:=subsets(toList(0..n-1),p);
      L2:=subsets(toList(0..n-1),q);
	--labels for basis of p+q exterior power:
      L3:=subsets(toList(0..n-1),p+q);
	--labels for the basis of the tensor product of 
	--p and q exterior powers:
      for i from 0 to b1 do (
		for j from 0 to b2 do (
		L=L|{(L1_i,L2_j)}; );
	);
      Mmut := mutableMatrix(ring(F),#L,#L3);
      for k from 0 to b3 do (
      for l from 0 to (#L)-1 do (
      	if (set(flatten(L_l))===set((L3)_k)) then (
		Mmut_(l,k) = (-1)^(numberOfInversions((L_l)_0,(L_l)_1)); 
	);
	);
	);
return map(exteriorPower(p,F)**exteriorPower(q,F),exteriorPower(p+q,F),promote(matrix Mmut, ring(F)));
)





--this function computes the contraction map
--sending wedge^r (F*) ** wedge^n (F) -> wedge^(n-r) F
extContract = (r,n,F) -> (D:=extComult(r,n-r,F);
M:=map((ring(F))^1,exteriorPower(r,dual F)**exteriorPower(r,F),matrix{flatten entries id_(exteriorPower(r,F))});
a:=id_(exteriorPower(r,dual(F)))**D;
b:=M**id_(exteriorPower(n-r,F));
C:=b*a;
return C;
)



--this function computes the comultiplication
--map on the divided power algebra 
--D_(p+q) F -> D_p F ** D_q F with respect
--to the GrRevLex basis order
divComult = (p,q,F) -> (n:=rank(F);
L:=compositions(n,p+q);
L1:=compositions(n,p);
L2:=compositions(n,q);
L3=new MutableList;
for i from 0 to (#L1)-1 do (
	for j from 0 to (#L2)-1 do (
		L3#(#L3)=(L1_i,L2_j);
	);
);
L3=toList L3;
Mmut=mutableMatrix(ring(F),#L3,#L);
for i from 0 to (#L)-1 do (
	for j from 0 to (#L3)-1 do (
		for k in L1 do (
			if (k,L_i-k)==L3_j then Mmut_(j,i)=1;
		);
	);
);
return map(symmetricPower(p,F)**symmetricPower(q,F),symmetricPower(p+q,F),matrix Mmut);
)


--this computes the natural map
--D_p F ** D_q F -> D_(p+q) F with
--respect to GrRevLex basis order
divProduct = (p,q,F) -> (n:=rank(F);
L:=compositions(n,p+q);
L1:=compositions(n,p);
L2:=compositions(n,q);
L3=new MutableList;
for i from 0 to (#L1)-1 do (
	for j from 0 to (#L2)-1 do (
		L3#(#L3)=(L1_i,L2_j);
	);
);
L3=toList L3;
Mmut:=mutableMatrix(ring(F),#L,#L3);
for i from 0 to (#L)-1 do (
	for j from 0 to (#L3)-1 do (
		for k in L1 do (
			if (k,L_i-k)==L3_j then Mmut_(i,j)=product(k,L_i-k,(l,m)->binomial(l+m,m));
		);
	);
);
return map(symmetricPower(p+q,F),symmetricPower(p,F)**symmetricPower(q,F),matrix Mmut);
)


--this computes the comultiplication
--map S_(p+q) F -> S_p F ** S_q F with
--respect to GrRevLex order.
symComult = (p,q,F) -> (n:=rank(F);
L:=compositions(n,p+q);
L1:=compositions(n,p);
L2:=compositions(n,q);
L3:={};
for i from 0 to (#L1)-1 do (
	for j from 0 to (#L2)-1 do (
		L3=L3|{(L1_i,L2_j)};
	);
);
Mmut:=mutableMatrix(ring(F),#L3,#L);
for i from 0 to (#L)-1 do (
	for j from 0 to (#L3)-1 do (
		for k in L1 do (
			if (k,L_i-k)==L3_j then Mmut_(j,i)=product(k,L_i-k,(l,m)->binomial(l+m,m));
		);
	);
);
return map(symmetricPower(p,F)**symmetricPower(q,F),symmetricPower(p+q,F),matrix Mmut);
)

--this computes the natural map
--S_p F ** S_q F -> S_(p+q) F with
--respect to GrRevLex basis order
symProduct = (p,q,F) -> (n:=rank(F);
L:=compositions(n,p+q);
L1:=compositions(n,p);
L2:=compositions(n,q);
L3:={};
for i from 0 to (#L1)-1 do (
	for j from 0 to (#L2)-1 do (
		L3=L3|{(L1_i,L2_j)};
	);
);
Mmut:=mutableMatrix(ring(F),#L,#L3);
for i from 0 to (#L)-1 do (
	for j from 0 to (#L3)-1 do (
		for k in L1 do (
			if (k,L_i-k)==L3_j then Mmut_(i,j)=1;
		);
	);
);
return map(symmetricPower(p+q,F),symmetricPower(p,F)**symmetricPower(q,F),matrix Mmut);
)

--this computes the contraction map
--S_r (F) ** D_n (F) --> D_(n-r) F.
divContract = (r,n,F) -> (D:=divComult(r,n-r,F);
M:=map((ring(F))^1,symmetricPower(r,dual F)**symmetricPower(r,F),matrix{flatten entries id_(symmetricPower(r,F))});
a:=id_(symmetricPower(r,dual(F)))**D;
b:=M**id_(symmetricPower(n-r,F));
C:=b*a;
return C;
)

--this computes the contraction map
--D_n(F) **S_r (F)  --> D_(n-r) F.
revdivContract = (n,r,F) -> (D:=divComult(n-r,r,F);
M:=map((ring(F))^1,symmetricPower(r,dual F)**symmetricPower(r,F),matrix{flatten entries id_(symmetricPower(r,F))});
a:=D**id_(symmetricPower(r,dual(F)));
b:=id_(symmetricPower(n-r,F))**M;
C:=b*a;
return C;
)

--====================================================================================================
-- Buchsbaum Eisenbud L and K complexes

--This computes the module L_b^a (F) (a schur module
--corresponding to the appropriate hook partition)
hook = (a,b,F) -> (C:=extComult(a-1,1,F);
M:=map(symmetricPower(b+1,F),symmetricPower(1,F)**symmetricPower(b,F),transpose divComult(1,b,F));
K:=map(exteriorPower(a-1,F)**symmetricPower(b+1,F),exteriorPower(a,F)**symmetricPower(b,F),(id_(exteriorPower(a-1,F))**M)*(C**id_(symmetricPower(b,F))));
L:= kernel K;
return L;
)

--lists the generators of the hook module
--in the order that M2 is using, I can explain what the output means
genList = (a,b,F) -> (L=matrix gens hook(a,b,F);
n = rank source L;
m=rank target L;
L1=flatten table(subsets(1..rank F,a),compositions(rank F,b),(u,v)->{u,v});
L2=new MutableList;
for i from 0 to n-1 do (
	for j from 0 to m-1 do (
		if not((L_{i})_(j,0)==0) then (L2#(#L2)=flatten {(L_{i})_(j,0),L1_j};);
	);
);
MatrixExpression pack(a+1,toList L2)
)

--this computes the map L_b^a(F) -> L_b^a-1 (F) induced by the map phi.
--phi here must be a map F -> R.
indMap= (a,b,phi) -> (F=source phi;
T:=(-1)^a*(koszul(a,phi)**id_(symmetricPower(b,F)))*inducedMap(exteriorPower(a,F)**symmetricPower(b,F),hook(a,b,F));
I:=(-1)^(a-1)*inducedMap(exteriorPower(a-1,F)**symmetricPower(b,F),hook(a-1,b,F));
N:=prune(T//I);
N)

--this gives the standard L-complex of Buchsbaum-Eisenbud. Here phi : F -> R;
--if grade (im phi) = rank F, this is a minimal free resolution of (im phi)^b.
Lcomplex = (b,phi) -> (F:=source phi;
L=new MutableList;
--L#1=matrix{rsort ((entries symmetricPower(b,phi))_0)};
for i from 2 to rank(F)+1 do (
	L#i=indMap(i-1,b,phi);
);
L#1=dual gens kernel (dual L#2);
return chainComplex((toList(L))_{1..#toList(L)-1});
)



--this computes the module K_b^a (F). It is the Weyl module corresponding
--to the appropriate hook partition.
dualhook = (a,b,F) -> (C:=extComult(a-1,1,F);
M:=map(symmetricPower(b-1,F),symmetricPower(1,F)**symmetricPower(b,F),divContract(1,b,F));
K:=map(exteriorPower(a-1,F)**symmetricPower(b-1,F),exteriorPower(a,F)**symmetricPower(b,F),(id_(exteriorPower(a-1,F))**M)*(C**id_(symmetricPower(b,F))));
L:= kernel K;
return L;
)

--this computes the induced map K_b^a (F) -> K_b^a-1 (F). phi is a map F -> R.
indMapDual= (a,b,phi) -> (F:= source phi;
    T:=(koszul(a,phi)**id_(symmetricPower(b,F)))*inducedMap(exteriorPower(a,F)**symmetricPower(b,F),dualhook(a,b,F));
I:=inducedMap(exteriorPower(a-1,F)**symmetricPower(b,F),dualhook(a-1,b,F));
N:=prune(T//I);
return N;)


--this computes the K-complexes of Buchsbaum-Eisenbud. These are dual to the L-complexes.
Kcomplex = (b,phi) -> (F:= source phi;
    K=new MutableList;
for i from 0 to rank F-1 do (
	K#(#K)=(-1)^i*indMapDual(i,b,phi);
);
K#(#K)=gens kernel K#(rank F-1 );
return chainComplex((toList(K))_{1..#(toList(K))-1});
)

--======================================================================================================
--The Eagon-Northcott family of complexes

--M is an arbitrary matrix, viewed as a map F -> G. This computes the induced map
--wedge^a F ** D_b (G^*) -> wedge^(a-1) F ** D_(b-1) (G^*).
ENback = (a,b,M) -> (F:=source M;
G:=target M;
C:=extComult(a-1,1,F);
D:=divContract(1,b,G);
H:=(id_(exteriorPower(a-1,F))**D)*(id_(exteriorPower(a-1,F))**M**symmetricPower(b,G))*(C**id_(symmetricPower(b,G)));
return H;
)

--M is an arbitrary matrix, viewed as a map F -> G. This computes the induced map
--wedge^a F ** S_b (G) -> wedge^(a-1) F ** S_(b+1) (G).
ENfront = (a,b,M) -> (F:=source M;
G:=target M;
C:=extComult(a-1,1,F);
D:=map(symmetricPower(b+1,G),G**symmetricPower(b,G),dual divComult(1,b,G));
H:=(id_(exteriorPower(a-1,F))**D)*(id_(exteriorPower(a-1,F))**M**symmetricPower(b,G))*(C**id_(symmetricPower(b,G)));
return H;
)

--M is an arbitrary matrix, viewed as a map F -> G. This computes the induced map
--wedge^a F -> wedge^(a-g) F ** 
ENbump = (a,M) -> (F:=source M;
G:=target M;
g:=rank(G);
C:=extComult(a-g,g,F);
D:=map(symmetricPower(0,G),exteriorPower(g,G)**exteriorPower(g,dual(G)),extContract(g,g,dual G));
H:=(id_(exteriorPower(a-g,F))**D)*(id_(exteriorPower(a-g,F))**exteriorPower(g,M))*(C);
return H;
)


--this computes the Eagon-Northcott family of complexes, as in Eisenbud's book
--on page 595. The inputs are an integer n and a matrix M.
EN = (n,M) -> (f:=rank(source M);
g:=rank(target M);
ENo = new MutableList;
if n>=0 then (
for i from 0 to f-g-n do (
	ENo#(f-g+1-i)=(-1)^(f-g+1-i)*ENback(f-i,f-g-n-i,M);
);
ENo#(1+n)=ENbump(g+n,M);
if n>0 then (
for i from 1 to n do (
	ENo#(i)=ENfront(i,n-i,M);
);
);
ENO=chainComplex((toList(ENo))_{1..#toList(ENo)-1});
);
if n<0 then (
    for i from 0 to f-g-n do (
	ENo#(#ENo) =(-1)^(f-g+1-i)*ENback(f-i,f-g-n-i,M);
	);
    ENO=chainComplex(reverse (toList(ENo))_{0..#toList(ENo)-2})[-n];
    );
return removeZeroTrailingTerms ENO;
)

--========================================================================================================
--Generalized Eagon-Northcott complexes coming from "Linear Strands of Determinantal Facet Ideals"

--EN basis list, P is the list of facets
--in wedge^a times sym_b
genENBlist = (a,b,P,G) -> (P1 := flatten apply(P,j->subsets(j,a));
g=rank(G);
P2:=compositions(g,b);
P3:=new MutableList;
for i in P1 do (
	for j in P2 do (
		P3#(#P3)=(i,j);
	);
);
unique toList P3
)


--this computes the projection from wedge^a times D_b to the
--free submodule generated by all the generalized EN basis elts.
genENprojection = (a,b,P,G) -> (f=max (P/max);
g= rank(G);
K1:=subsets(1..f,a);
K2:=compositions(g,b);
K3=new MutableList;
for i in K1 do (
	for j in K2 do (
		K3#(#K3)=(i,j);
	);
);
K3=toList K3;
bList = genENBlist(a,b,P,G);
mutMat = mutableMatrix(ring G,length bList, length K3);
for j from 0 to length(K3)-1 do (
	if member(K3_j,set(bList)) then (
	mutMat_(position(bList,i->(i==K3_j)),j)=1;
	);
);
matrix mutMat
)

--generalized Eagon-Northcott differential
genENdiff = (a,b,M,P) -> (if rank(genENprojection(a-1,b-1,P,target M))>0 and rank(genENprojection(a,b,P,target M))>0 then (
matrix entries transpose (matrix entries transpose(genENprojection(a-1,b-1,P,target M)*ENback(a,b,M))//(matrix entries transpose genENprojection(a,b,P,target M)))
)
)

--this is the generalized Eagon-Northcott complex associated to a determinantal facet ideal.
--the integer i specifies the size of the minors, and P should be entered as a list of lists,
--where each entry is a clique of the associated simplicial complex
genEN = (i,P) -> (M=genMat(i,max (P/max));
f:=rank(source M);
g:=rank(target M);
ENo = new MutableList;
for i from 0 to f-g do (
	ENo#(#ENo)=genENdiff(f-i,f-g-i,M,P);
);
ENO=chainComplex(delete(null,reverse toList ENo));
return ENO;
)

--=================================================================================================================
--Complexes of "Canonical Complexes Associated to a Matrix" by Andy Kustin

--phi is a map G -> F. This is the induced map
--appearing in the complex L_phi^(a,b) of Andy Kustin's paper
--"Canonical Complexes Associated to a Matrix" 
indMapn= (k,a,b,phi) -> (F:=target phi;
G:=source phi;
m1:=(id_(exteriorPower(k,F))**(dual divContract(1,1,G))**id_(exteriorPower(a,G)**symmetricPower(b,G)));
m2:=(id_(exteriorPower(k,F))**phi**(reverseFactors(G,exteriorPower(a,G),0,0))**id_(symmetricPower(b,G)));
m3:=(wedgeProduct(k,1,F)**id_(exteriorPower(a,G))**symProduct(1,b,G));
T:=m3*m2*m1*inducedMap(exteriorPower(k,F)**exteriorPower(a,G)**symmetricPower(b,G),exteriorPower(k,F)**hook(a,b,G));
I:=inducedMap(exteriorPower(k+1,F)**exteriorPower(a,G)**symmetricPower(b+1,G),exteriorPower(k+1,F)**hook(a,b+1,G));
N:=prune (T//I);
return N;)

--This is the complex L_phi^(a,b) of the aforementioned paper
genLcomplex = (N,a,phi) -> (L:=new MutableList;
for i from 1 to rank(target phi)-N-1 do (
	L#(i-1)=indMapn(N+i,a,i,phi);
);
return chainComplex(reverse (toList(L)));
)


--this computes the contraction map
--S_r (F) ** D_n (F) --> D_(n-r) F.
divContract = (r,n,F) -> (D:=divComult(r,n-r,F);
M:=map((ring(F))^1,symmetricPower(r,dual F)**symmetricPower(r,F),matrix{flatten entries id_(symmetricPower(r,F))});
a:=id_(symmetricPower(r,dual(F)))**D;
b:=M**id_(symmetricPower(n-r,F));
C:=b*a;
return C;
)




--this computes the Weyl module associated to a hook partition in a different way.
--the resulting module is off from the output of the other dualhook code by an index shift,
--but this is how Andy computes the modules in his paper.
ndualhook = (a,b,F) -> (m1:=id_(exteriorPower(a,F))**dual(extContract(1,1,F))**id_(symmetricPower(b,F));
m2:=wedgeProduct(a,1,F)**divContract(1,b,F);
return kernel (m2*m1)
)

--this is the induced map appearing in the K_phi^(a,b) complex.
indMapnDual= (k,a,b,phi) -> (F:=target phi;
G:=source phi;
m1:=(id_(exteriorPower(k,F))**(dual divContract(1,1,G))**id_(exteriorPower(a,G)**symmetricPower(b,G)));
m2:=(id_(exteriorPower(k,F))**phi**(reverseFactors(G,exteriorPower(a,G),0,0))**id_(symmetricPower(b,G)));
m3:=(wedgeProduct(k,1,F)**id_(exteriorPower(a,G))**divContract(1,b,G));
T:=m3*m2*m1*inducedMap(exteriorPower(k,F)**exteriorPower(a,G)**symmetricPower(b,G),exteriorPower(k,F)**ndualhook(a,b,G));
I:=inducedMap(exteriorPower(k+1,F)**exteriorPower(a,G)**symmetricPower(b-1,G),exteriorPower(k+1,F)**ndualhook(a,b-1,G));
N:=prune (T//I);
return N;
)

--this is the complex K_phi^(a,b).
genKcomplex = (N,a,phi) -> (K:=new MutableList;
for i from 0 to N do (
	K#i=indMapnDual(i,a,N-i,phi);
);
return chainComplex(reverse (toList(K)));
)


indBackdiff = (a,q,b,l,phi) -> (G:=source phi;
F:=target phi;
I1:=inducedMap(ambient(hook(a,q,G)),hook(a,q,G))**inducedMap(ambient(hook(b,l,F)),hook(b,l,F));
I2:=inducedMap(ambient(hook(a,q-1,G)),hook(a,q-1,G))**inducedMap(ambient(hook(b+1,l,F)),hook(b+1,l,F));
m1:=id_(exteriorPower(a,G)**symmetricPower(q,G))**(dual divContract(1,1,G))**id_(exteriorPower(b,F)**symmetricPower(l-1,F));
m2:=id_(exteriorPower(a,G))**revdivContract(q,1,G)**phi**id_(exteriorPower(b,F)**symmetricPower(l-1,F));
m3:=id_(exteriorPower(a,G)**symmetricPower(q-1,G))**wedgeProduct(1,b,F)**id_(symmetricPower(l-1,F));
phiI:=m3*m2*m1;
return prune ((phiI*I1)//I2)
)

indFrontdiff = (a,q,b,l,phi) -> (G:=source phi;
F:=target phi;
I1:=inducedMap(ambient(hook(a,q,G)),hook(a,q,G))**inducedMap(ambient(hook(b,l,F)),hook(b,l,F));
I2:=inducedMap(ambient(hook(a,q+1,G)),hook(a,q+1,G))**inducedMap(ambient(hook(b+1,l,F)),hook(b+1,l,F));
m1:=id_(exteriorPower(a,G)**symmetricPower(q,G))**dual(extContract(1,1,G))**id_(exteriorPower(b,F)**symmetricPower(l-1,F));
m2:=id_(exteriorPower(a,G))**symProduct(q,1,G)**phi**id_(exteriorPower(b,F)**symmetricPower(l-1,F));
m3:=id_(exteriorPower(a,G)**symmetricPower(q+1,G))**wedgeProduct(1,b,F)**id_(symmetricPower(l-1,F));
phiI:=m3*m2*m1;
return prune ((phiI*I1)//I2)
)

--for the definition of the next two, see Definition 7.2 of the "Canonical Complexes ..."
--paper by Kustin.
--map connecting K_phi to the appropriate exterior power
frontBump = (b,a,l,phi) -> (G:=source phi;
F:=target phi;
I1:=inducedMap(ambient(hook(a,l,F)),hook(a,l,F));
I2:=id_(exteriorPower(b,G))**inducedMap(ambient(hook(a+b,l,F)),hook(a+b,l,F));
m1:=dual(extContract(1,1,exteriorPower(b,G)))**id_(exteriorPower(a,F)**symmetricPower(l,F));
m2:=id_(exteriorPower(b,G))**exteriorPower(b,phi)**id_(exteriorPower(a,F)**symmetricPower(l,F));
m3:=id_(exteriorPower(b,G))**wedgeProduct(b,a,F)**id_(symmetricPower(l,F));
phiI:=m3*m2*m1;
return prune ((phiI*I1)//I2)
)

--map connecting the appropriate exterior power to L_phi
backBump = (b,a,l,phi) -> (G:=source phi;
F:=target phi;
I1:=id_(exteriorPower(b,G))**inducedMap(ambient(hook(a,l,F)),hook(a,l,F));
I2=inducedMap(ambient(hook(a+b,l,F)),hook(a+b,l,F));
m1:=exteriorPower(b,phi)**id_(exteriorPower(a,F)**symmetricPower(l,F));
m2:=wedgeProduct(b,a,F)**id_(symmetricPower(l,F));
phiI:=m2*m1;
return prune ((phiI*I1)//I2)
)

--this is a generalized version of the C complex of the paper by Kustin.
--phi: G -> F, and rank(G) <= rank(F). The integers i, a, and l are parameters
--that give different complexes as they are varied.
genCcomplex = (i,a,l,phi) -> (F:=target phi;
G:=source phi;
f:=rank F;
g:=rank G;
C:=new MutableList;
for j from 1 to i do (C#(#C)=indFrontdiff(g-a,i+1-j,f-j,l,phi););
C#(#C)=frontBump(g+1-a,f-g+a-i-1,l,phi);
C#(#C)=backBump(a,f-g-i-1,l,phi);
return chainComplex(toList C);
)


--this is the C complex of the paper by Kustin.
--phi: G -> F, and rank(G) <= rank(F). The integers i and a are parameters
--that give different complexes as they are varied.
Ccomplex = (i,a,phi) -> (F:=target phi;
G:=source phi;
f:=rank F;
g:=rank G;
C:=new MutableList;
for j from 1 to i do (
	C#(j-1)=indMapn(f-j,g-a,i+1-j,phi);
);
n1:=id_(exteriorPower(f-g+a-i-1,F))**(dual divContract(1,1,exteriorPower(g+1-a,G)));
l1:=prune ((extComult(g-a,1,G))//inducedMap(exteriorPower(g-a,G)**G,hook(g-a,1,G)));
n2:=id_(exteriorPower(f-g+a-i-1,F))**exteriorPower(g+1-a,phi)**l1;
n3:=wedgeProduct(f-g+a-i-1,g+1-a,F)**id_(prune hook(g-a,1,G));
C#(i)=n3*n2*n1;
C#(i+1)=(wedgeProduct(f-g-i-1,a,F))*(id_(exteriorPower(f-g-i-1,F))**exteriorPower(a,phi));
for j from i+2 to f-g+1 do (
	C#(j)=indMapnDual(f-g-j,a,j-i-1,phi);
);
return chainComplex((toList(C)));
)

--============================================================================================================
--Tate-like complex associated to the resolution of a cyclic module

--this code subtracts 1 from L in a way that mimics the
--rank drops in the tate-like complex
sub1s = L -> (L1 := new MutableList;
for i from 0 to #L-1 do (
	L2=new MutableList from L;
	if L2#i>0 and i>0 then (
	L2#i=L2#i-1;
	L2#(i-1)=L2#(i-1)+1;
	L1#(#L1)=toList L2; );
	if L2#i>0 and i==0 then (
	L2#i=L2#i-1;
	L1#(#L1)=toList L2; ); 
);
return rsort toList(L1);
)


--this gives the ranks for the tate-like complex
tateRanks = (L,n) -> (Ln=new MutableList;
Ln#(#Ln)={L};
if n>1 then (
for i from 2 to n do (
	Ln#(#Ln) = unique flatten (apply(Ln#(#Ln-1),i->sort sub1s i));
);
);
return toList Ln;
)


--this gives the tate-like module
tateModule = (F,L)->(M:=exteriorPower(L_0,F_1);
for i from 1 to #L-1 do (if even(i) then (M=M**exteriorPower(L_i,F_(i+1)));
			if odd(i) then (M=M**symmetricPower(L_i,F_(i+1)));
);
return M;
)


--this is just the koszul differential
kosDiff = (n,d) -> ((d**id_(exteriorPower(n-1,source d))*extComult(1,n-1,source d)
)

--computes induced differential 
--D : D_n (F_k) -> F_(k-1)**D_(n-1) (F_k)
divDiff = (n,d) -> ((d**id_(symmetricPower(n-1,source d))*divComult(1,n-1,source d)
)

--gives the appropriate differential between tate-like modules
tateDiff = (L,F) -> (Ln=new MutableList;
for i from 0 to length(L)-1 do (
	if i==0 then Ln#(#Ln) =  ((F.dd_(i+1)**id_(exteriorPower(L_i-1,F_(i+1))))*extComult(1,L_i-1,F_(i+1)) ,(length(F):0),replace(i,0,L),replace(i,L_i-1,L));
	if (i>0 and even(i)) then (m1=divProduct(L_(i-1),1,F_(i))**id_(exteriorPower(L_i-1,F_(i+1)));
			m2=id_(symmetricPower(L_(i-1),F_i))**F.dd_(i+1)**id_(exteriorPower(L_i-1,F_(i+1)));
			m3=id_(symmetricPower(L_(i-1),F_i))**extComult(1,L_i-1,F_(i+1));
			difsgn=sum toList apply(0..i-1,l->(l+1)*L_l);
			Lv=replace(i-1,0,replace(i,0,L));
			Lb1=new MutableList from (length(F):0);
			Lb2=new MutableList from (length(F):0);
			for t from 0 to i-1 do (
				Lb1#t=Lv_t;
			);
			for t from 0 to (length(F)-i) do (
				Lb2#(length(F)-t-1)=Lv_(length(F)-t-1);
			);
			Ln#(#Ln)=((-1)^difsgn*m1*m2*m3,toList Lb1,toList Lb2,replace(i-1,L_(i-1)+1,replace(i,L_i-1,L)));
		);
	if odd(i) then (m1=wedgeProduct(L_(i-1),1,F_(i))**id_(symmetricPower(L_i-1,F_(i+1)));
			m2=id_(exteriorPower(L_(i-1),F_i))**F.dd_(i+1)**id_(symmetricPower(L_i-1,F_(i+1)));
			m3=id_(exteriorPower(L_(i-1),F_i))**divComult(1,L_i-1,F_(i+1));
			difsgn=sum toList apply(0..i-1,l->(l+1)*L_l);
			Lv=replace(i-1,0,replace(i,0,L));
			Lb1=new MutableList from (length(F):0);
			Lb2=new MutableList from (length(F):0);
			for t from 0 to i-1 do (
				Lb1#t=Lv_t;
			);
			for t from 0 to (length(F)-i) do (
				Lb2#(length(F)-t-1)=Lv_(length(F)-t-1);
			);
			Ln#(#Ln)=((-1)^difsgn*m1*m2*m3,toList Lb1,toList Lb2, replace(i-1,L_(i-1)+1,replace(i,L_i-1,L)));
		);
);
apply(Ln,(l,j,k,k1)->(k1,(id_(tateModule(F,j))**l**id_(tateModule(F,k)))))
)



--Build a Tate like complex on a chain complex
--F, with F_0 = R. L is a list of length = length F specifying
--the appropriate exterior/divided powers for the
--topmost term of the complex. n is the homological 
--degree you want the topmost term to appear in.
tateLike = method()
tateLike(ChainComplex,List,ZZ) := (F,L,n) -> (n=n+1;
D:=tateRanks(L,n);
D1=new MutableList;
for l from 0 to length(D)-2 do (
	M1=new MutableList;
	for i from 0 to length(D_l)-1 do (
		L1=new MutableList from (length(D_(l+1)):0);
		Lp=tateDiff((D_l)_i,F);
		for j from 0 to length(toList Lp)-1 do (
			for k from 0 to length(D_(l+1))-1 do (
				if (D_(l+1))_k==(Lp#j)_0 then (L1#(k)=transpose ((Lp#j)_1);
				);
			);
		);
		Col=(toList L1);
		M1#(#M1)=Col;
	);
	D1#(#D1)=matrix toList M1;
);
return chainComplex((reverse toList(D1))/transpose);
)

--puts the Tatelike complex in the standard homological degree; disadvantage is that
--we don't always want the entire Tatelike complex
tateLike(ChainComplex,List) := (F,L) -> (n=0;
    for i from 1 to length(L) do n=n+i*L_(i-1);
    tateLike(F,L,n)
    )

--===============================================================================================================
--Gulliksen-Negard complex

--this outputs the cofactor matrix of a matrix M.
cofactor = M -> (n:=rank source M;
return determinant(M)*id_(Q^n)//M;
)

--give a matrix M, this outputs the corresponding vector
mat2Vec = M -> (matrix {flatten entries M})

--given a vector v of length n^2, this outputs the corresponding square matrix
vec2Mat = (v,n) -> (matrix pack(n,flatten entries v))


--the next few codes are needed for defining Gulliksen-Negard.
psiGN = M -> (n:=rank source M;
M1=new MutableList;
for i from 0 to n-1 do (
	for j from 0 to n-1 do (
		F=mutableMatrix(ring M,n,n);
		F_(i,j)=1;
		M1#(#M1)=(mat2Vec(M*matrix(F)))|(mat2Vec(matrix(F)*M));
	);
);
return matrix entries transpose matrix pack(1,toList M1);
)

phiGN = M -> (n:=rank source M;
M1=new MutableList;
for i from 0 to n-1 do (
	for j from 0 to n-1 do (
		F=mutableMatrix(ring M,n,n);
		F_(i,j)=1;
		M1#(#M1)=mat2Vec(matrix(F)*M);
	);
);
for i from 0 to n-1 do (
	for j from 0 to n-1 do (
		F=mutableMatrix(ring M,n,n);
		F_(i,j)=-1;
		M1#(#M1)=mat2Vec(M*matrix(F));
	);
);
return matrix entries transpose matrix pack(1,toList M1);
)

d1GN = M-> (n:=rank source M;
Mtilde := cofactor(M);
M1=new MutableList;
for i from 0 to n-1 do (
	for j from 0 to n-1 do (
		F=mutableMatrix(ring M,n,n);
		F_(i,j)=1;
		M1#(#M1)=trace(Mtilde*matrix(F));
	);
);
return matrix entries transpose matrix pack(1,toList M1);
)

--Given an arbitrary square matrix M, this outputs
--the classical Gulliksen-Negard complex. The construction used here
--is taken directly from Bruns-Vetter "Determinantal Rings".
gulliksenNegard = M -> (n:=rank source M;
i := transpose (mat2Vec(id_(Q^n))|mat2Vec(id_(Q^n)));
pii:=(mat2Vec(id_(Q^n))|mat2Vec(-id_(Q^n)));
I1:=inducedMap(ambient image psiGN(M),image psiGN(M));
I2:=inducedMap(ambient kernel pii,kernel pii);
ind1:=I1//I2;
ind2:=inducedMap((ker pii)/(image i),ker pii);
d3:=prune (ind2*ind1);
ind3:=phiGN(M)*inducedMap(ambient ker pii,ker pii);
d2:=matrix entries transpose(transpose(prune ind3)//transpose(prune ind2));
d4:=transpose mat2Vec(cofactor(M));
return chainComplex({d1GN(M),d2,d3,d4});
)


--=======================================================================================================
--Some useful codes that I use in the examples

--Define the coordinate ring/generic nxn skew symmetric matrix
pfaffs = n -> (L={};
for i from 1 to n-1 do (L=L|{x_(i,i+1)..x_(i,n)};
);
Q=QQ[splice L];
genericSkewMatrix(Q,n)
--F=resBE(M);
)

--define coordinate ring/generic nxm matrix
genMat = (n,m) -> (L={};
for i from 1 to n do ( L=L|{x_(i,1)..x_(i,m)};
);
Q=QQ[splice L,MonomialOrder=>Lex];
Mmut:=mutableMatrix(Q,n,m);
for i from 0 to n-1 do (
	for j from 0 to m-1 do (
		Mmut_(i,j)=x_(i+1,j+1);
	);
);
matrix Mmut
)

--=============================================================================================================
--Examples


--Running free module:
F=QQ^4;
M1=extComult(2,1,F) 
M2=transpose wedgeProduct(2,1,F)
M1==M2
M1=symProduct(1,1,F)
M2=divComult(1,1,F)
M1==transpose M2
M1=symComult(1,2,F)
M2=divProduct(1,2,F)
M1==transpose M2


--L/K complex examples

Q=QQ[x_1..x_3];
--minimal free resolution of homogeneous maximal ideal cubed
Lcomplex(3,vars Q)
genList(2,3,QQ^3)
genList(1,3,QQ^3)
--The L/K complexes are dual to each other
extend(Lcomplex(2,vars Q),(dual Kcomplex(1,vars Q))[-3],id_(Q^1)) --is an isomorphism

--examples for Eagon-Northcott
M=genMat(2,4)
--the standard classical Eagon-Northcott complex
EN(0,M)
--the classical Buchsbaum-Rim complex
EN(1,M)
extend(EN(0,M),(dual EN(2,M))[-3],id_(Q^1))
--a class of self dual complexes
EN(2,genMat(2,6))
EN(3,genMat(3,9))
EN(1,genMat(3,5))


--generalized Eagon-Northcott complexes
genEN(2,{{1,2,3,4},{2,3,4,5}})
prune HH oo
--the fact that homology is 0 shows that the complex has no
--minimal nonfaces

--canonical complexes examples (this is example 7.7 of Kustin's paper)
M=genMat(4,3)
C=Ccomplex(0,2,M)
C.dd

--Tatelike complexes
Q=QQ[x_1..x_3]
F=Lcomplex(2,vars Q);
TL1 = tateLike(F,{2,0,0})
--let's use this to find a DG algebra structure on F
(extend(F,TL1,id_(Q^1)))_2 --this gives the product between deg 1 elements
genList(0,2,QQ^3)
genList(1,2,QQ^3)
genList(2,2,QQ^3)
TL2 = tateLike(F,{1,1,0},3)
(extend(F,TL2,id_(Q^1)))_3 --this is the product between deg 1 and 2 elements


--Add in some computations of Buchsbaum-Eisenbud multipliers

--=======================================================================
--Multilinear algebra codes

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


--this function computes the comultiplication
--map on the divided power algebra 
--D_(p+q) F -> D_p F ** D_q F with respect
--to the GrRevLex basis order
divComult = (p,q,F) -> (n:=rank(F);
L:=compositions(n,p+q);
L1:=compositions(n,p);
L2:=compositions(n,q);
L3:=new MutableList;
for i from 0 to (#L1)-1 do (
	for j from 0 to (#L2)-1 do (
		L3#(#L3)=(L1_i,L2_j);
	);
);
L3=toList L3;
Mmut:=mutableMatrix(ring(F),#L3,#L);
for i from 0 to (#L)-1 do (
	for j from 0 to (#L3)-1 do (
		for k in L1 do (
			if (k,L_i-k)==L3_j then Mmut_(j,i)=1;
		);
	);
);
return map(symmetricPower(p,F)**symmetricPower(q,F),symmetricPower(p+q,F),matrix Mmut);
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
genList = (a,b,F) -> (L:=matrix gens hook(a,b,F);
n:= rank source L;
m:=rank target L;
L1:=flatten table(subsets(1..rank F,a),compositions(rank F,b),(u,v)->{u,v});
L2:=new MutableList;
for i from 0 to n-1 do (
	for j from 0 to m-1 do (
		if not((L_{i})_(j,0)==0) then (L2#(#L2)=flatten {(L_{i})_(j,0),L1_j};);
	);
);
MatrixExpression pack(a+1,toList L2)
)

--this computes the map L_b^a(F) -> L_b^a-1 (F) induced by the map phi.
--phi here must be a map F -> R.
indMap= (a,b,phi) -> (F:=source phi;
T:=(-1)^a*(koszul(a,phi)**id_(symmetricPower(b,F)))*inducedMap(exteriorPower(a,F)**symmetricPower(b,F),hook(a,b,F));
I:=(-1)^(a-1)*inducedMap(exteriorPower(a-1,F)**symmetricPower(b,F),hook(a-1,b,F));
N:=prune(T//I);
N)

--this gives the standard L-complex of Buchsbaum-Eisenbud. Here phi : F -> R;
--if grade (im phi) = rank F, this is a minimal free resolution of (im phi)^b.
Lcomplex = (b,phi) -> (F:=source phi;
L:=new MutableList;
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
    K:=new MutableList;
for i from 0 to rank F-1 do (
	K#(#K)=(-1)^i*indMapDual(i,b,phi);
);
K#(#K)=gens kernel K#(rank F-1 );
return chainComplex((toList(K))_{1..#(toList(K))-1});
)

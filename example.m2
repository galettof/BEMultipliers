-- example: Koszul complex
restart
needsPackage "BEMultipliers"
A=QQ[x,y,z]
K=koszul vars A
-- compute multipliers, returns a list
a=bem(K)
-- can also see individual ones
bem(K,2)
-- let's check the structure theorem
exteriorPower(rank K.dd_2,K.dd_2) ==
bem(K,2) * dualMultiplier(K,3)
exteriorPower(rank K.dd_1,K.dd_1) ==
bem(K,1) * dualMultiplier(K,2)


-- example: Eagon-Northcott of 2x4 matrix
restart
needsPackage "BEMultipliers"
A=QQ[x_(1,1)..x_(2,4)]
G=transpose genericMatrix(A,4,2)
I=minors(2,G)
RI=res I
-- compute multipliers, returns a list
a=bem(RI);
-- can also see individual ones
bem(RI,2)
-- let's check the structure theorem
exteriorPower(rank RI.dd_2,RI.dd_2) ==
bem(RI,2) * dualMultiplier(RI,3)
exteriorPower(rank RI.dd_1,RI.dd_1) ==
bem(RI,1) * dualMultiplier(RI,2)


-- example: all squarefree monomials of degree 2
restart
needsPackage "BEMultipliers"
A=QQ[x_1..x_4]
I=ideal(subsets(gens A,2)/product)
RI=res I
-- compute multipliers, returns a list
a=bem(RI);
-- can also see individual ones
bem(RI,2)
-- let's check the structure theorem
exteriorPower(rank RI.dd_2,RI.dd_2) ==
bem(RI,2) * dualMultiplier(RI,3)
exteriorPower(rank RI.dd_1,RI.dd_1) ==
bem(RI,1) * dualMultiplier(RI,2)


-- example: B-E's resolution of 6x6 pfaffians of 7x7 skew matrix
restart
needsPackage "BEMultipliers"
A=QQ[x_(1,2)..x_(1,7),x_(2,3)..x_(2,7),x_(3,4)..x_(3,7),
    x_(4,5)..x_(4,7),x_(5,6)..x_(5,7),x_(6,7)]
G=genericSkewMatrix(A,7)
I=pfaffians(6,G);
RI=res I
-- compute multipliers, returns a list
elapsedTime a=bem(RI);
-- can also see individual ones
bem(RI,2)
-- let's check the structure theorem
exteriorPower(rank RI.dd_2,RI.dd_2) ==
bem(RI,2) * dualMultiplier(RI,3)
exteriorPower(rank RI.dd_1,RI.dd_1) ==
bem(RI,1) * dualMultiplier(RI,2)

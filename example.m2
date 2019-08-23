-- example with Koszul complex
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
bem(K,2) * dualMultiplier(K,3) * exteriorDuality(K,2)
exteriorPower(rank K.dd_1,K.dd_1) ==
bem(K,1) * dualMultiplier(K,2) * exteriorDuality(K,1)


-- example with Eagon-Northcott of 2x4 matrix
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
bem(RI,2) * dualMultiplier(RI,3) * exteriorDuality(RI,2)
exteriorPower(rank RI.dd_1,RI.dd_1) ==
bem(RI,1) * dualMultiplier(RI,2) * exteriorDuality(RI,1)


-- example with monomial ideal
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
bem(RI,2) * dualMultiplier(RI,3) * exteriorDuality(RI,2)
exteriorPower(rank RI.dd_1,RI.dd_1) ==
bem(RI,1) * dualMultiplier(RI,2) * exteriorDuality(RI,1)


-- example with larger Eagon-Northcott
restart
needsPackage "BEMultipliers"
A=QQ[x_(1,1)..x_(4,6)]
G=transpose genericMatrix(A,6,4)
I=minors(4,G)
RI=res I
-- compute multipliers, returns a list
a=bem(RI);
-- can also see individual ones
bem(RI,2)
-- let's check the structure theorem
exteriorPower(rank RI.dd_2,RI.dd_2) ==
bem(RI,2) * dualMultiplier(RI,3) * exteriorDuality(RI,2)
exteriorPower(rank RI.dd_1,RI.dd_1) ==
bem(RI,1) * dualMultiplier(RI,2) * exteriorDuality(RI,1)

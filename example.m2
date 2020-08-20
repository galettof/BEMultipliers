-- example: Koszul complex
restart
needsPackage "BuchsbaumEisenbudMultipliers"
A=QQ[x,y,z]
K=koszul vars A
-- get last multiplier with full rank computation
elapsedTime aMultiplier(1,K,ComputeRanks=>true)
-- remove cache
remove(K.cache,aMultiplier)
-- get last multiplier finding rank with Buchsbaum-Eisenbud
elapsedTime aMultiplier(1,K)
-- other multipliers
aMultiplier(3,K)
aMultiplier(2,K)
-- let's check the structure theorem
E2=exteriorPower(rank K_2,K_2)
exteriorPower(rank K.dd_2,K.dd_2) == aMultiplier(2,K) * ((dual aMultiplier(3,K))**E2) * exteriorDuality(rank K.dd_2,K_2)
E1=exteriorPower(rank K_1,K_1)
exteriorPower(rank K.dd_1,K.dd_1) == aMultiplier(1,K) * ((dual aMultiplier(2,K))**E1) * exteriorDuality(rank K.dd_1,K_1)

-- let's compute all lower order multipliers
aMultiplier(1,3,K)
lowerBEM(0,1,K)
lowerBEM(0,2,K)
lowerBEM(0,3,K)
lowerBEM(1,3,K)

-- example: Eagon-Northcott of 2x4 matrix
restart
needsPackage "BuchsbaumEisenbudMultipliers"
A=QQ[x_(1,1)..x_(2,4)]
G=transpose genericMatrix(A,4,2)
I=minors(2,G)
RI=res I
-- get last multiplier with full rank computation
elapsedTime aMultiplier(1,RI,ComputeRanks=>true)
-- remove cache
remove(RI.cache,aMultiplier)
-- get last multiplier finding rank with Buchsbaum-Eisenbud
elapsedTime aMultiplier(1,RI)
-- other multipliers
aMultiplier(3,RI)
aMultiplier(2,RI)
-- let's check the structure theorem
E2=exteriorPower(rank RI_2,RI_2)
exteriorPower(rank RI.dd_2,RI.dd_2) == aMultiplier(2,RI) * ((dual aMultiplier(3,RI))**E2) * exteriorDuality(rank RI.dd_2,RI_2)
E1=exteriorPower(rank RI_1,RI_1)
exteriorPower(rank RI.dd_1,RI.dd_1) == aMultiplier(1,RI) * ((dual aMultiplier(2,RI))**E1) * exteriorDuality(rank RI.dd_1,RI_1)


-- example: all squarefree monomials of degree 2
restart
needsPackage "BuchsbaumEisenbudMultipliers"
A=QQ[x_1..x_4]
I=ideal(subsets(gens A,2)/product)
RI=res I
-- get last multiplier with full rank computation
elapsedTime aMultiplier(1,RI,ComputeRanks=>true)
-- remove cache
remove(RI.cache,aMultiplier)
-- get last multiplier finding rank with Buchsbaum-Eisenbud
elapsedTime aMultiplier(1,RI)
-- other multipliers
aMultiplier(3,RI)
aMultiplier(2,RI)
-- let's check the structure theorem
E2=exteriorPower(rank RI_2,RI_2)
exteriorPower(rank RI.dd_2,RI.dd_2) == aMultiplier(2,RI) * ((dual aMultiplier(3,RI))**E2) * exteriorDuality(rank RI.dd_2,RI_2)
E1=exteriorPower(rank RI_1,RI_1)
exteriorPower(rank RI.dd_1,RI.dd_1) == aMultiplier(1,RI) * ((dual aMultiplier(2,RI))**E1) * exteriorDuality(rank RI.dd_1,RI_1)


-- example: B-E's resolution of 6x6 pfaffians of 7x7 skew matrix
restart
needsPackage "BuchsbaumEisenbudMultipliers"
A=QQ[x_(1,2)..x_(1,7),x_(2,3)..x_(2,7),x_(3,4)..x_(3,7),
    x_(4,5)..x_(4,7),x_(5,6)..x_(5,7),x_(6,7)]
G=genericSkewMatrix(A,7)
I=pfaffians(6,G);
RI=res I
-- get last multiplier with full rank computation
elapsedTime aMultiplier(1,RI,ComputeRanks=>true)
-- remove cache
remove(RI.cache,aMultiplier)
-- get last multiplier finding rank with Buchsbaum-Eisenbud
elapsedTime aMultiplier(1,RI)
-- other multipliers
aMultiplier(3,RI)
aMultiplier(2,RI)
-- let's check the structure theorem
E2=exteriorPower(rank RI_2,RI_2)
exteriorPower(rank RI.dd_2,RI.dd_2) == aMultiplier(2,RI) * ((dual aMultiplier(3,RI))**E2) * exteriorDuality(rank RI.dd_2,RI_2)
E1=exteriorPower(rank RI_1,RI_1)
exteriorPower(rank RI.dd_1,RI.dd_1) == aMultiplier(1,RI) * ((dual aMultiplier(2,RI))**E1) * exteriorDuality(rank RI.dd_1,RI_1)

-- example: bigger Koszul complex
restart
needsPackage "BuchsbaumEisenbudMultipliers"
A=QQ[x_1..x_4]
K=koszul vars A
-- compute multipliers, returns a list
a=bem(K)
-- can also see individual ones
bem(2,K)
-- let's check the structure theorem
exteriorPower(rank K.dd_2,K.dd_2) ==
bem(2,K) * dualMultiplier(3,K)
exteriorPower(rank K.dd_1,K.dd_1) ==
bem(1,K) * dualMultiplier(2,K)

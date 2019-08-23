-- example with Koszul complex
restart
needsPackage "BEMultipliers"
A=QQ[x,y,z]
K=koszul vars A
--a=bem(K,2)
a=bem(K)
w2d2=exteriorPower(2,K.dd_2)
c2=exteriorDuality(K,2)
b3=dualMultiplier(K,3)
a_2*b3*c2 == w2d2


-- example with Eagon-Northcott of 2x4 matrix
restart
needsPackage "BEMultipliers"
A=QQ[x_(1,1)..x_(2,4)]
G=transpose genericMatrix(A,4,2)
I=minors(2,G)
RI=res I
a=bem(RI);
w2d2=exteriorPower(5,RI.dd_2)
c2=exteriorDuality(RI,2)
b3=dualMultiplier(RI,3)
a_2*b3*c2 == w2d2
-- this is still wrong
-- I think it's because I'm not using the dual of the
-- exterior duality in constructing the multipliers

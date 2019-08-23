-- example with Koszul complex
restart
needsPackage "BEMultipliers"
A=QQ[x,y,z]
K=koszul vars A
a=bem(K,2)
bem(K)

-- example with Eagon-Northcott of 2x4 matrix
restart
load "BEmult.m2"
A=QQ[x_(1,1)..x_(2,4)]
G=transpose genericMatrix(A,4,2)
I=minors(2,G)
RI=res I
a=BEmult(RI);
E=exteriorPower(5,RI.dd_2);
M=a_1*(transpose a_2)*inverse(wedgeIso(3,8));


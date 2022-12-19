function modyield=NelsonSiegelYield(startPar,T)


%to make the code more readable
b0=startPar(1); b1=startPar(2); b2=startPar(3); Lum=startPar(4);

modyield=b0+(b1+b2).*(1-exp(-T./Lum))./(T./Lum)-b2.*exp(-T./Lum);

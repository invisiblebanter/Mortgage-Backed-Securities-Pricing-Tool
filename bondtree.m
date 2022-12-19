function r=bondtree(shortTree,T,startParam,dt)
%function [modPrice,modVol]=bondtree(shortTree,T,startParam)
B=ones(T+1,T+1); 
%for replicating a new short tree
sT=shortTree;
%for i=2:T
    for j=1:T
        if j==1
            sT(T,j)=startParam(1);
        else
            sT(T,j)=sT(T,j-1)*exp(2*startParam(2)*sqrt(dt));
        end
    end
%end
% for iterating around the bond tree
for i=T:-1:1
    for j=1:i
        B(i,j)=(0.5*B(i+1,j)+0.5*B(i+1,j+1))*exp(-sT(i,j)*dt);
    end
end
%output for the bondtree function
modPrice=B(1,1);
upVol=-(1/(T-1))*log(B(2,2));
downVol=-(1/(T-1))*log(B(2,1));
modVol=0.5/sqrt(dt)*log(upVol/downVol);
r=[modPrice,modVol];







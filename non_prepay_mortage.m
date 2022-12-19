function output=non_prepay_mortage(principal,mortageRate,n,frequency,shortTree)

r=mortageRate/frequency;            % frequency is refering to monthly or yearly payment
T=n*frequency;                      % n is referring to Years in this case    
dt=1/frequency;                     % single time step  
non_prepay_table=nan(T+1,T+1);      % initialize the table
c=principal*r/(1-1/(1+r)^T);        % calculate each single payment

non_prepay_table(T+1,:)=c;          % the last payment is c

%% -------------------------start iteration the tree------------------------
for i=T:-1:1
    for j=1:i
        if i==1
            non_prepay_table(i,j)=0.5*(non_prepay_table(i+1,j)+non_prepay_table(i+1,j+1))*exp(-shortTree(i,j)*dt);
        else
            non_prepay_table(i,j)=0.5*(non_prepay_table(i+1,j)+non_prepay_table(i+1,j+1))*exp(-shortTree(i,j)*dt)+c;
        end        
    end
end

p0=non_prepay_table(1,1);      
pUp=non_prepay_table(2,2);     
pDown=non_prepay_table(2,1);   
pUU=non_prepay_table(3,3);
pUD=non_prepay_table(3,2);
pDD=non_prepay_table(3,1);

yUU=shortTree(3,3);          %the new short rate tree
yUD=shortTree(3,2);
yDD=shortTree(3,1);
yUp=shortTree(2,2) ;
yDown=shortTree(2,1);

duration=-1/p0*(pUp-pDown) /(yUp-yDown);
duraUp=-1/pUp*(pUU-pUD) / (yUU-yUD);
duraDown=-1/pDown*(pUD-pDD)/(yUD-yDD);
convexity=duration^2-(duraUp-duraDown)/(yUp-yDown);
output=[p0 duration convexity];
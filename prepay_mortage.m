function output=prepay_mortage(principal,mortageRate,n,frequency,shortTree,fixed_table)

r=mortageRate/frequency;        % the input of mortage rate is annual
T=n*frequency;
dt=1/frequency;
c=principal*r/(1-1/(1+r)^T);    % calculate fixed payment amount

%% build up the optimal tree

optimal_Tree=nan(T+1,T+1);
optimal_Tree(end,1:end)=c;
continuation_tree=nan(361,361);                   % continuation tree is the same as the continuation value in the excel file
continuation_tree(end,1:end)=c;
for i=T:-1:1
    for j=i:-1:1
        if i==360
            continuation_tree(i,j)=0.5*(continuation_tree(i+1,j)+continuation_tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(continuation_tree(i,j),fixed_table(i)) + c;   %compare and obtain the smallest value
        elseif i==1
            continuation_tree(i,j)=0.5*(optimal_Tree(i+1,j)+optimal_Tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(continuation_tree(i,j),fixed_table(i));
        else
            continuation_tree(i,j)=0.5*(optimal_Tree(i+1,j)+optimal_Tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(continuation_tree(i,j),fixed_table(i)) + c;
        end
    end
end


p0=optimal_Tree(1,1);       
pUp=optimal_Tree(2,2);      
pDown=optimal_Tree(2,1);    
pUU=optimal_Tree(3,3);
pUD=optimal_Tree(3,2);
pDD=optimal_Tree(3,1);

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
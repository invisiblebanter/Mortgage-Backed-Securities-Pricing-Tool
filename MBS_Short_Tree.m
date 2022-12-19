clear;clc
%% Set parameters

%End date, in years
T=30;

%Number of steps in the binomial tree
N=360;

%set time step size
dt=T/N;

%m for vector of maturities
m=(dt:dt:T)';

%import data for the treasury yield curve

yieldData=dataset('File','yieldAndVolatility2000.csv','Delimiter',',')

%par vector contains, respectively, beta0, beta1, beta2, and tau

startPar=[1,1,1,1];

%Alias the error function with observed data and parameters.

errorNelsonSiegelYield=@(Par) sum((NelsonSiegelYield(Par,yieldData.m(1:end))-yieldData.y(1:end)).^2);
errorNelsonSiegelVol=@(Par) sum((NelsonSiegelVol(Par,yieldData.m(1:end))-yieldData.vol(1:end)).^2);

%Run optimization
[parMinYield,errorMin] = fminsearch(errorNelsonSiegelYield,startPar);
[parMinVol,errorMin] = fminsearch(errorNelsonSiegelVol,startPar);

%Alias new function fittedYield that takes as input some maturity m and 
%returns fitted yields at the estimated parameters
fittedFunctionYield= @(m) NelsonSiegelYield(parMinYield,m);
fittedFunctionVol  = @(m) NelsonSiegelVol(parMinVol,m);

observedData=dataset(m);
observedData.yield=fittedFunctionYield(m);
observedData.price=1./((1+observedData.yield/2).^(2*m))
observedData.volatility=fittedFunctionVol(m);

%Create empty tree. the (i,j) node represents time i and node j of the
%short rate tree
shortTree=NaN(N,N);

%Starting Guesses for fminsearch
startR=.01;
startParam=[.01 .01]

%% Create first node in tree
i=1;
thisPrice=observedData.price(i);
thisMaturity=i;

%No volatility matching for first node
thisError=@(startR) (thisPrice-exp(-startR*dt))^2;

thisR=fminsearch(thisError,startR);

shortTree(i,1)=thisR;

%% Moving on in the tree

for i = 2:N
        i/N
        thisPrice=observedData.price(i);
        thisMaturity=i;
        thisVolatility=observedData.volatility(i);
        thisPriceAndVol=[thisPrice,thisVolatility];
        %T=i; 

        thisError=@(startParam) sum((thisPriceAndVol-bondtree(shortTree,i,startParam,dt)).^2);       
        thisParam=fminsearch(thisError,startParam);
        startParam=thisParam;
        
        %Build next step in tree with the right parameters
        thisBottomR=thisParam(1);
        thisSigma=thisParam(2);
        treeStep=exp(2*thisSigma*sqrt(dt));

        %populate next step in shortTree
        shortTree(i,1)=thisBottomR;
        for j=2:i
            shortTree(i,j)=shortTree(i,j-1)*treeStep;
        end

end

%% ------------------------------------------------------------------------
%------------------------------------building mortage----------------------
mtg=100000000;                                                %mortage 
mtg_r=0.07;                                                   %mortage yearly rate
mtg_mth_r=mtg_r/12;                                           %mortage monthly rate
%mtg_yr_pmt=(mtg*(mtg_r))/(1-(1/(1+mtg_r)^30));               %mortage yearly payment
mtg_mth_pmt=(mtg*(mtg_mth_r))/(1-(1/(1+mtg_mth_r)^360));      %mortage monthly payment

% new_m for vector of maturities and also new table
new_m=(dt:dt:T)';
new_m(2:361)=new_m;
new_m(1)=0;

% creating a fixed table
fixed_table=dataset(new_m);
fixed_table.fixed_monthly_payment=mtg_mth_pmt*ones(361,1);
fixed_table.oustanding_principal=nan(361,1);
fixed_table.int_pmt=nan(361,1);
fixed_table.principal_pmt=nan(361,1);
fixed_table.oustanding_principal(1)=mtg;

for i=2:length(fixed_table)
    fixed_table.int_pmt(i)=mtg_mth_r * fixed_table.oustanding_principal(i-1);
    fixed_table.principal_pmt(i)=mtg_mth_pmt - fixed_table.int_pmt(i);
    fixed_table.oustanding_principal(i)=fixed_table.oustanding_principal(i-1) - fixed_table.principal_pmt(i);
end

% non_prepay table for never refinance
non_prepay_table=nan(361,361);
non_prepay_table(end,1:end)=mtg_mth_pmt;

for i=360:-1:1
    for j=1:i
        if i==1
            non_prepay_table(i,j)=0.5*(non_prepay_table(i+1,j)+non_prepay_table(i+1,j+1))*exp(-shortTree(i,j)*dt);
        else
            non_prepay_table(i,j)=0.5*(non_prepay_table(i+1,j)+non_prepay_table(i+1,j+1))*exp(-shortTree(i,j)*dt)+mtg_mth_pmt;
        end        
    end
end

%% optimal tree with refinance
% no prepay
non_prepay=nan(361,361);
non_prepay(end,1:end)=mtg_mth_pmt;

for i=360:-1:1
    for j=1:i
        non_prepay(i,j)=0.5*(non_prepay(i+1,j)+non_prepay(i+1,j+1))*exp(-shortTree(i,j)*dt);       
    end
end

optimal_Tree=nan(361,361);
optimal_Tree(end,1:end)=mtg_mth_pmt;
continuation_tree=nan(361,361);                   % continuation tree is the same as the continuation value in the excel file
continuation_tree(end,1:end)=mtg_mth_pmt;
for i=360:-1:1
    for j=i:-1:1
        if i==360
            continuation_tree(i,j)=0.5*(continuation_tree(i+1,j)+continuation_tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(continuation_tree(i,j),fixed_table.oustanding_principal(i)) + mtg_mth_pmt;   %compare and obtain the smallest value
        elseif i==1
            continuation_tree(i,j)=0.5*(optimal_Tree(i+1,j)+optimal_Tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(continuation_tree(i,j),fixed_table.oustanding_principal(i));
        else
            continuation_tree(i,j)=0.5*(optimal_Tree(i+1,j)+optimal_Tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(continuation_tree(i,j),fixed_table.oustanding_principal(i)) + mtg_mth_pmt;
        end
    end
end






    











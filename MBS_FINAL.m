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
observedData.price=1./((1+observedData.yield/2).^(2*m));
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
%-----------------building mortage with no prepayment----------------------
mtg=100000000;                                                %mortage 
mtg_r=0.07;                                                   %mortage yearly rate
mtg_mth_r=mtg_r/12;                                           %mortage monthly rate
mtg_mth_pmt=(mtg*(mtg_mth_r))/(1-(1/(1+mtg_mth_r)^360));      %mortage monthly payment

% new_m for vector of maturities and also new table
new_m=(dt:dt:T)';
new_m(2:361)=new_m;                                           %move one step backwards
new_m(1)=0;                                                   %the first time spot is t=0

% creating a fixed table
fixed_table=dataset(new_m);
fixed_table.fixed_monthly_payment=mtg_mth_pmt*ones(361,1);
fixed_table.outstanding_principal=nan(361,1);
fixed_table.int_pmt=nan(361,1);
fixed_table.principal_pmt=nan(361,1);
fixed_table.outstanding_principal(1)=mtg;

for i=2:length(fixed_table)
    fixed_table.int_pmt(i)=mtg_mth_r * fixed_table.outstanding_principal(i-1);
    fixed_table.principal_pmt(i)=mtg_mth_pmt - fixed_table.int_pmt(i);
    fixed_table.outstanding_principal(i)=fixed_table.outstanding_principal(i-1) - fixed_table.principal_pmt(i);
end

rep_fix=fixed_table.outstanding_principal(:);                %we only need to use the outstanding principal for the parts after this

%% using function built for non prepayment and prepayment
% initialize the input values we need
mortgageRate = 0.07;
principal = 100000000;
frequency=12;                                                % monthly payment
n=30;                                                        % time span is 30 years

multiplier=0.4:0.025:2;

for i =1:length(multiplier)
    mortgageValueWithRefinancing(i,:)=prepay_mortage(principal,mortgageRate,n,frequency,shortTree.*multiplier(i),rep_fix);          
    mortgageValueWithOutRefinancing(i,:)=non_prepay_mortage(principal,mortgageRate,n,frequency,shortTree.*multiplier(i));
end

%% plotting part
subplot(3,2,1);
plot(yieldData.m, yieldData.y,'*',observedData.m,observedData.yield,'LineWidth',2)
title('Fitted yield curve')
grid on
subplot(3,2,2);
grid on
plot(yieldData.m,yieldData.vol,'*',observedData.m,observedData.volatility,'LineWidth',2) 
title('Fitted volatility curve');
grid on

subplot(3,2,3)
plot(multiplier*shortTree(1,1),mortgageValueWithRefinancing(:,1),multiplier*shortTree(1,1),mortgageValueWithOutRefinancing(:,1), "LineWidth",2);
legend("Mtge value with prepayment option","Mtge value without prepayment option");
title("Value");
grid on

subplot(3,2,4)
plot(multiplier*shortTree(1,1),mortgageValueWithRefinancing(:,2),multiplier*shortTree(1,1),mortgageValueWithOutRefinancing(:,2), "LineWidth",2);
legend("Mtge duration with prepayment option","Mtge duration without prepayment option")
title("Duration");
grid on

subplot(3,2,5)
plot(multiplier*shortTree(1,1),mortgageValueWithRefinancing(:,3),multiplier*shortTree(1,1),mortgageValueWithOutRefinancing(:,3),'LineWidth',2);
legend("Mtge convexity with prepayment option","Mtge convexity without prepayment option")
title("Convexity");
grid on

fprintf("\n\n\n%s %23.2f",'the fixed monthly payment is: ',mtg_mth_pmt);

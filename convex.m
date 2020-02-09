clc;clear;close all;

N=4;
sigma=1;
alpha=[5;6;7;8];
pmax=[0.1 1 10];

mu = mu1(alpha,sigma, pmax) % find mu for the first set of alpha


alpha=[1;2;3;4];
mu2 = mu1(alpha,sigma, pmax) % find mu for the second set of alpha

%Comparing
figure(1)
hold on
plot(pmax,mu)
plot(pmax,mu2)
legend('\alpha=[5;6;7;8]','\alpha=[1;2,;3,;4]')
xlabel('pmax')
ylabel('mu')
hold off



function mu = mu1(alpha,sigma, pmax) %calculating mu
A=10; % First point for bisection
B=0; % Second point for bisection

for i=1:length(pmax)
    mu(i)=bisection1(pmax(i),alpha,sigma,A,B); % find mu for each pmax by bisection method
end
end
%%

function mu=bisection1(pmax,alpha,sigma,A,B)
mu=(A+B)/2; % find midle point
while abs(pmax-pn(mu,alpha,sigma))>=0.0000001% repeat the cycle body while difference between Pmax and sum of pn bigger then 0.000001  
    if ((pmax -pn(mu,alpha,sigma)).*(pmax-pn(A,alpha,sigma)))<0 % Check the condition if f(C)*f(A) < 0
        B=mu; % then B = C
    elseif ((pmax -pn(mu,alpha,sigma)).*(pmax-pn(A,alpha,sigma)))>=0% Check the condition if f(C)*f(A) > 0
        A=mu; % then A = C
    end
    mu=(A+B)/2; % find new C
end
end
function p = pn(mu,alpha,sigma)% find sum of  pn 
    p=0;
    for i = 1:4
    p = p + max(1/mu-sigma/alpha(i),0);
    end
end
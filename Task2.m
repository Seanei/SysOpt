clc;clear;close all;
%% Initial data
rng(1);
M=16;  % # of antennas
K=4;  % # of receivers
sigma=1;
H=1/sqrt(2)*(randn(M,K)+j*randn(M,K)); % generating of random channels 
... for K receivers with M antennas  
nu = [1 0.1 0.01 0.001 0.0001];
  % threshold of interference
SNR=0:5:30;   % Signal to noise ratio
pmax=10.^(SNR/10); 
%% Method 1
for s = 1:length(pmax)
    pk = pmax(s)/K;  % power for each receiver
    R = 0;         

for k = 1:K
    
    A= 10e6;  % First point for bisection
    B=10e-6;  % Second point for bisection
    lambda = bisection1(pk,H,k,A,B,M); % calc lambda via bisection method
    h = zeros(M,M);
    for i = 1:4
        if i ~= k
            h = H(:,i)*H(:,i)' + h; 
        end
    end
    vk = (inv(h + lambda * eye(M)))*H(:,k); % calc beamforming vector
    V(:,k,s)=vk; % collect beamforming vectiors
    Lam(k,s) = lambda;  % collect coresponding lambda
end
    for k = 1:K
        vk = V(:,k,s); 
        h=0;
        for i = 1:K
            if k ~=i
            h= (abs(H(:,k)'*V(:,i,s)))^2 + h; 
            end
        end
        gamma = ((abs(H(:,k)'*vk))^2)/(h+sigma); % calc of SINR 
        ... at the k-th receiver
        R = R + log2(1+gamma); % cell sum rate
    end
Rs(s)=R;
end
V2=V;
for n = 1:length(nu)
for s = 1:length(pmax)
 nun = nu(n)
[R, v2]=method2(Lam(:,s), nu(n), V2(:,:,s),K, H, sigma); % calc cell 
... sum rate and correspondig beamformers via method 2
V2(:,:,s)= v2; % save matrix of v
Rs2(s)=R;
end
Rsn(n,:)=Rs2; % collect cell sum rate for diff nu

end

%% Zero-forcing

for s = 1:length(pmax)
    pk = pmax(s)/K;  % power for each receiver
    V_zf = pinv(H'); 
    R=0;
    for k = 1:K
        vk = V_zf(:,k);
        vk = (vk/sqrt(vk'*vk))*sqrt(pk);
        h=0;
        for i = 1:K
            if k ~=i
            h= (abs(H(:,k)'*V_zf(:,i)))^2 + h; 
            end
        end
        gamma = ((abs(H(:,k)'*vk))^2)/(h+sigma); % calc of SINR 
        ... at the k-th receiver
        R = R + log2(1+gamma); % cell sum rate
    end
Rzf(s) = R;
end
%% Graphs for Task2 c) - i) iii) and ii)
%i) and iii)
figure(1)
hold on
plot(0:5:30,Rs2,'-r*','LineWidth',2)
plot(0:5:30,Rs,'b','LineWidth',2)
plot(0:5:30,Rzf,'g','LineWidth',2)

legend('Method 2','Method 1','Zero - Forcing')
xlabel('SNR [dB]')
ylabel('Cell sum rate')
grid on;
hold off

% ii)
figure(2)
hold on
plot(0:5:30,Rs,'b')
plot(0:5:30,Rsn(1,:),'-o')
plot(0:5:30,Rsn(2,:),'.-')
plot(0:5:30,Rsn(3,:),'-*')
plot(0:5:30,Rsn(4,:),'-o')
plot(0:5:30,Rsn(5,:),'-')

legend('Method1','Method2 nu=1', 'Method2 nu=10^-^1', 'Method2 nu=10^-^2', 'Method2 nu=10^-^3', 'Method2 nu=10^-^4')
xlabel('SNR [dB]')
ylabel('Cell sum rate')
grid on
hold off

function lambda=bisection1(pk,H,k,A,B,M)   % func. for bisection algorithm
lambda=(A+B)/2;     % choose lambda as middle of segment
h = zeros(M,M);      
    for i = 1:4
        if i ~= k
            H(:,i);
            h = H(:,i)*H(:,i)' + h; 
        end
    end
    vkl = (inv(h + lambda * eye(M)))*H(:,k);  % calc vk through Method 1 for lambda
    vkA = (inv(h + A * eye(M)))*H(:,k);       % calc vk through Method 1 for A
    vkB = (inv(h + B * eye(M)))*H(:,k);        %calc vk through Method 1 for B
    if er(B,vkB,pk)*er(A,vkA,pk)<0  % Check the condition if f(B)*f(A) < 0
        if er(lambda,vkl,pk)*er(A,vkA,pk)<0  % Check the condition if f(lambda)*f(A) < 0
            B=lambda;   % shift end of segment to lambda
        elseif er(lambda,vkl,pk)*er(A,vkA,pk) >0   % Check the condition if f(Blambda)*f(A) > 0
            A=lambda; % shift beggin of segment to lambda
        end
        if abs(er(lambda,vkl,pk))>=0.00001 % repeat the cycle body while
    ...difference between Pmax and sum of pn bigger then 0.00001   
            lambda = bisection1(pk,H,k,A,B,M);  
        end
    else 
        error('error')
    end
end
  
%% METHOD 2
function [R,V2]=method2(Lam, nu, V2,K, H, sigma)
    mu=ones(K,K);% initial mu is a set of 1
    lambda = Lam';
for k = 1:K

    for i = 1:4
        if i ~= k
            A = 10e6;%first point for bisection
            B = 10e-6;% second point for bisection          
            h = zeros(16,16);
            vk = V2(:,k);
            if abs(H(:,i)'*vk)^2>nu%if level of side channel bigger then target level of interference... 
               [mu_opt, vkm] = bisection2(A, B, H, lambda(k), k, i, nu, K, mu(k,:));%calling bisection method for assuming optimal vk
                V2(:,k) = vkm;%renew v matrix
            end
        end
    end
end
    R = 0;
    for k = 1:K% for each channel
            h=0;
        for i = 1:4% calculate sum of abs(hk'*vj)^2
            if k ~=i
            h= (abs(H(:,k)'*V2(:,i)))^2 + h;
            end
        end
        gamma = ((abs(H(:,k)'*V2(:,k)))^2)/(h+sigma);% calculate gamma
        R = R + log2(1+gamma);% calculate of sum rate
    end

end

function [mu_opt, vkm] = bisection2(A, B, H, lambda, k, j, nu, K, mu_opt)
    h_j = H(:,j); 
    h_k = H(:,k);
    M=size(H,1);

    mu_optA = mu_opt;
    mu_optA(j) = A;% set of mu for the first point
    mu_optB = mu_opt;
    mu_optB(j) = B;% set of mu for the second point
    vkm = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k;%Calculate vk for initial mu
    vkA = inv(hs(K, H, mu_optA,k)+ lambda*eye(M))*h_k;% Calculate vk for first point mu
    vkB = inv(hs(K, H, mu_optB,k)+ lambda*eye(M))*h_k;%Calculate vk for second point mu
    if err2(h_j, vkA,nu)*err2(h_j, vkB, nu)<=0%if function of mu crosses zero...
        mu_j = (A+B)/2;% middle between A and B
        mu_opt(j) = mu_j;% set of mu for middle point
        vkm = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k;% calculate vk for middle point set of mu
        if err2(h_j, vkA, nu)*err2(h_j, vkm, nu)<0% Condition of bisection 
            B = mu_j;
        else err2(h_j, vkA, nu)*err2(h_j, vkm, nu)>0
            A = mu_j;
        end
        if abs(err2(h_j, vkm,nu))> 0.00001% if error bigger then 0.00001...
            [mu_opt, vkm] = bisection2(A, B, H, lambda, k, j, nu, K, mu_opt);% call this function again
        end
    else %if function of mu doesn't cross zero...
        error('Error')
    end
    
    
end


function h = hs(K, H, mu,k)% function for calculating sum of muj*hj*hj'
h = zeros(size(H,1));
for j = 1: K
    if j ~= k
        h =h + mu(j) .* (H(:,j)*H(:,j)');
    end  
end
end

function err = er(lambda,vk, pk) %function for calculation error method 1

    err = lambda.*(vk'*vk - pk); 
end

function err = err2(h_j,vk, nu) %function for calculation error method 2
err = abs(h_j'*vk)^2-nu;
end
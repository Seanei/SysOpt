clc;clear;close all;
rng(5);
M=16;
K=4;
sigma=1;

H=1/sqrt(2)*(randn(K,M)+j*randn(K,M)); % random channel
nu=10^-2;
%%

for s = 0:1:6
    pk = 10^(s*5/10)/K;
    R = 0;

for k = 1:K
    
    A= 10e6;
    B=10e-6;
    [vk, lambda] = bisection1(pk,H,k,A,B,M);
    V(:,k)=vk;
    Lam(k) = lambda
end
    for k = 1:K
        h=0;
        for i = 1:4
            if k ~=i
            h= (abs(H(k,:)*V(:,i)))^2 + h;
            end
        end
        gamma = ((abs(H(k,:)*vk))^2)/(h+sigma);
        R = R + log2(1+gamma);
    end

Rs(s+1)=R;
end

for s = 0:6
    
Rs2=method2(Lam, nu, V,K, H, sigma,s);
end




%%
figure(1)
plot( 0:5:30,Rs)
figure(2)
 plot( 0:5:30,Rs2)  
function [vkl, lambda]=bisection1(pk,H,k,A,B,M)

lambda=(A+B)/2;
h = zeros(M,M);
    for i = 1:4
        if i ~= k
            h = H(i,:)'*H(i,:) + h;
        end
    end

    vkl = (inv(h + lambda * eye(M)))*H(k,:)';
    vkA = (inv(h + A * eye(M)))*H(k,:)';
    vkB = (inv(h + B * eye(M)))*H(k,:)';
    if er(B,vkB,pk).*er(A,vkA,pk)<0
        if er(lambda,vkl,pk).*er(A,vkA,pk)<0
            B=lambda;
        elseif er(lambda,vkl,pk).*er(A,vkA,pk) >0
            A=lambda;
        end
    
        if abs(er(lambda,vkl,pk))>=0.00001
            [vkl,lambda] = bisection1(pk,H,k,A,B,M);
        end
    else 
        fprintf('Error %d');
    end
end
  

function Rs=method2(Lam, nu, V,K, H, sigma,s);
for k = 1:K
    vk = V(:,k);
    M = size(H,2);
    lambda = Lam(k);
    mu=ones(K,K);
    for i = 1:4
        if i ~= k
            A = 10e6;
            B = 10e-6;
            mu_opt = bisection2(A, B, H, lambda, k, i, nu, K, mu(k,:))
            mu
            mu(k,:) = mu_opt;
            
           
        end
        vk = inv(hs(K, H, mu(k,:),k)+ lambda*eye(M))*H(k,:)'
        V(:,k)=vk;
    end

   
end

    for k = 1:K
        h=0;
        R = 0;
        for i = 1:4
            if k ~=i
            h= (abs(H(k,:)*V(:,i)))^2 + h;
            end
        end
        
        
        gamma = ((abs(H(k,:)*vk))^2)/(h+sigma)
        R = R + log2(1+gamma);
    end
Rs(s+1)=R;

end

function mu_opt = bisection2(A, B, H, lambda, k, j, nu, K, mu_opt)
    h_j = H(j,:);
    h_k = H(k,:);
    M=size(H,2);
    
    mu_opt(j) = A;
    vkA = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k';
    mu_opt(j) = B;
    vkB = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k';
    
    if err2(h_j, vkA,nu)*err2(h_j, vkB, nu)<0
        mu_j = (A+B)/2;
        mu_opt(j) = mu_j;
        vkm = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k';
        if err2(h_j, vkA, nu)*err2(h_j, vkm, nu)<0
            B = mu_j;
        else err2(h_j, vkA, nu)*err2(h_j, vkm, nu)>0
            A = mu_j;
        end
        if abs(err2(h_j, vkm,nu))> 0.01
            mu_opt = bisection2(A, B, H, lambda, k, j, nu, K, mu_opt);
        end
    else
        error('Error')
    end
    
    
end


function h = hs(K, H, mu,k)
h = zeros(size(H,2));
for j = 1: K
    if j ~= k
        h =h + mu(j) .* (H(j,:)'*H(j,:));
    end  
end

end
function err = er(lambda,vk, pk);

    err = lambda.*(vk'*vk - pk);
end

function err = err2(h_j,vk, nu)
err = abs(h_j*vk)^2-nu;
end
clc;clear;close all;
rng(1);
M=16;
K=4;
sigma=1;

H=1/sqrt(2)*(randn(M,K)+j*randn(M,K)); % random channel
nu=10e0;
%%

for s = 0:1:6
    pk = 10^(s*5/10)/K;
    R = 0;

for k = 1:K
    
    A= 10e6;
    B=10e-6;
    lambda = bisection1(pk,H,k,A,B,M);
    
    h = zeros(M,M);
    for i = 1:4
        if i ~= k
            h
            h = H(:,i)*H(:,i)' + h;
        end
    end
    h
    vk = (inv(h + lambda * eye(M)))*H(:,k);
    V(:,k)=vk;
    Lam(k,s+1) = lambda
end
    for k = 1:K
        vk = V(:,k)
        h=0;
        for i = 1:4
            if k ~=i
            h= (abs(H(:,k)'*V(:,i)))^2 + h;
            end
        end
        (abs(H(:,k)'*vk))
        gamma = ((abs(H(:,k)'*vk))^2)/(h+sigma);
        R = R + log2(1+gamma);
    end

Rs(s+1)=R;
end
for s = 1:7
V2=V;

[R, V2]=method2(Lam(:,s), nu, V2,K, H, sigma,s)
Rs2(s)=R
end




%%
figure(1)
hold on
plot(0:5:30,Rs2,'r')
plot(0:5:30,Rs,'b')
legend('2','1')
hold off

function lambda=bisection1(pk,H,k,A,B,M)

lambda=(A+B)/2;
h = zeros(M,M);
    for i = 1:4
        if i ~= k
            H(:,i);
            h = H(:,i)*H(:,i)' + h;
        end
    end
    vkl = (inv(h + lambda * eye(M)))*H(:,k);
    vkA = (inv(h + A * eye(M)))*H(:,k);
    vkB = (inv(h + B * eye(M)))*H(:,k);
    if er(B,vkB,pk).*er(A,vkA,pk)<0
        if er(lambda,vkl,pk).*er(A,vkA,pk)<0
            B=lambda;
        elseif er(lambda,vkl,pk).*er(A,vkA,pk) >0
            A=lambda;
        end
    
        if abs(er(lambda,vkl,pk))>=0.00001
            lambda = bisection1(pk,H,k,A,B,M);
        end
    else 
        error('error')
    end
end
  
%% METHOD 2
function [R,V2]=method2(Lam, nu, V,K, H, sigma,s);
    mu=ones(K,K);
    V2=V
    lambda = Lam'
for k = 1:K

    
    for i = 1:4
        if i ~= k
            A = 10e6;
            B = 10e-6;
            t=0
            
            h = zeros(16,16);
        for l = 1:4
            if l ~= k
                h = H(:,l)*H(:,l)' + h;
            end
        end
            vk = V(:,k);
            if abs(H(:,i)'*vk)^2>nu
                
                 lambda
                [mu_opt, vkm] = bisection2(A, B, H, lambda(k), k, i, nu, K, mu(k,:));
                V2(:,k) = vkm;
                mu(k,:) = mu_opt;
            
            end
        end

    end

   
end
    R = 0;
    h=0;
    for k = 1:K
        for i = 1:4
            if k ~=i

            (abs(H(:,k)'*V2(:,i)))^2
            h= (abs(H(:,k)'*V2(:,i)))^2 + h;
            end
        end
        
        gamma = ((abs(H(:,k)'*V2(:,k)))^2)/(h+sigma);
        R = R + log2(1+gamma)
    end
    R


end

function [mu_opt, vkm] = bisection2(A, B, H, lambda, k, j, nu, K, mu_opt)
    h_j = H(:,j);
    h_k = H(:,k);
    M=size(H,1);

    mu_optA = mu_opt;
    mu_optB = mu_opt;
    vkm = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k;
    mu_optA(j) = A;
    
    

    vkA = inv(hs(K, H, mu_optA,k)+ lambda*eye(M))*h_k;
    mu_optB(j) = B;
    vkB = inv(hs(K, H, mu_optB,k)+ lambda*eye(M))*h_k;

    if err2(h_j, vkA,nu)*err2(h_j, vkB, nu)<=0
        mu_j = (A+B)/2;
        mu_opt(j) = mu_j;
        vkm = inv(hs(K, H, mu_opt,k)+ lambda*eye(M))*h_k;
        if err2(h_j, vkA, nu)*err2(h_j, vkm, nu)<0
            B = mu_j;
        else err2(h_j, vkA, nu)*err2(h_j, vkm, nu)>0
            A = mu_j;
        end
        if abs(err2(h_j, vkm,nu))> 0.00001
            [mu_opt, vkm] = bisection2(A, B, H, lambda, k, j, nu, K, mu_opt);
        end
    else
        error('Error')
    end
    
    
end


function h = hs(K, H, mu,k)
h = zeros(size(H,1));

for j = 1: K
    if j ~= k
        
        h =h + mu(j) .* (H(:,j)*H(:,j)');
    end  
end

end
function err = er(lambda,vk, pk);

    err = lambda.*(vk'*vk - pk);
end

function err = err2(h_j,vk, nu)
err = abs(h_j'*vk)^2-nu;
end
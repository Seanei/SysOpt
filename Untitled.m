clc;clear;close all;
rng(2);
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
    V(:,k)=vk
    h=0;
end
    for k = 1:K
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
%%
plot( 0:5:30,Rs)
    
function [vkl, lambda]=bisection1(pk,H,k,A,B,M)

lambda=(A+B)/2;
h = zeros(M,M);
    for i = 1:4
        if i ~= k
            h = H(i,:)'*H(i,:) + h;
        end
    end

    vkl = (inv(h + lambda * eye(M)))*H(k,:)';
    vkA = (inv(h + A * eye(M)))*H(k,:)'
    vkB = (inv(h + B * eye(M)))*H(k,:)'
    if er(B,vkB,pk).*er(A,vkA,pk)<0
        fprintf('OK');
        if er(lambda,vkl,pk).*er(A,vkA,pk)<0
            B=lambda
        elseif er(lambda,vkl,pk).*er(A,vkA,pk) >0
            A=lambda
        end
    
        if abs(er(lambda,vkl,pk))>=0.00001
            [vkl,lambda] = bisection1(pk,H,k,A,B,M);
        end
    %else 
        %fprintf('Error %d');
    end
end
  
% 
% function vk=bisection2(pk,H,k,A,B,M, lambda, nu);
% 
% mu=(A+B)/2;
% h = zeros(M,M);
%     for i = 1:4
%         if i ~= k
%         h = H(i,:)'*H(i,:) + h;
%         end
%     end
% 
% vk = (inv(h + lambda * eye(M)))*H(k,:)';
% 
%     if er2(mu, h,vk,nu).*er2(mu, h,vk,nu)<0
%         B=lambda;
%     elseif er2(mu, h,vk,nu).*er2(mu, h,vk,nu) >=0
%         A=lambda;
%     end
%     
%     if abs(er2(mu, h,vk,nu))>=0.0000001;
%         vk = bisection1(pk,H,k,A,B,M);
%     end
% %     pn(mu);
% end
% 
% 
% function err = er2(mu, h,vk,nu)
% 
% 
% 
% end
function err = er(lambda,vk, pk);

    err = lambda.*(vk'*vk - pk)
end


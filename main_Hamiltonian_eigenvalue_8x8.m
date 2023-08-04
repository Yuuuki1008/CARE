clc;
clear;
rng(1)
format long e
% create H_tilde---------------------------------------------------
a_11 = [0 0.1; -0.1 0];
a_21 = zeros(2);
a_12 = zeros(2);
a_22 = [10 5; -5 10];
A = [a_11 a_12; a_21 a_22];
B = randi(5,4);
G = B*B'*10^(-2);
H = G;
H_t = [A G; H -A'];
%-------------------------------------------------------------------
%Compute eigenvalues
Lambda = eig(H_t);
figure(1)
plot(Lambda,'o')
title('Plot of eigenvalues')
%Choose eigenvalues
a = real(Lambda(1));
b = real(Lambda(7));
h = imag(Lambda(1));
l = imag(Lambda(7));

%find the case of gamma formula-------------------
C = (l^2-h^2)/(b-a);
test = (l^2/b)-a;
test_left = -(b-a)+2*h^2/a;
test_right = (b-a)+2*l^2/b;
%Compute gamma
if C <= test
    if C < test_left
        gamma = sqrt(a^2+h^2);
    elseif C > test_right
        gamma = sqrt(b^2+l^2);
    else
        gamma = sqrt(a*b+(a*l^2-b*h^2)/(b-a));
    end
else
    gamma = sqrt(b^2+l^2);
end
gamma_opt=gamma;
%-------------------------------------------------
%Compute A0 G0 H0 of SDA
A_gamma = A-gamma*eye(4);
W_gamma = A_gamma'+H*inv(A_gamma)*G;

A_prime = eye(4)+2*gamma*(inv(W_gamma))';
G_prime = -2*gamma*(inv(W_gamma))'*G*(inv(A_gamma))';
H_prime = -2*gamma*inv(W_gamma)*H*inv(A_gamma);
%---------------------------------------------------
% Slove SDA, output the step and solution X.
[it1,X] = SDA(A_prime,G_prime,H_prime)
rsdl = norm(-X*G*X-A'*X-X*A+H)
%---------------------------------------------------
% Change gamma and repeat the process.
% inital gamma = 0.01
gamma_step = 0.01;
it_sum = zeros(1,1000);
rsd_sum = zeros(1,1000);
gamma_sum =zeros(1,1000);
for k=1:1000
    gamma = 0.01 + gamma_step * k;
    A_gamma = A-gamma*eye(4);
    W_gamma = A_gamma'+H*inv(A_gamma)*G;
    A_prime = eye(4)+2*gamma*(inv(W_gamma))';
    G_prime = -2*gamma*(inv(W_gamma))'*G*(inv(A_gamma))';
    H_prime = -2*gamma*inv(W_gamma)*H*inv(A_gamma);
    [it,X] = SDA(A_prime,G_prime,H_prime);
    rsd = norm(-X*G*X-A'*X-X*A+H);
    it_sum(k) =it;
    rsd_sum(k) = rsd;
    gamma_sum(k) = gamma;
end
figure(2)
plot(gamma_sum,it_sum,'-')
hold on 
plot(gamma_opt,it1,'o')
hold off
xlabel('\gamma')
ylabel('iteration')
figure(3)
plot(gamma_sum,rsd_sum,'-')
hold on 
plot(gamma_opt,rsdl,'o')
hold off
xlabel('\gamma')
ylabel('residual')

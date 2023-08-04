clc;
clear;
rng(1)
format long e
% create H_tilde---------------------------------------------------
a_11 = [0 0.1; -0.1 0];
a_22 = [0 2; -2 0];
a_33 = [5 3; -3 5];
a_44 = [10 5; -5 10];
A = [a_11 zeros(2) zeros(2) zeros(2);
    zeros(2) a_22 zeros(2) zeros(2);
    zeros(2) zeros(2) a_33 zeros(2);
    zeros(2) zeros(2) zeros(2) a_44];
B = randi(5,8);
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
a = real(Lambda(3));
h = imag(Lambda(3));
x = real(Lambda(13));
y = imag(Lambda(13));
m = (y-h)/(x-a);
b = real(Lambda(15));
l = m*(b-a)+h;

%find the case of gamma formula-------------------
C = (l^2-h^2)/(b-a);
test = (l^2/b)-a;
test_left = -(b-a)+2*(h^2/a);
test_right = (b-a)+2*(l^2/b);

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

gamma_opt= gamma;
%-------------------------------------------------
%Compute A0 G0 H0 of SDA
A_gamma = A-gamma*eye(8);
W_gamma = A_gamma'+H*inv(A_gamma)*G;

A_prime = eye(8)+2*gamma*(inv(W_gamma))';
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
gamma_sum = zeros(1,1000);
for k=1:2000
    gamma = 0.01 + gamma_step * k;
    A_gamma = A-gamma*eye(8);
    W_gamma = A_gamma'+H*inv(A_gamma)*G;
    A_prime = eye(8)+2*gamma*(inv(W_gamma))';
    G_prime = -2*gamma*(inv(W_gamma))'*G*(inv(A_gamma))';
    H_prime = -2*gamma*inv(W_gamma)*H*inv(A_gamma);
    [it,X] = SDA(A_prime,G_prime,H_prime);
    rsd = norm(-X*G*X-A'*X-X*A+H);
    it_sum(k) =it;
    rsd_sum(k) = rsd;
    gamma_sum(k) = gamma;
end
%Compare with rectangle-------------------
test_rec = a*(b-a)/2;
if h^2 <= test_rec
    gamma_rec = sqrt(a*b-h^2);
else 
    gamma_rec = sqrt(a^2+h^2);
end
gamma = gamma_rec;
A_gamma = A-gamma*eye(8);
W_gamma = A_gamma'+H*inv(A_gamma)*G;

A_prime = eye(8)+2*gamma*(inv(W_gamma))';
G_prime = -2*gamma*(inv(W_gamma))'*G*(inv(A_gamma))';
H_prime = -2*gamma*inv(W_gamma)*H*inv(A_gamma);
% Slove SDA, output the step and solution X.
[it_rec,X_rec] = SDA(A_prime,G_prime,H_prime)
rsd_rec = norm(-X*G*X-A'*X-X*A+H)

% Plot figure of iteration and residual.
figure(2)
plot(gamma_sum,it_sum,'-')
hold on 
plot(gamma_opt,it1,'bo')
plot(gamma_rec,it_rec,'ro')
hold off
xlabel('\gamma')
ylabel('iteration')
legend('iteration','Trapezoidal','rectangle')
figure(3)
plot(gamma_sum,rsd_sum,'-')
hold on 
plot(gamma_opt,rsdl,'bo')
plot(gamma_rec,rsd_rec,'ro')
legend('residual','Trapezoidal','rectangle')
hold off
xlabel('\gamma')
ylabel('Residual')

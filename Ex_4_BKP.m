clearvars; close all; clc

% Define RK4 coefficients
b11 = 1;
b21 = 0.5;
b22 = b21;
a1 = 1;
a2 = 1;

% Request number 3
alfa = linspace(pi, 0.000001,1000);

F_RK2 = @(h,a) eye(size(A(a))) + a2*h*A(a)*(b21+b22) + h^2*a2*b22*b11*A(a)^2;
F_RK4 = @(h,a) (A(a)*h)^4/24 + (A(a)*h)^3/6 + (A(a)*h)^2/2 + A(a)*h + eye(2);

h = zeros(length(alfa),1);
lambdas = zeros(length(alfa),1);

% Compute the stability region
for i=1:length(alfa)
    try
        T = fzero(@(x) max(abs(eig(F_RK2(x,alfa(i))))) - 1,3);
        h(i)=T; 
    catch
        
    end
    
    lambdas(i,1:2) = eig(A(alfa(i)),'vector');
end

hlambda = h.*lambdas;
plot(hlambda)
grid on
axis equal

fill(real(hlambda),imag(hlambda), [0.8 0.8 0.8]) ;
hold on
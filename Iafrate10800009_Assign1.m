%% Ex 1
clearvars; close all; clc

f = @(x) cos(x) - x;
fplot(f), grid on;
a = 0;
b = 1.5;
Tol = 1e-8;

tic
[root,i] = bisection(f,a,b,Tol);
t_BI = toc();

tic
[root_1,i_1] = secant(f,a,b,Tol);
t_SEC = toc();

tic
[root_2, i_2] = regulafalsi(f,a,b,Tol);
t_RF=toc();

%% Ex 2 todo: how do we compute the accuracy??
clearvars; close all; clc

x = linspace(-2,2,100);
y = x;

f = @(x) [x(1).^2 - x(1) - x(2); x(1).^2/16 + x(2).^2 - 1];


J = @(x) [2*x(1)-1, -1; x(1)/8, 2*x(2)];           % Analytical Jacobian
x0 = [1;1];
Tol = 1e-8;

tic
[roots,fevals] = newton_analytical(f,x0,J,Tol);
t_analytical = toc();

tic
[roots1,fevals1] = newton_FD(f,x0,Tol);
t_FD = toc();

% Compute exact roots
syms x y real
[solx,soly] = solve(x^2-x-y==0,x^2/16+y^2-1);

% Compare the two solutions

%% Ex 3
clearvars; close all; clc

f = @(t,x) x - t^2 + 1;
x_an = @(t) t.^2 + 2*t +1 - 0.5*exp(t);
tspan = [0 2];
x0 = 0.5;

h = [0.5 0.2 0.05 0.01];

% Memory allocation
solRK2 = zeros(size(h));
timeRK2 = solRK2;

figure (1)
hold on

% RK2 errors and timings
for i=1:length(h)
    [sol,~,times] = RK2(f,tspan,x0,h(i));
    solRK2(i) = sol(:,end);
    plot(times,sol)
    g = @() RK2(f,tspan,x0,h(i));
    timeRK2(i) = timeit(g);
end
plot(0:0.5:2,x_an(0:0.5:2))
grid on
legend('h=0.5','h=0.2','h=0.05','h=0.01','analytical')

% RK4 errors and timings
solRK4 = zeros(size(h));
timeRK4 = solRK4;

figure(2)
hold on
for i = 1:length(h)

    [sol,fevals(i),times] = RK4(f,tspan,x0,h(i));
    solRK4(i) = sol(:,end);
    plot(times,sol)
    g = @() RK4(f,tspan,x0,h(i));
    timeRK4(i) = timeit(g);
end
plot(0:0.01:2,x_an(0:0.01:2))
grid on
legend('h=0.5','h=0.2','h=0.05','h=0.01','analytical')


errorRK4 = abs(solRK4-x_an(2));

% Plotting
figure(3)
errorRK2 = abs(solRK2-x_an(2));
fig_RK2=loglog(timeRK2,errorRK2);
fig_RK2.LineWidth=2;
grid on
xlabel('Integration time [s]')
ylabel('Solution error')

figure(4)
fig_RK4=loglog(timeRK4,errorRK4);
fig_RK4.LineWidth=2;
xlabel('Integration time [s]')
ylabel('Solution error')
grid on

%% Ex 4 
clearvars; close all; clc

% Request number 3
alfa = linspace(pi, 0,100);

algorithm = 4;             %  Algorithm selector for ffwd fcn. (2 = RK2)

% Determine the stepsize to have discrete eigenvalues on the unitary circle
h = zeros(length(alfa),1);
problems = [];

for i = 1:length(alfa)
    a = alfa(i);
    try
        h(i) = fzero(@(x) max(abs(eig(ffwd(A(a),x,algorithm))))-1,5);
        
    catch mess
        warning(mess.message)
    end
end

% Compute continous system eigenvalues
lambdas = zeros(length(alfa),1);

for i=1:length(alfa)
 lambdas(i) = max(eig(A(alfa(i))));
end

hlambda = h.*lambdas;
plot(hlambda)
hold on
plot(real(hlambda),-imag(hlambda))
grid on
axis equal

fill(real(hlambda),imag(hlambda), [0.8 0.8 0.8]) ;
fill(real(hlambda),-imag(hlambda), [0.8 0.8 0.8]) ;
hold on

% Given that lambda = 1
h = [0.5 0.2 0.05 0.01];

plot(h,'*')

%% Ex 5
clearvars; close all; clc

x0 = [1;1];
tspan = [0 1];

alfa = linspace(pi, 0,50);
x_an = @(a) expm(A(a))*x0;        % Final time is t=1, no need to specify it
tol = [1e-3 1e-4 1e-5 1e-6];

LineSpec = {'LineWidth',2};

for i=1:length(alfa)
 lambdas(i) = max(eig(A(alfa(i))));
end

guess = [tol;
    0.01 0.01 0.01 0.01;
    0.5 0.5 0.5 0.5];


for i=1:3
    figure(i)
    hold on
    axis equal
    for j=1:length(tol)
        for k=1:length(alfa)
            a = alfa(k);
            h(k) = fzero(@(x) norm(wrapper(@(t,y) A(a)*y,tspan,x0,abs(x),i) ...
                    - x_an(a),inf) - tol(j),guess(i,j));
        end
        hl = abs(h).*lambdas;
        plot([real(hl) flip(real(hl))],...
            [imag(hl) -flip(imag(hl))],LineSpec{:})
        legend_entries{j}=sprintf("tol=%0.0e",tol(j));

    end
    grid on
    xlabel('$Re\{h\lambda\}$','Interpreter','latex')
    ylabel('$Im\{h\lambda\}$','Interpreter','latex')

    legend(legend_entries{:})
end
%% Ex 6
clearvars; close all; clc

% alfa angles vector
alfa = linspace(pi, 0, 100);

% select the F operator of BI2
algorithm = 5;

% Compute continous system eigenvalues
lambdas = zeros(length(alfa),1);

for i=1:length(alfa)
 lambdas(i) = max(eig(A(alfa(i))));
end

h = zeros(length(alfa),1);

% Req.2 plot the stability region of BI2 for θ = 0.4
th=0.4;

for i=1:length(alfa)
    a=alfa(i);
    
    try
        g = @(x) max(abs(eig(ffwd(A(a),x,algorithm,th))))-1;
        h(i) = fzero(g,10);
        
    catch mess
        warning(mess.message)
    end
    
end

figure()
hold on
LineSpecs = {'LineWidth',2};

hl = abs(h).*lambdas;
X=[real(hl); flip(real(hl))];
Y = [imag(hl); -flip(imag(hl))];
plot(X,Y,LineSpecs{:})

legend_entries = {'$BI2_{0.4}$'};
grid on
axis equal
axis padded
legend(legend_entries{:},'Interpreter','latex','Location','best')

xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')

% plot the stability region for θ = [0.1 0.3 0.7 0.9]
theta = [0.1 0.3 0.7 0.9];

figure()
hold on

for k = 1:length(theta)
    th=theta(k);
    
    for i=1:length(alfa)
        a=alfa(i);
        
        try
            g = @(x) max(abs(eig(ffwd(A(a),x,algorithm,th))))-1;
            h(i) = fzero(g,5);
        
        catch mess
            warning(mess.message)
        end

    end

    hl = abs(h).*lambdas;
    X=[real(hl); flip(real(hl))];
    Y = [imag(hl); -flip(imag(hl))];
    plot(X,Y,LineSpecs{:})
%     f = fill(X,Y,'r');
%     set(f,'facealpha',.5)
end

legend_entries = {'$BI2_{0.1}$','$BI2_{0.3}$','$BI2_{0.7}$','$BI2_{0.9}$'};
grid on
axis equal
axis padded
legend(legend_entries{:},'Interpreter','latex','Location','best')

xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')

%% Ex 7
clearvars; close all; clc
B = [-180.5 219.5; 179.5 -220.5];       % Stiff system

x0 = [1;1];
tspan = [0 5];
h = 0.1;

xt = @(t) expm(B.*t)*x0;

% RK4 integration (no bueno)
numerical = RK4(@(t,x) B*x,tspan,x0,h);

% BI2_0.4 integration (decent)
theta = 0.1;
F = ffwd(B,h,5,theta);

Nsteps = diff(tspan)/h;
xx=zeros(2,Nsteps);
analy = xx;

xx(:,1)=x0;

for k=1:Nsteps
    xx(:,k+1) = F*xx(:,k);
end

% Discussion
times = tspan(1):h:tspan(end);
for i=1:Nsteps+1
    analy(:,i)=xt(times(i));
end

% Analytical plot
figure(1)
plot(times,analy(1,:),times,analy(2,:))
hold on

% Numerical plots
plot(times,xx(1,:),times,xx(2,:))
grid on

% figure(2)
% plot(times,analy(1,:))
% hold on
% % Numerical plot
% plot(times,numerical(1,:))
% grid on

% Eigenvalues of the dynamical system (stiff) have very different values
hlambdas = h*eig(B);

figure()

hold on
[R_BI2,X,Y] = easystability(5,theta);
[~,figBI2] = contour(X,Y,-R_BI2,[-1 -1]);
figBI2.LineWidth=1.5;

R_RK4 = easystability(4);
[~,figRK4] = contour(X,Y,-R_RK4,[-1 -1],'k');
figRK4.LineWidth=1.5;

plot(real(hlambdas),imag(hlambdas),'*r')
grid on
axis padded

%% Functions
function [root,fevals] = bisection(f,a,b,Tol)

%Initialization
x = (a+b)/2;
fevals = 0;

while abs(b-a)/2>=Tol
    x = (a+b)/2;
    
    if f(x)*f(a)<0
        b = x;
    else
        a = x;
    end
    
    fevals = fevals + 1;
    
end
root = x;
end

function [root,fevals] = secant(f,x0,x1,Tol)
xk = x1;
xk_1 = x0;

fevals = 2;

while abs(xk-xk_1)>=Tol
    fk = f(xk);
    fk_1 = f(xk_1);
    
    xk1 = xk - (xk - xk_1)*fk/(fk - fk_1);
    fevals = fevals + 2;
    
    xk_1=xk;
    xk = xk1;
end
root = xk1;
end

function [root,fevals] = regulafalsi(f,a,b,Tol)

%Initialization
fevals = 0;

while abs(b-a)/2>=Tol
    fa = f(a);
    fb = f(b);
    x = b - (b-a)*fb/(fb-fa);
    fx = f(x);
    
    if fx>0
        b = x;
        
    else
        a = x;
        
    end
    
    fevals = fevals + 3;
    
end

root = x;
end

function [roots,fevals] = newton_analytical(fun,x0,J,Tol)
x = x0;
f = fun(x);
fevals = 0;
Nmax = 1000;

while abs(f)>=Tol % && fevals<Nmax
    f = fun(x);
    dX = -J(x)\f;
    x = x +dX;
    fevals = fevals+1;
end
roots = x;
end

function [roots,fevals] = newton_FD(fun,x0,Tol)
x = x0;
f = fun(x);
fevals = 0;
while abs(f)>=Tol
    f = fun(x);
    invJ = inv(FD_Jacobian(x,fun));     % HOW CAN I AVOID COMPUTING THE INVERSE?
    x = x - invJ*f;
    fevals = fevals+1;
end
roots = x;
end

function J = FD_Jacobian(X,f)

J = zeros(size(X,1));
F0 = f(X);
hj = max(sqrt(eps),sqrt(eps)*abs(X));

for i=1:size(X)
    EPS_i = zeros(size(X));
    EPS_i(i) = hj(i);
    Fi = f(X + EPS_i);
    J(:,i) = (Fi - F0)'/hj(i);
end

end

function [x,fevals,times] = RK1(f,tspan,x0,h)

x = x0;
times = tspan(1):h:tspan(end);
fevals = 0;
t = tspan(1);

if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

hs = diff(times);
for i=1:(length(times)-1)
    
    h = hs(i);
    fk = f(t,x);
    fevals = fevals + 1;
    x = x + h*fk;
    t = t + h;
end
end

function [x,fevals,times] = RK2(f,tspan,x0,h)

x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);

fevals = 0;

if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

hs = diff(times);

for i=1:(length(times)-1)
    xk = x(:,i);
    
    fk = f(t,xk);
    
    % Shrink the last timestep
    h = hs(i);
    
    % Predictor step
    xp = xk + h*fk;
    t = t + h;
    fp = f(t,xp);
    fevals = fevals + 1;
    
    % Corrector step
    x(:,i+1) = xk + h*(0.5*fk + 0.5*fp);
    
end

end

function [x,fevals,times] = RK4(f,tspan,x0,h)

x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);

fevals = 0;

if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

hs = diff(times);

for k=1:(length(times)-1)
    h = hs(k);
    
    xk = x(:,k);
    
    % Predictor steps
    % #1
    k1 = f(t,xk);
    
    % #2
    t2 = t + h/2;
    x2= xk + k1*h/2;
    k2 = f(t2,x2);
    
    % #3
    t3 = t + h/2;
    x3 = xk + k2*h/2;
    k3 = f(t3,x3);
    
    % #4
    t4 = t + h;
    x4 = xk + k3*h;
    k4 = f(t4,x4);
    
    % Corrector step
    x(:,k+1) = xk + (k1 + 2*k2 + 2*k3 + k4)*h/6;
    
    t = t+h;
    fevals = fevals+4;
    
end
end

function A_A = A(alfa)
A_A =[0 1; -1 2*cos(alfa)];
end

function F = ffwd(A,h,algor,th)

% Function that returns the F operator s.t. x_(k+1) = F*xk
% algor:
%   + 1 RK1
%   + 2 RK2
%   + 4 RK4

[~,n]=size(A);
I = eye(n);

switch algor
    case 1      % Forward euler (RK1)
        F = I + h*A;
        
    case 2      % Runge-Kutta (RK2)
        F = I + h*A + 0.5*(h*A)^2;
        
    case 4      % Runge-Kutta (RK4)
        F =  I + A*h + (A*h)^2/2 + (A*h)^3/6 + (A*h)^4/24;
        
    case 5      % θ-methods
        F = inv(I - A*(1-th)*h + 0.5*(A*(1-th)*h)^2)*...
            (I + A*th*h + 0.5*(A*th*h)^2);
        
end
end

function [Rhat,X,Y]=easystability(algorithm,th)

% INPUT:
%   algorithm       variable defining the type of algorithm to plot
%                       2: RK2
%                       4: RK4
%                       5: BI2_theta
[X,Y] = meshgrid(-5:0.01:5,-5:0.01:5);
hl = X+1i*Y;

switch algorithm
    case 2 % RK2
        R = 1 + hl + .5*hl.^2;
    case 4 % RK4
        R = 1 + hl + .5*hl.^2 + (1/6)*hl.^3 + (1/24)*hl.^4;
    case 5 % BI2_theta
        R = 1./(1 - (1-th).*hl+ 0.5*((1-th)*hl).^2).*...
            (1 + th*hl + 0.5*(th*hl).^2);
end

Rhat = abs(R);
end

function z_end=wrapper(f,tspan,x0,h,algor)
% Wrapper function, returning the state at the final integration step

switch algor
    case 1
        Z = RK1(f,tspan,x0,h);
    case 2
        Z = RK2(f,tspan,x0,h);
    case 3
        Z = RK4(f,tspan,x0,h);
end
z_end=Z(:,end);
end
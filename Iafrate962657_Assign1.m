% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Ex 1
clearvars; close all; clc

% Plot the function to estimate the zero interval
xx=linspace(-5,5,1000);
f = @(x) cos(x) - x;
plot(xx,f(xx),'LineWidth',2), grid on;
hold on
plot(xx,zeros(size(xx)),'--r','LineWidth',1.5)
xlabel('x','Interpreter','latex')
ylabel('f(x)','Interpreter','latex')

% zero interval from plot
a = 0.5;
b = 1.1;

% Digits accuracy
accuracy = 8;       
N = 10^-(accuracy+1);

% Root computations
tic
[root,i_0] = bisection(f,a,b,N);
t_BI = toc();

tic
[root1,i_1] = secant(f,a,b,N);
t_SEC = toc();

tic
[root2, i_2] = regulafalsi(f,a,b,N);
t_RF=toc();

%% Ex 2
clearvars; close all; clc

% Plot of the vector field
f = @(x) [x(1).^2 - x(1) - x(2); x(1).^2/16 + x(2).^2 - 1];
[X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
U = X.^2 - X - Y;
V =  X.^2/16 + Y.^2 - 1;
q=quiver(X,Y,U,V,2.5);
q.MaxHeadSize = 0.05;

labelopts = {'Interpreter','latex','FontSize',14};
titleopts = {'Interpreter','latex','FontSize',14};

title('$\vec{f}(\vec{x})$',titleopts{:})
xlabel('$x_1$',labelopts{:})
ylabel('$x_2$',labelopts{:})

% Initial guesses for the zeros, from the plot
x01 = [1.1;1];                                        
x02 = [-1.1;1];

J = @(x,unused) [2*x(1)-1, -1; x(1)/8, 2*x(2)];      % Analytical Jacobian

N = 5;      % Number of iterations

tic
jacobian_type = 1;          % 1: analytical jacobian
roots(:,1) = newton(f,x01,N,jacobian_type,J);
roots(:,2) = newton(f,x02,N,jacobian_type,J);
t_analytical = toc();

tic
jacobian_type = 2;          % 2: forward difference jacobian
rootsFD(:,1) = newton(f,x01,N,jacobian_type);
rootsFD(:,2) = newton(f,x02,N,jacobian_type);
t_FD = toc();

tic
jacobian_type = 3;          % 3: central difference jacobian
rootsCD(:,1) = newton(f,x01,N,jacobian_type);
rootsCD(:,2) = newton(f,x02,N,jacobian_type);
t_CD = toc();

% Compute exact roots
syms x y real
[solx,soly] = solve(x^2-x-y==0, x^2/16 +y^2-1==0);
SOL = double([solx';soly']);

% Compare the solutions to the exact one
error_an = vecnorm(roots-SOL);
error_FD = vecnorm(rootsFD-SOL);
error_CD = vecnorm(rootsCD-SOL);

%% Ex 3
clearvars; close all; clc

f = @(t,x) x - t^2 + 1;     % Analytical solution
x_an = @(t) t.^2 + 2*t +1 - 0.5*exp(t);
tspan = [0 2];
x01 = 0.5;

h = [0.5 0.2 0.05 0.01];

% Memory allocation
solRK2 = zeros(size(h));
timeRK2 = solRK2;

% Setup plots style
legend_options = {'Interpreter','latex','Location','best',...
    'FontSize',12};
label_options = {'Interpreter','latex','FontSize',14};

figure (1)
hold on

errorRK2 = zeros(size(h));
% RK2 computations
for i=1:length(h)
    [sol,~,times] = RK2(f,tspan,x01,h(i));
    solRK2(i) = sol(:,end);
    plot(times,sol,'--')
    errorRK2(i) = max(abs(sol-x_an(0:h(i):tspan(end))));
    g = @() RK2(f,tspan,x01,h(i));
    timeRK2(i) = timeit(g);
end

timesRK2 = tspan(1):0.5:tspan(2);

plot(timesRK2,x_an(timesRK2),'LineWidth',1.2)
grid on

legend_entries = {'h=0.5','h=0.2','h=0.05','h=0.01','analytical'};
legend(legend_entries{:},legend_options{:})
xlabel('$time [s]$',label_options{:})
ylabel('$x_{RK2}(t)$',label_options{:})

% RK4 errors and timings
solRK4 = zeros(size(h));
timeRK4 = solRK4;

figure(2)
hold on
fevals = zeros(size(h));
errorRK4 = zeros(size(h));

for i = 1:length(h)
    
    [sol,fevals(i),times] = RK4(f,tspan,x01,h(i));
    solRK4(i) = sol(:,end);
    plot(times,sol,'--')
    errorRK4(i) = max(abs(sol-x_an(0:h(i):tspan(end))));
    g = @() RK4(f,tspan,x01,h(i));
    timeRK4(i) = timeit(g);
end
timesRK4 = tspan(1):0.01:tspan(2);
plot(timesRK4,x_an(timesRK4),'LineWidth',1.2)
grid on
legend(legend_entries{:},legend_options{:})
xlabel('$time [s]$',label_options{:})
ylabel('$x_{RK4}(t)$',label_options{:})
    
% Plotting
figure(3)
fig_RK2=loglog(timeRK2,errorRK2);
fig_RK2.LineWidth=2;
grid on
xlabel('Integration time [s]',label_options{:})
ylabel('Solution error',label_options{:})

figure(4)
fig_RK4=loglog(timeRK4,errorRK4);
fig_RK4.LineWidth=2;
xlabel('Integration time [s]',label_options{:})
ylabel('Solution error',label_options{:})
grid on

%% Ex 4 
clearvars; close all; clc

% Request number 3
alfa = linspace(pi, 0,100);

algorithm = [2 4];             %  Algorithm selector for ffwd fcn. (2 = RK2)

% Determine the stepsize to have discrete eigenvalues on the unitary circle
h = zeros(length(alfa),1);
hs = [0.5 0.2 0.05 0.01];   % Specified timesteps of exercise 3

for k = 1:2
    figure(k)
    for i = 1:length(alfa)
        a = alfa(i);
        try
            h(i) = fzero(@(x) max(abs(eig(ffwd(A(a),x,algorithm(k)))))-1,5);
        catch mess
            warning(mess.message)
        end
    end

    % Compute continous system eigenvalues
    lambdas = zeros(length(alfa),1);

    for i=1:length(alfa)
     lambdas(i) = max(eig(A(alfa(i))));
    end

    hl = h.*lambdas;
    hold on
    
    X=[real(hl); flip(real(hl))];
    Y = [imag(hl); -flip(imag(hl))];
    
    fill(X,Y, [0.8 0.8 1]) ;

    plot(X,Y,'HandleVisibility','off')
    grid on
    axis equal
    axis padded
    
    
    xlabel('$Re\{h\lambda\}$','Interpreter','latex')
    ylabel('$Im\{h\lambda\}$','Interpreter','latex')

    % Given that lambda = 1
    
    plot(hs,zeros(size(hs)),'*r')
    
    legend('Stable Region','$h_i\lambda$','Interpreter','latex')
end

%% Ex 5
% (TAKES A COUPLE OF MINUTES TO RUN)

clearvars; close all; clc

% Initial conditions and timespan
x01 = [1;1];
tspan = [0 1];

alfa = linspace(pi, 0,100);       % angles vector of the eigenvalues
x_an = @(a) expm(A(a))*x01;       % analytical solution, at final time t=1

% tolerances vector
tol = [1e-3 1e-4 1e-5 1e-6];

LineSpec = {'LineWidth',2};

% Compute continous system eigenvalues
for i=1:length(alfa)
    lambdas(i) = max(eig(A(alfa(i))));
end

% Initial guesses for fzero
guess = [tol;                   % RK1 guesses
        0.01 0.01 0.01 0.01;    % RK2
        0.5 0.5 0.5 0.5];       % RK4

tic
 
for i=1:3       % three integration methods
    figure(i)
    hold on
    axis equal
    
    for j=1:length(tol)     % for each tolerance
        guessRK1 = tol(j);      % Initial guess for RK1
        options = optimset('TolX',1e-2*tol(j));
        
        for k=1:length(alfa)
            
            a = alfa(k);
            AA = A(a);
            
            if i==1
                h(k) = fzero(@(x) norm(wrapper(@(t,y) AA*y,tspan,x01,abs(x),i) ...
                    - x_an(a),inf)-tol(j),guessRK1,options);
                
                % guess update for next step
                guessRK1=h(k);
                
            else
                h(k) = fzero(@(x) norm(wrapper(@(t,y) AA*y,tspan,x01,abs(x),i) ...
                    - x_an(a),inf) - tol(j),guess(i,j));
            end
        end
        AA = A(pi);
        [~,fevals(i,j)] = wrapper(@(t,y) AA*y,tspan,x01,abs(h(1)),i);
            
        % Plotting
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

time=toc();

figure()
loglog(tol,fevals,'LineWidth',2)
grid on
xlabel('Tolerance','Interpreter','latex','FontSize',12)
ylabel('Function evaluations','Interpreter','latex','FontSize',12)
legenditems = {'RK1','RK2','RK4'};
legend(legenditems{:},'Interpreter','latex','FontSize',12)


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

legend(legend_entries{:},'Interpreter','latex','Location','best')

xlabel('$Re\{h\lambda\}$','Interpreter','latex')
ylabel('$Im\{h\lambda\}$','Interpreter','latex')

levelsopts = {'ShowText','on','LineStyle','--','HandleVisibility','off',...
    'LineWidth',1.5};
        
axis padded

[R_BI2_09,X,Y] = easystability(5,theta(4));
contour(X,Y,R_BI2_09,[0.8 1.2],...
                    'LineColor',[0.49,0.18,0.56],levelsopts{:});

[R_BI2_01,X,Y] = easystability(5,theta(1));
contour(X,Y,R_BI2_01,[0.8 1.2],...
                    'LineColor',[0.00,0.45,0.74],levelsopts{:});
axis([-5.5 5.5 -4 4])

%% Ex 7

clearvars; close all; clc
B = [-180.5 219.5; 179.5 -220.5];       % System dynamics (Stiff system)

x01 = [1;1];
tspan = [0 5];
h = 0.1;

% Analytical solution
xt = @(t) expm(B.*t)*x01;

% RK4 integration (no bueno)
numerical = RK4(@(t,x) B*x,tspan,x01,h);

% BI2_0.4 integration (decent)
theta = 0.1;
F = ffwd(B,h,5,theta);

Nsteps = diff(tspan)/h;
xx=zeros(2,Nsteps);
analy = xx;

xx(:,1)=x01;

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
plot(times,analy(1,:),times,analy(2,:),'LineWidth',1.5)
hold on

% BI2 plots
plot(times,xx(1,:),times,xx(2,:),'LineWidth',1.5)
grid on

legend_options = {'Interpreter','latex','Location','best',...
    'FontSize',12};
label_options = {'Interpreter','latex','FontSize',14};

legend_entries = {'$x_1 - analytical$','$x_2 - analytical$',...
                    '$x_1 - BI2_{0.4}$','$x_2 - BI2_{0.4}$'};
                
legend(legend_entries{:},legend_options{:})
xlabel('$time [s]$',label_options{:})
ylabel('$x_1\;;x_2$',label_options{:})

% RK4 plots
figure()
semilogy(times,[abs(numerical(1,:)); abs(numerical(2,:))],'LineWidth',2)
grid on
legend_entries = {'$x_1 - RK4$','$x_2 - RK4$'};

legend(legend_entries{:},legend_options{:})
xlabel('$time [s]$',label_options{:})
ylabel('$log(x_1)\;;log(x_2)$',label_options{:})

% Eigenvalues plots in the stability region
hlambdas = h*eig(B);

figure()

hold on
[R_BI2_09,X,Y] = easystability(5,theta);
[~,figBI2] = contour(X,Y,-R_BI2_09,[-1 -1]);
figBI2.LineWidth=1.5;

R_RK4 = easystability(4);
[~,figRK4] = contour(X,Y,-R_RK4,[-1 -1],'k');
figRK4.LineWidth=1.5;
fig_RK4.Fill='on';

plot(real(hlambdas),imag(hlambdas),'*r','LineWidth',2)
grid on
axis padded

legend_entries = {'$BI2_{0.1} stability$','$RK4 \; stability$',...
                    '$h\lambda_i$'};
legend(legend_entries{:},legend_options{:})

xlabel('$Re\{h\lambda\}$',label_options{:})
ylabel('$Im\{h\lambda\}$',label_options{:})

%% Functions
function [root,fevals] = bisection(f,a,b,Tol)
%BISECTION  compute the root of an equation using the *bisection* method
% 
% Author:
%       Davide Iafrate
% 
% INPUT:
%   f       function handle of the objective function
%   a       left extreme of the search interval
%   b       right extreme of the search interval
%   Tol     tolerance of the solution
% 
% OUTPUT:
%   root    estimated root of the function
%   fevals  number of function evaluations

%Initialization
fevals = 1;
fa = f(a);

while abs(b-a)/2>=Tol
    c = (a+b)/2;
    
    fc = f(c);
    
    if fc*fa<0
        b = c;
    else
        a = c;
        fa = fc;
    end
    
    fevals = fevals + 1;
    
end
root = c;
end

function [root,fevals] = secant(f,x0,x1,Tol)
%SECANT  compute the root of an equation using the *secant* method
% 
% Author:
%       Davide Iafrate
% 
% INPUT:
%   f       function handle of the objective function
%   x0       left extreme of the search interval
%   x1       right extreme of the search interval
%   Tol     tolerance of the solution
% 
% OUTPUT:
%   root    estimated root of the function
%   fevals  number of function evaluations

xk = x1;
xk_1 = x0;

fevals = 0;

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
%REGULAFALSI  compute the root of an equation using the *regula falsi*
%               method
% 
% Author:
%       Davide Iafrate
% 
% INPUT:
%   f       function handle of the objective function
%   a       left extreme of the search interval
%   b       right extreme of the search interval
%   Tol     tolerance of the solution
% 
% OUTPUT:
%   root    estimated root of the function
%   fevals  number of function evaluations

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

function [roots] = newton(fun,x0,Nmax,type,J)
%newton  compute the roots of function *fun* using Newton's algorithm, with
%           either the analytical Jacobian or its numerical approximation
%
% Author
%   Davide Iafrate
% 
% Inputs:
%   fun:  function handle representing the RHS of the system
%   x0 : initial guess
%   Nmax: number of iterations
%   jacobian_type:  scalar value to choose what Jacobian to use
%                      |-> 1: analytical
%                      |-> 2: Forward Difference numerical
%                      |-> 3: Central Difference numerical
%
% Outputs:
%   roots  :    roots of the system

x = x0;

if nargin<5
    JAC = {NaN,@FD_Jacobian,@CD_Jacobian};
elseif nargin==5
    JAC = {J,@FD_Jacobian,@CD_Jacobian};
else
    error("Too many arguments")
end

for i = 1:Nmax
    f = fun(x);
    Jac = JAC{type}(x,fun);
    DX = - Jac\f;
    x = x + DX;
end
roots = x;
end

function J = FD_Jacobian(X,f)
%FD_Jacobian  return the Jacobian of function f computed numerically using
%   the Forward Difference scheme
%
%
% Author
%   Davide Iafrate
% 
% Inputs:
%   X : [N,1] state of the system in which we want to evaluate the jacobian
%   f:  function handle representing the RHS of the system
%
% Outputs:
%   J[NxN]  :    System Jacobian in X

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

function J = CD_Jacobian(X,f)
%CD_Jacobian  return the Jacobian of function f computed numerically using
%   the Central Difference scheme
%
%
% Author
%   Davide Iafrate
% 
% Inputs:
%   X : [N,1] state of the system in which we want to evaluate the jacobian
%   f:  function handle representing the RHS of the system
%
% Outputs:
%   J[NxN]  :    System Jacobian in X

J = zeros(size(X,1));

% Compute the perturbation
hj = max(sqrt(eps),sqrt(eps)*abs(X));

for i=1:size(X)
    % Perturb the function in each state direction
    EPS_i = zeros(size(X));
    EPS_i(i) = hj(i);
    Fi = f(X + EPS_i);
    F_i = f(X - EPS_i);
    
    % Compute the i-th column of the Jacobian
    J(:,i) = (Fi - F_i)'/(2*hj(i));
end

end

function [x,fevals,times] = RK1(f,tspan,x0,h)
%RK1  integrate non-linear system dynamics using the Runge-Kutta 1 (Forward
%       Euler) integration
%
%
% Author
%   Davide Iafrate 
% 
% Inputs:
%   f : function handle of the RHS of the system dynamics, f(t,x)
%   tspan[1x2]  : vector of time interval extremes
%   x0[Nx1]: initial condition of the state
%   h[1x1]:  time-step size
%
% Outputs:
%   x[NxY]  :    System state propagated for Y time istants
%   fevals[1x1] : number of function evaluations
%   times[Yx1]:  time istants of state computation vector

% Initialize state and time variables
x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);

% Check for non-integer number of steps, adapting the final step

if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

% Compute each timestep size (necessary for non-integer number of steps)

hs = diff(times);

for i=1:(length(times)-1)    
    h = hs(i);
    fk = f(t,x);
    x = x + h*fk;
    t = t + h;
end

fevals = i;

end

function [x,fevals,times] = RK2(f,tspan,x0,h)
%RK2  integrate non-linear system dynamics using the Runge-Kutta 2 (Heun's
%method) integration
%
%
% Author
%   Davide Iafrate 
% 
% Inputs:
%   f : function handle of the RHS of the system dynamics, f(t,x)
%   tspan[1x2]  : vector of time interval extremes
%   x0[Nx1]: initial condition of the state
%   h[1x1]:  time-step size
%
% Outputs:
%   x[NxY]  :    System state propagated for Y time istants
%   fevals[1x1] : number of function evaluations
%   times[Yx1]:  time istants of state computation vector

% Initialize state and time variables
x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);

% Check for non-integer number of steps, adapting the final step

if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

% Compute each timestep size (necessary for non-integer number of steps)

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
    
    % Corrector step
    x(:,i+1) = xk + h*(0.5*fk + 0.5*fp);
    
end

fevals = 2*i;

end

function [x,fevals,times] = RK4(f,tspan,x0,h)
%RK4  integrate non-linear system dynamics using the Runge-Kutta 4
%integrator
%
%
% Author
%   Name: DAVIDE 
%   Surname: IAFRATE
% 
% Inputs:
%   f : function handle of the RHS of the system dynamics, f(t,x)
%   tspan[1x2]  : vector of time interval extremes
%   x0[Nx1]: initial condition of the state
%   h[1x1]:  time-step size
%
% Outputs:
%   x[NxY]  :    System state propagated for Y time istants
%   fevals[1x1] : number of function evaluations
%   times[Yx1]:  time istants of state computation vector

% Initialize state and time variables
x = x0;
times = tspan(1):h:tspan(end);
t = tspan(1);

% Check for non-integer number of steps, adapting the final step
if mod(diff(tspan),h)~=0
    times=[times tspan(end)];
end

% Compute each timestep size (necessary for non-integer number of steps)
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
    
    % Update time
    t = t + h;    
end

% 4 function evaluations per step are required
fevals = 4*k;

end

function A_A = A(alfa)
%A  return the system dynamics matrix as function of the angle alfa
%
%
% Author:
%    Davide Iafrate
% 
% Inputs:
%   alfa : [1,1] scalar angle   [rad]
%
% Outputs:
%   A_A[2x2]  :    System dynamics matrix

A_A =[0 1; -1 2*cos(alfa)];
end

function F = ffwd(A,h,algor,th)
%ffwd Function that returns the F operator of the selected algorithm s.t.
%       x_(k+1) = F*xk
% 
% 
% AUTHOR:
%   Davide Iafrate
% 
% Inputs:
%   A[NxN]          system dynamics matrix
%   h[1x1]          timestep size
%   algor[1x1]      scalar value to select the integration algorithm
%                    |-> 1: Runge-Kutta 1 (Forward Euler)
%                    |-> 2: Runge-Kutta 2 (Heun)
%                    |-> 4: Runge-Kutta 4
%   th : [1,1]  *optional* scalar value of the theta parameter for the 
%               theta-methods
% 
% Outputs:
%   F[NxN]          integration method F-matrix

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
        F = (I - A*(1-th)*h + 0.5*(A*(1-th)*h)^2)\...
            (I + A*th*h + 0.5*(A*th*h)^2);
        
end

end

function [Rhat,X,Y] = easystability(algorithm,th)
%EASYSTABILITY      function computing the value of the he absolute value 
%   of the mapped eigenvalues for the chosen integration method, 
%   in the range real(hl)=[-5 5], imag(hl)=[-5 5]
%   The output can be directly fed into the contour function as
%       contour(X,Y,Rhat) in order to produce the marginal stability
%       isolines of the method
%
%
% Author:
%    Davide Iafrate
% 
% Inputs:
%   algorithm : [1,1] scalar value to select the alcorithm   [-]
%                       |-> 2: Runge-Kutta 2 (Heun)
%                       |-> 4: Runge-Kutta 4
%                       |-> 5: BI2_theta
% 
%   th : [1,1]  *optional* scalar value of the theta parameter for the 
%               theta-methods
%
% Outputs:
%   Rhat[1001x1001]  :    method's operator absolute value matrix
%   X[1001x1001]     :    matrix where each row is the vector -5:0.01:5
%   Y[1001x1001]     :    matrix where each row is the vector -5:0.01:5

[X,Y] = meshgrid(-5:0.01:5,-5:0.01:5);
% Make the mesh the imaginary plane
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
% absolute value of the complex-valued R matrix
Rhat = abs(R);

end

function [x_end,fevals] = wrapper(f,tspan,x0,h,algor)
%WRAPPER  Wrapper function for the integrators, returning the state at the
%         final integration step
%
% Author:
%   Davide Iafrate 
% 
% Inputs:
%   f : function handle of the RHS of the system dynamics, f(t,x)
%   tspan[1x2]  : vector of time interval extremes
%   x0[Nx1]: initial condition of the state
%   h[1x1]:  time-step size
%
% Outputs:
%   x[Nx1]  :    System state propagated for Y time istants
%   fevals[1x1] : number of function evaluations

% Integrators selection
INT = {@RK1,@RK2,@RK4};

[xx,fevals,~] = INT{algor}(f,tspan,x0,h);
x_end=xx(:,end);
end
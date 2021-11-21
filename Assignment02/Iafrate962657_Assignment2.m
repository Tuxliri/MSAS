% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Ex 1
clearvars; close all; clc

k = 3.1416;
b = 2.7120;

% Set initial conditions for the state
y0 = [0; 0; 0; 0];
tspan = [0 10];

[tt,xx] = ode45(@(t,y) jck(t,y,b,k),tspan,y0);

dy = jck(tt,xx',b,k);
ddth1 = dy(3,:);
ddth2 = dy(4,:);

figure()
plot(tt,ddth1,tt,ddth2,'LineWidth',2)
grid minor
ylabel('Angular acceleration [rad/s^2]')
xlabel('time [s]')
legend('$\ddot{\theta}_1$','$\ddot{\theta}_2$','Interpreter','latex','FontSize',14)
F = open('data.mat');
data = F.data;

% parameters estimation
x0 = [2;3];
sol = lsqnonlin(@(y) costFcn(y,data),x0,[0;0],[inf;inf]);


b = sol(1);
k = sol(2);

tspan = linspace(0,10,1001);
[tt,xx] = ode45(@(t,y) jck(t,y,b,k),tspan,y0);
dy = jck(tt,xx',b,k);
ddth1 = dy(3,:);
ddth2 = dy(4,:);

figure()
semilogy(tt,abs(ddth1'-data(:,2)),tt,abs(ddth2'-data(:,3)),'LineWidth',1) 
ylabel('$|\ddot{\theta}_{1,2}^N - y_{1,2}^N|$','Interpreter','latex','FontSize',12)
xlabel('time [s]')
grid on
legend('$|\ddot{\theta}_1^N-y_1^N|$','$|\ddot{\theta}_2^N-y_2^N|$','Interpreter','latex','FontSize',14)

%% Ex 2
clearvars; close all; clc

% Values from text
% Isothermal expansion of gas
PmInf = 2.5*10^6;
VmInf = 10*(0.1)^3;
P0 = 21e6;

% for the isothermal transformation pV=const.
V0 = VmInf*PmInf/P0;

%%%%%%%%% Data input %%%%%%%%
% Fluid
data.skydrol.rho = 890;

% Accumulator
data.accumulator.V0 = V0;
data.accumulator.P0 = P0;
data.accumulator.gamma = 1.2;
data.accumulator.kA = 1.12;
data.accumulator.Vmax = 0.0024;

% Actuator
data.actuator.Atop = 1/4*pi*(50e-3)^2;
data.actuator.Abottom = 1/4*pi*((50e-3)^2-(22e-3)^2);
data.load.F0 = 1000;
data.load.k = 120e3;
data.actuator.cc_max = 0.2;
data.actuator.m = 2;

% Distributor
data.distributor.d0 = 5e-3;
data.distributor.kd = 12;
data.distributor.P4_0 = data.load.F0/data.actuator.Atop;
data.distributor.P5_0 = 0.1e6; %data.load.F0/data.actuator.Abottom;

% Delivery line
data.checkValve.kCV = 2;
data.deliveryLine.L = 2;
data.deliveryLine.D = 18e-3;
data.deliveryLine.f = 0.032;

% Return line
data.returnLine.L = 15;
data.returnLine.D = 18e-3;
data.returnLine.f = 0.035;

% Tank
data.tank.P_T = 0.1e6;
data.tank.kT = 1.12;

% Command
data.command.t1 = 1;
data.command.t2 = 1.5;

% Integration
Vt0 = 1e-3;
x0 = 0;
vx0 = 0;
Vacc0 = 0;

% System state initial conditions

x0 = [Vt0;
      x0;
      vx0;
      Vacc0];

% Final integration time
tf = 3;

% Integration timespan
tspan = [0 tf];

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);


opts = odeset(opts,'event',@(t,y) actuator_event(t,y,data));
[tt1,xx1,te,xxe,ie] = ode15s(@(t,x) exercise2ODE(t,x,data),tspan,x0,opts);

% Integrate from te to tf setting the velocity to 0
x0 = xxe;
x0(3) = 0;
tspan = [te tf];
[tt2,xx2] = ode15s(@(t,x) exercise2ODE(t,x,data),tspan,x0,opts);
tt = [tt1;tt2];
xx = [xx1; xx2];

% Retrieve system response at each time instant
for i =1:length(tt)
    [dy,parout(i,:)] = exercise2ODE(tt(i),xx(i,:),data);
end

disp('--> parameters evaluation done...')

% Extract circuit parameters

Q1 = parout(:,1);
Q2 = parout(:,2);
Q3 = parout(:,3);
Q4 = parout(:,4);
Q5 = parout(:,5);
Q6 = parout(:,6);
Q7 = parout(:,7);
Qacc = parout(:,8);
p1 = parout(:,9);
p2 = parout(:,10);
p3 = parout(:,11);
p4 = parout(:,12);
p5 = parout(:,13);
p6 = parout(:,14);
p7 = parout(:,15);
pAcc = parout(:,16);
Av = parout(:,17);
u = parout(:,18);

figure()
plot(tt,xx(:,2),'LineWidth',2)
ylabel('Piston position $[m]$','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')
grid minor

figure()
plot(tt,xx(:,3),'LineWidth',2)
ylabel('Piston velocity $[m/s]$','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')
grid minor

figure()
plot(tt,-xx(:,4),'LineWidth',2)
hold on
plot(tt,xx(:,1),'LineWidth',2)
xlabel('time [s]','Interpreter','latex')
ylabel('Volume $[m^3]$','Interpreter','latex')
grid minor
legend('$V_{accumulator}$','$V_{tank}$','Interpreter','latex')

figure()
hold on, grid minor, box on
plot(tt,Q3.*tt*1e3,tt,Q5.*tt*1e3)
xlabel('time [s]')
ylabel('voumetric flow rate $\mathrm{\left[\frac{l}{s}\right]}$', 'interpreter', 'latex')
legend('$\mathrm{Q_3}$','$\mathrm{Q_6}$', 'interpreter', 'latex')

figure()
hold on, grid minor, box on
[ylab,~,~] = plotyy(tt,u,tt,Av*1e6);
xlabel('time [s]')
ylabel(ylab(1),'valve position [-]')
ylabel(ylab(2),'valve orifice opening [mm^2]')

figure()
hold on, grid minor, box on
plot(tt,pAcc*1e-6,tt,p1*1e-6,...
    tt,p2*1e-6,tt,p3*1e-6,tt,p4*1e-6,tt,p5*1e-6)
xlabel('time [s]')
ylabel('pressure [MPa]')
legend('$\mathrm{p_{Acc}}$','$\mathrm{p_1}$','$\mathrm{p_2}$',...
    '$\mathrm{p_3}$','$\mathrm{p_4}$','$\mathrm{p_5}$','interpreter', 'latex')

figure()
hold on, grid minor, box on
plot(tt,p6*1e-6,tt,p7*1e-6,tt,data.tank.P_T*ones(size(tt))*1e-6,'-.')
xlabel('time [s]')
ylabel('pressure [MPa]')
legend('$\mathrm{p_6}$','$\mathrm{p_7}$','$\mathrm{p_{tank}}$', 'interpreter', 'latex')

% figure()
% hold on, grid minor, box on
% plot(tt,(pAcc-p1)*1e-6,...
%     tt,(pAcc-p2)*1e-6,tt,(pAcc-p3)*1e-6,tt,(pAcc-p4)*1e-6)%,tt,(pAcc-p5)*1e-6)
% xlabel('time [s]')
% ylabel('pressure [MPa]')
% legend('$\mathrm{p_{Acc}}$','$\mathrm{p_1}$','$\mathrm{p_2}$',...
%     '$\mathrm{p_3}$','$\mathrm{p_4}$','$\mathrm{p_5}$','interpreter', 'latex')

%% Ex 3
clearvars; close all; clc

% Circuit parameters
R1 = 1000;
R2 = 100;
L = 1e-3;
C = 1e-3;

% System dynamics matrix
A = [0 1;
    -(1/L/C)/(1+R2/R1) -(1/R1/C+R2/L)/(1+R2/R1)];

eig(A)

% Initial conditions
tStart = tic;
y0 = [1; 0];
[tt,yy] = ode23s(@RLC,[0 1],y0);
plot(tt,yy(:,1),'LineWidth',1.5)
ylabel('$V_c(t)$ [V]','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
t_free = toc(tStart);

grid minor

y0 = [1; 0];

% Compare the computational time of the non-stiff and stiff integrators
tStart = tic;     
[tt,yy] = ode45(@forcedRLC,[0 10],y0);
t_ode45 = toc(tStart);

tStart = tic;     
[tt,yy] = ode23s(@forcedRLC,[0 10],y0);
t_ode23 = toc(tStart);

figure()
plot(tt,yy(:,1),'LineWidth',1.5)
grid minor
ylabel('$V_c(t)$ [V]','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')

%% Ex 4
clearvars; close all; clc

y0 = 20*(ones(6,1));
t=0;

tspan = [0:0.2:0.9 1:1:5 10];
for i=1:length(tspan)
    legendentries{i}=sprintf('$t=%0.5g$',tspan(i));
end

[tt,xx] = ode23s(@nozzle,tspan,y0);

% Plot temperature profiles
y = [0.4 0.4];    


figure(1)
hold on

X = [0 0.5 3 5.5 5.5 8.1 10.5 11];
yy = [20+980*min(tt,1) xx 20*ones(length(tt),1)];
plot(X,yy(:,:),'LineWidth',1.5)
xline(X,'--r')
annotation('textarrow',[0.2 0.3],[0.2 0.9],'String','$t$',...
    'Interpreter','latex','LineWidth',1,'FontSize',14)

legend(legendentries{:},'Interpreter','latex')
ylabel('$temperature\; T_i(t)\; [^{\circ}C]$','Interpreter','latex') 
xlabel('$x\; [mm]$','Interpreter','latex') 
xlim([X(1) X(end)])

% 2 nodes per element
y0 = 20*(ones(8,1));
t=0;
[tt1,xx1] = ode23s(@nozzle2,tspan,y0);


figure(2)

X = [0 0.5:(5/3):5.5 5.5:(5/3):10.5 11];
yy1 = [20+979*min(tt1,1) xx1 20*ones(length(tt1),1)];

plot(X,yy1(:,:),'LineWidth',1.5)
annotation('textarrow',[0.2 0.3],[0.2 0.9],'String','$t$',...
    'Interpreter','latex','LineWidth',1,'FontSize',14)
xline(X,'--r')

legend(legendentries{:},'Interpreter','latex')
ylabel('$temperature\; T_i(t)\; [^{\circ}C]$','Interpreter','latex') 
xlabel('$x\; [mm]$','Interpreter','latex') 
xlim([X(1) X(end)])

%% Functions
function dy = jck(~,y,b,k)
%JCK function computing the dynamics of the mechanical system in exercise 1
J1 = 0.2;
J2 = 0.1;
T0 = 0.1;

th1 = y(1,:);
th2 = y(2,:);
dth1 = y(3,:);
dth2 = y(4,:);

ddth1(1,:) = k.*(th2-th1)/J1;
ddth2(1,:) = (k.*(th1-th2)-b.*sign(dth2).*dth2.^2+T0)/J2;

dy = [dth1;
      dth2;
      ddth1;
      ddth2];
  
end

function cost = costFcn(y,data)
    b = y(1);
    k = y(2);
    
    % Set initial conditions for the state
    y0 = [0; 0; 0; 0];
    tspan = linspace(0,10,1001);

    [tt,xx] = ode45(@(t,y) jck(t,y,b,k),tspan,y0);

    J1 = 0.2;
    J2 = 0.1;
    T0 = 0.1;
    ddth1 = k*(xx(:,2)-xx(:,1))/J1;
    ddth2 = (k*(xx(:,1)-xx(:,2))-b*sign(xx(:,4)).*xx(:,4).^2+T0)/J2;

    cost = [ddth1-data(:,2);
            ddth2-data(:,3)];
end

function dy = RLC(t,y)
R1 = 1000;
R2 = 100;
L = 1e-3;
C = 1e-3;

Vc = y(1);
dVc = y(2);

dy = [dVc;
     1/(1+R2/R1)*(-(1/R1 + R2*C/L)*dVc - 1/L*Vc)];
end

function dy = forcedRLC(t,y)
R1 = 1000;
R2 = 100;
L = 1e-3;
C = 1e-3;
f = 5;

Vc = y(1);
dVc = y(2);

vt = sin(2*pi*f*t)*atan(t);
dvt = 2*pi*f*cos(2*pi*f*t)*atan(t) + sin(2*pi*f*t)/(1 + t^2);

dy = [y(2);
     1/(C*(1+R2/R1))*(-1/L*Vc - (1/R1 + R2*C/L)*dVc + dvt/R1 + vt/L)];
end

function dy = nozzle(t,y)
% segments length
Ti = 20 + 980*min(t,1);
To = 20;

Ac = 1;
L =1;

% thermal conductivities
k = [401 45 45 100 0.56 0.56 401];
C = [0 468 468 468 100 100 100 0];
rho = [0 8900 8900 8900 600 600 600 0];

% Extract temperatures from state vector
T2 = y(1);
T3 = y(2);
T4 = y(3);
T5 = y(4);
T6 = y(5);
T7 = y(6);% Generate temperatures vector
T = [Ti; T2; T3; T4; T5; T6; T7; To];
DX = [0.5 2.5 2.5 0.1 2.5 2.5 0.5]*1e-3;
DXmass = [0 5/4 5/2 5/4 5/4 5/2 5/4 0]*1e-3;
m=rho.*DXmass;
R = DX./k;
R(4) = 1/5000;

dT = zeros(length(y)+2,1);
for i=2:7
    dT(i) = 1/(m(i)*C(i))*((T(i-1)-T(i))/R(i-1) + (T(i+1)-T(i))/R(i));
end

dy=dT(2:7);
end

function dy = nozzle2(t,y)
% segments length
Ti = 20 + 980*min(t,1);
To = 20;

Ac = 1;
L =1;

% thermal conductivities
k = [401 45 45 45 100 0.56 0.56 0.56 401];
C = [0 468 468 468 468 100 100 100 100 0];
rho = [0 8900 8900 8900 8900 600 600 600 600 0];

% Extract temperatures from state vector
T2 = y(1);
T3 = y(2);
T4 = y(3);
T5 = y(4);
T6 = y(5);
T7 = y(6);
T8 = y(7);
T9 = y(8);

% Generate temperatures vector
T = [Ti; T2; T3; T4; T5; T6; T7; T8; T9; To];
DX = [0.5 5/3 5/3 5/3 0.1 5/3 5/3 5/3 0.5]*1e-3;

DXmass = [0 5/6 5/3 5/3 5/6 5/6 5/3 5/3 5/6 0]*1e-3;

m = rho.*DXmass;
R = DX./k;
R(5) = 1/5000;

dT = zeros(length(y)+2,1);
for i=2:9
    dT(i) = 1/(m(i)*C(i))*((T(i-1)-T(i))/R(i-1) + (T(i+1)-T(i))/R(i));
end

dy=dT(2:9);
end

function [dy,outpar] = exercise2ODE(t,y,data)

% Extract variables from state
x = y(2);
vx = y(3);
Vacc = y(4);

% pilot piston boundary conditions
x_max = data.actuator.cc_max;

if x < 0
    x = 0;
end
if x <= 0 && vx < 0
    x = 0;
    vx = 0;

end
if x > x_max
    x = x_max;
end
if x >= x_max && vx > 0
    x = x_max;
    vx = 0;
end

m = data.actuator.m;

% flows due to piston
Atop = data.actuator.Atop;
Abottom = data.actuator.Abottom;

Q4 = Atop*vx;
Q5 = -Abottom*vx;

% DISTRIBUTOR model
d0 = data.distributor.d0;

% control history
u = command(t,data);
alpha = 2*acos(1 - 2*abs(u));

Av = alpha*d0^2/8 - d0*(0.5 - u)*d0/2*sin(alpha/2);

% Flows
if u == 0
    Q3 = 0;
    Q6 = 0;
end
if u > 0
    Q3 = Q4;
    Q6 = Q5;
end
if u < 0
    Q3 = Q5;
    Q6 = Q4;
end

% Fluid data
rho = data.skydrol.rho;

% ACCUMULATOR
V0 = data.accumulator.V0;
P0 = data.accumulator.P0;
gamma = data.accumulator.gamma;
Pacc = P0*(V0/(V0+Vacc))^gamma;

% Delivery line
kA = data.accumulator.kA;
deliveryArea = 1/4*pi*(data.deliveryLine.D)^2;

P1 = Pacc - 0.5*kA*rho*Q3/deliveryArea*abs(Q3/deliveryArea);

% CHECK valve
kCV = data.checkValve.kCV;

P2 = P1 - 0.5*kCV*rho*Q3/deliveryArea*abs(Q3/deliveryArea);

% Delivery Line losses
v23 = Q3/deliveryArea;
L23 = data.deliveryLine.L;
D23 = data.deliveryLine.D;
f23 = data.deliveryLine.f;

P3 = P2 - 0.5*f23*L23/D23*rho*v23*abs(v23);

% TANK
returnArea = 1/4*pi*(data.returnLine.D)^2;
P_T = data.tank.P_T;
kT = data.tank.kT;
v67 = Q6/returnArea;

P7 = P_T - 0.5*kT*rho*Q6/returnArea*abs(Q6/returnArea);

% Return Line losses
L67 = data.returnLine.L;
D67 = data.returnLine.D;
f67 = data.returnLine.f;

P6 = P7 - 0.5*f67*L67/D67*rho*v67*abs(v67);

%

% Pressure drops
kd = data.distributor.kd;
P4_0 = data.distributor.P4_0;
P5_0 = data.distributor.P5_0;

if u == 0
    P4 = P4_0;
    P5 = P5_0;
elseif u > 0
    if Av < 1e-9
        P4 = P4_0;
        P5 = P5_0;
    else             
        P4 = P3 - kd*.5*rho*abs(Q3)*Q3/Av^2;
        P5 = P6 - kd*.5*rho*abs(Q6)*Q6/Av^2;
    end
% else
%     if a_0 < 1e-9
%         p_4 = p_eq0;
%         p_6 = data.system.p_i;
%     else
%         p_6 = p_5 - data.closed_centre_valve.k*.5*data.skydrol.rho*abs(Q_5)*Q_5/a_0^2;
%         p_4 = p_8 - data.closed_centre_valve.k*.5*data.skydrol.rho*abs(Q_8)*Q_8/a_0^2;
%     end
end

% Load
F0 = data.load.F0;
k = data.load.k;

Fx = F0 + k*x;

% States derivatives
dVT = -Q6;
dx = vx;
ddx = 1/m*(P4*Atop - P5*Abottom - Fx);
dVacc = Q3;

% accumulator phyhsical constraints
if Vacc >= data.accumulator.Vmax
    dVacc = 0;
end

% pilot piston physical constraints
if (x >= x_max && ddx > 0) || (x <= 0 && ddx < 0)
    ddx = 0;
end

dy = [dVT;
      dx;
      ddx;
      dVacc];
  
% Parameters output
Q2 = Q3;
Q1 = Q2;
Qacc = Q1;
Q7 = Q6;

outpar = [Q1 Q2 Q3 Q4 Q5 Q6 Q7 Qacc P1 P2 P3 P4 P5 P6 P7 Pacc Av u];

end

function u = command(t,data)
% Distributor control command

t1 = data.command.t1;
t2 = data.command.t2;

if t <= t1
    u = 0;
elseif t > t1 && t < t2
    u = (t - t1)/(t2 - t1);
elseif t >= t2
    u = 1;
end

end

function [value,isterminal,direction] = actuator_event(t,y,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
xmax = data.actuator.cc_max;
x = y(2);
value = xmax-x;
isterminal = 1;
direction = -1;
end

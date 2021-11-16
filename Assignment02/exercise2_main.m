clc 
clear all 
close all

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
data.load.k = 120;
data.actuator.cc_max = 0.2;
data.actuator.m = 2;

% Distributor
data.distributor.d0 = 5e-3;
data.distributor.kd = 12;
data.distributor.P4_0 = 0;
data.distributor.P5_0 = 0; %data.load.F0/data.actuator.Abottom;

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

%% Integration
Vt0 = 1e-3;
x0 = 0;
vx0 = 0;
Vacc0 = 0;

x0 = [Vt0;
      x0;
      vx0;
      Vacc0];
  
tspan = [0 3];

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[tt,xx] = ode23s(@(t,x) exercise2ODE(t,x,data),tspan,x0,opts);


% opts = odeset(opts,'event',@(t,y) actuator_event(t,y,data));
% [tt,xx,te,xxe,ie] = ode23s(@(t,x) exercise2ODE(t,x,data),tspan,x0,opts);

figure()
plot(tt,xx(:,2))
title('x piston')

figure()
plot(tt,xx(:,3))
title('v piston')


figure()
plot(tt,xx(:,4))
title('V acc')

figure()
plot(tt,xx(:,1))
title('V tank')
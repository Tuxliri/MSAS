function [dy,outpar] = exercise2ODE(t,y,data)

%% Extract variables from state
V_T = y(1);
x = y(2);
vx = y(3);
Vacc = y(4);

%% flows due to piston
Atop = data.actuator.Atop;
Abottom = data.actuator.Abottom;

Q4 = Atop*vx;
Q5 = Abottom*vx;

%% pilot piston boundary conditions
cc_max = data.actuator.cc_max;

if x < 0
    x = 0;
end
if x <= 0 && vx < 0
    x = 0;
    vx = 0;
    Q4 = 0;
    Q5 = 0;
end
if x > cc_max
    x = cc_max;
end
if x >= cc_max && vx > 0
    x = cc_max;
    vx = 0;
    Q4 = 0;
    Q5 = 0;
end

m = data.actuator.m;

%% DISTRIBUTOR model
d0 = data.distributor.d0;

% control history
u = command(t,data);
alpha = 2*acos(1 - abs(u));
Av = d0^2/4*(alpha - sin(alpha));

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

%% Fluid data
rho = data.skydrol.rho;

%% ACCUMULATOR
V0 = data.accumulator.V0;
P0 = data.accumulator.P0;
gamma = data.accumulator.gamma;
Pacc = P0*(V0/(V0-Vacc))^gamma;

%% Delivery line
kA = data.accumulator.kA;
deliveryArea = 1/4*pi*(data.deliveryLine.D)^2;

P1 = Pacc - 0.5*kA*rho*Q3/deliveryArea*abs(Q3/deliveryArea);

%% CHECK valve
kCV = data.checkValve.kCV;

P2 = P1 - 0.5*kCV*rho*Q3/deliveryArea*abs(Q3/deliveryArea);

%% Delivery Line losses
v23 = Q3/deliveryArea;
L23 = data.deliveryLine.L;
D23 = data.deliveryLine.D;
f23 = data.deliveryLine.f;

P3 = P2 - 0.5*f23*L23/D23*rho*v23*abs(v23);

%% TANK
returnArea = 1/4*pi*(data.returnLine.D)^2;
P_T = data.tank.P_T;
kT = data.tank.kT;
v67 = Q6/returnArea;

P7 = P_T + 0.5*kT*rho*Q6/returnArea*abs(Q6/returnArea);

%% Return Line losses
L67 = data.returnLine.L;
D67 = data.returnLine.D;
f67 = data.returnLine.f;

P6 = P7 + 0.5*f67*L67/D67*rho*v67*abs(v67);

%%

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
        P5 = P6 + kd*.5*rho*abs(Q6)*Q6/Av^2;
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

%% Load
F0 = data.load.F0;
k = data.load.k;

Fx = F0 + k*x;
%% States derivatives
dVT = Q6;
dx = vx;
ddx = 1/m*(P4*Atop - P5*Abottom - Fx);
dVacc = Q3;

% accumulator phyhsical constraints
if Vacc >= data.accumulator.Vmax
    dVacc = 0;
end

% pilot piston physical constraints
if (x >= cc_max && ddx > 0) || (x <= 0 && ddx < 0)
    ddx = 0;
end

dy = [dVT;
      dx;
      ddx;
      dVacc];
  
end


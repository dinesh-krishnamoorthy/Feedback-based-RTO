function [sys,par,F] = WilliamOtto6state(Ts)
%
% Simplified Model with 2 reactions, ignoring the intermediate compnent C.
% Structural mismatch with WilliamOtto(Ts) function.
%
% Model taken from Zhang, Y. and Forbes, J.F., 2000. 
% Extended design cost: a performance criterion for 
% real-time optimization systems. Comput & Chem Eng, 24(8), pp.1829-1841.
%
% Written by: Dinesh Krishnamoorthy, Sept 2019 NTNU

import casadi.*


xa = MX.sym('xa',1);
xb = MX.sym('xb',1);
xc = MX.sym('xc',1);
xp = MX.sym('xp',1);
xe = MX.sym('xe',1);
xg = MX.sym('xg',1);

Fa = MX.sym('Fa',1);
Fb = MX.sym('Fb',1);
Tr = MX.sym('Tr',1);

Vr = 2105;% 2.63; % Reactor mass holdup - constant

CC.tau1 = 150;
CC.tauC = 150;
CC.Kp = CC.tau1/(0.00517*CC.tauC);

k01 = 1.6599e6;
k02 = 7.2117e8;
k03 = 2.6745e12;

B1 = 6666.7; % degK
B2 = 8333.3;
B3 = 11111;

k1 = k01*exp(-B1/(Tr));
k2 = k02*exp(-B2/(Tr));
k3 = k03*exp(-B3/(Tr));

dxa = (Fa - (Fa+Fb)*xa - Vr*xa*xb*k1)/Vr;
dxb = (Fb - (Fa+Fb)*xb - Vr*xa*xb*k1 - Vr*xb*xc*k2)/Vr;
dxc = -(Fa+Fb)*xc/Vr + 2*xa*xb*k1 - 2*xb*xc*k2 - xc*xp*k3;
dxp = -(Fa+Fb)*xp/Vr + xb*xc*k2 - 0.5*xp*xc*k3;
dxe = -(Fa+Fb)*xe/Vr + 2*xb*xc*k2;
dxg = -(Fa+Fb)*xg/Vr + 1.5*xp*xc*k3;

sys.diff = vertcat(dxa,dxb,dxc,dxp,dxe,dxg);
sys.x = vertcat(xa,xb,xc,xp,xe,xg,Tr);
sys.d = vertcat(Fa);
sys.u = vertcat(Fb,Tr);

sys.L = -(1043.38*xp*(Fa+Fb)+20.92*xe*(Fa+Fb) - 79.23*Fa - 118.34*Fb);

sys.nlcon = vertcat(xa,xg);
sys.lb = vertcat(-Inf,-Inf);
sys.ub = vertcat(0.12,0.08);

sys.y = sys.x;

ode = struct('x',sys.x,'p',vertcat(sys.d,sys.u),'ode',sys.diff,'quad',sys.L);
opts = struct('tf',Ts);

% create IDAS integrator
F = integrator('F','cvodes',ode,opts);

par.lbx = zeros(6,1);
par.ubx = 1.*ones(6,1);
par.lbu = [2;20+273];
par.ubu = [7;473];
par.dx0 = [1;1;1;1.0;1.0;1.0];
par.u0 = [0;373];
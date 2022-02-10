function [sys,par] = WilliamOtto2reaction
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
xp = MX.sym('xp',1);
xe = MX.sym('xe',1);
xg = MX.sym('xg',1);

Fa = MX.sym('Fa',1);
Fb = MX.sym('Fb',1);
Tr = MX.sym('Tr',1);

Vr = 2105;% 2.63; % Reactor mass holdup - constant

k01 = 1.655e8;
k02 = 2.611e13;

B1 = 8077.6; % degK
B2 = 12438.5;

k1 = k01*exp(-B1/Tr);
k2 = k02*exp(-B2/Tr);

CC.tau1 = 150;
CC.tauC = 150;
CC.Kp = CC.tau1/(0.00517*CC.tauC);
    
dxa = (Fa - (Fa+Fb)*xa - Vr*xa*xb^2*k1 - Vr*xa*xb*xp*k2)/Vr;
dxb = (Fb - (Fa+Fb)*xb - 2*Vr*xa*xb^2*k1 - Vr*xb*xb*xp*k2)/Vr;
dxp = -(Fa+Fb)*xp/Vr + xa*xb^2*k1 - xa*xb*xp*k2;
dxe = -(Fa+Fb)*xe/Vr + 2*xa*xb^2*k2;
dxg = -(Fa+Fb)*xg/Vr + 3*xa*xb*xp*k2;
dTr = -CC.Kp*dxg + CC.Kp/CC.tau1*(0.08-xg);

sys.diff = vertcat(dxa,dxb,dxp,dxe,dxg,dTr);
sys.x = vertcat(xa,xb,xp,xe,xg,Tr);
sys.d = vertcat(Fa);
sys.u = vertcat(Fb);
sys.y = sys.x;

sys.L = -(1043.38*xp*(Fa+Fb)+20.92*xe*(Fa+Fb) - 79.23*Fa - 118.34*Fb);

sys.nlcon = vertcat(xa,xg);
sys.lb = vertcat(-Inf,-Inf);
sys.ub = vertcat(0.12,0.08);

par.lbx = [zeros(5,1);20+273];
par.ubx = [1.*ones(5,1);473];
par.lbu = [2];
par.ubu = [7];
par.dx0 = [1;1;1.0;1.0;1.0;373];
par.u0 = [0];
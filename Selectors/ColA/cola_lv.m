function [sys,F_integrator,f,par] = cola_lv(Ts)
%
% colamod - This is a nonlinear model of a distillation column with
%           NT-1 theoretical stages including a reboiler (stage 1) plus a
%           total condenser ("stage" NT). The liquid flow dynamics are
%           modelled by a simple linear relationship.
%           Model assumptions: Two components (binary separation); constant
%           relative volatility; no vapor holdup; one feed and two products;
%           constant molar flows (same vapor flow on all stages); 
%           total condenser
%
%           The model is based on column A in Skogestad and Postlethwaite
%           (1996). The model has 82 states.
% 
%           Re-written using CasADi v3.5.1 by D. Krishnamoorthy (2019)
%
% 

import casadi.*
global pF pB pD 

%------------------------------------------------------------
% The following data need to be changed for a new column.
% These data are for "colmn A".
% Number of stages (including reboiler and total condenser: 
    NT=41; 
% Location of feed stage (stages are counted from the bottom):
    NF=21;
% Relative volatility
    alpha=1.5;
% Nominal liquid holdups
    M0(1)=0.5;      	% Nominal reboiler holdup (mol)
    i=2:NT-1; M0(i)=0.5.*ones(1,NT-2);% Nominal stage (tray) holdups (mol)
    M0(NT)=0.5;      	% Nominal condenser holdup (mol)
% Data for linearized liquid flow dynamics (does not apply to reboiler and condenser):
    taul= 0.063;     	% time constant for liquid dynamics (s)
    F0=1;               % Nominal feed rate (mol/s) 
    qF0 = 1;            % Nominal fraction of liquid in feed 
    L0=2.70629;     	% Nominal reflux flow (from steady-state data)
    L0b=L0 + qF0*F0;	% Nominal liquid flow below feed (mol/s)
    lambda=0;           % Effect of vapor flow on liquid flow ("K2-effect")
    V0=3.20629;V0t=V0+(1-qF0)*F0;% Nominal vapor flows - only needed if lambda is nonzero 
% End data which need to be changed

% Prices
  
    pV  = MX.sym('pV');           
    
%------------------------------------------------------------

% Differential states
x = MX.sym('x',NT);         % Liquid composition from btm to top
M = MX.sym('M',NT);         % Liquid hold up from btm to top

% Inputs and disturbances

LT = MX.sym('LT');          % Reflux
VB = MX.sym('VB');          % Boilup

F  = MX.sym('F');           % Feedrate

zF = 0.5;% MX.sym('zF');          % Feed composition
% alpha = MX.sym('alpha');   % relative volatility
qF = 1; % MX.sym('qF');          % Feed liquid fraction

% Vapor-liquid equilibria
i=1:NT-1;    y=alpha*x(i)./(1+(alpha-1)*x(i));

% Vapor Flows assuming constant molar flows
i=1:NT-1;    V=VB*ones(1,NT-1);
i=NF:NT-1;   V(i)=V(i) + (1-qF)*F;

% Liquid flows assuming linearized tray hydraulics with time constant taul
% Also includes coefficient lambda for effect of vapor flow ("K2-effect").
L = 0;
for i=2:NF
    L = [L; L0b + (M(i)-M0(i)')./taul + lambda.*(V(i-1)-V0)'];
end

for i=NF+1:NT-1
    L = [L; L0  + (M(i)-M0(i)')./taul + lambda.*(V(i-1)-V0t)'];
end
L= [L;LT];

   % P-Controllers for control of reboiler and condenser hold up.
    par.KcB=10;  par.KcD=10;         % controller gains
    par.MDs=0.5; par.MBs=0.5;        % Nominal holdups - these are rather small
    par.Ds=0.5;  par.Bs=0.5;          % Nominal flows
    MB=M(1);  MD=M(NT);      % Actual reboiler and condenser holdup
    D=par.Ds+(MD-par.MDs)*par.KcD;       % Distillate flow
    B=par.Bs+(MB-par.MBs)*par.KcB;       % Bottoms flow
    
    
% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Reboiler (assumed to be an equilibrium stage)
dMdt = L(2)      - V(1)      - B;
dMxdt= L(2)*x(2) - V(1)*y(1) - B*x(1);

% Column
for i=2:NT-1
dMdt = [dMdt; L(i+1)         - L(i)       + V(i-1)'         - V(i)'];
dMxdt= [dMxdt; L(i+1).*x(i+1) - L(i).*x(i) + V(i-1)'.*y(i-1) - V(i)'.*y(i)];
end

% Total condenser (no equilibrium stage)
dMdt = [dMdt; V(NT-1)         - LT       - D];
dMxdt = [dMxdt; V(NT-1)*y(NT-1) - LT*x(NT) - D*x(NT)];

% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(NF) =  dMdt(NF)  + F;
dMxdt(NF)= dMxdt(NF) + F*zF;

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
dxdt = [];
for i=1:NT
dxdt = [dxdt;(dMxdt(i) - x(i).*dMdt(i) )./M(i)];
end

T = 100-20*x; % Tray temperatures

% Build ODE
diff = vertcat(dxdt,dMdt);
x_var = vertcat(x,M);
d_var = vertcat(F,pV);
p_var = vertcat(LT,VB);

J = pF*F + pV*VB - pB*B - pD*D; 

nlcon = [];% vertcat(x_var(41),1-x_var(1));
lb = [];%[0.95;0.292];
ub = [];%[1;1];

f = Function('f',{x_var,p_var,d_var},{diff,J},{'x','p','d'},{'xdot','qj'});

ode = struct('x',x_var,'p',vertcat(p_var,d_var),'ode',diff,'quad',J);
opts = struct('tf',Ts);

% create CVODES integrator
F_integrator = integrator('F','cvodes',ode,opts);


par.NT = NT;

sys.x = x_var;
sys.u = p_var;
sys.d = d_var;
sys.dx = diff;
sys.y = vertcat(T,F,LT,VB,sys.x(41),1-sys.x(1),pV);

sys.L = J;

sys.nlcon = nlcon;
sys.lb = lb;
sys.ub = ub;

par.lbx = 1e-5*ones(2*par.NT,1);
par.ubx = 2*ones(2*par.NT,1); par.lbx(par.NT) = 0.95;
par.dx0 = 0.5*ones(2*par.NT,1);
par.lbu = [0;0;0];
par.ubu = [10;4.008;6];
par.u0  = [2.706;3.206;1];

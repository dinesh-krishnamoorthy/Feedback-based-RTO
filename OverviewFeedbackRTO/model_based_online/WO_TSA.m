clc
clear
import casadi.*

Ts = 1; % Sample time

% Load William-Otto reactor model
[sys,par,F] = WilliamOtto(Ts);
d_val = 1.3; % Disturbance Fa

% Compute the cost Jacobian function
Ju = Jacobian(sys);

% Initialize WO reactor
u_in = 3;
[xf,exitflag] = solveODE(sys,par,d_val,u_in);


%% Configure EKF

EKF = prepareEKF(sys,Ts);

ny = numel(sys.y);
nxEKF = numel(sys.x);
xk_hat = xf;
uEKF = u_in;
Pk = 1e3.*eye(nxEKF);
Qk = 1e3.*eye(nxEKF);
Rk = 1e0.*eye(ny);
K = 0;

%% Simulation

nIter = 4*3600;
h = waitbar(0,'Simulation in Progress...');
GC.u = u_in(1);
GC.err = 0;
GC.err0 = 0;
Ju_hat = 0;
for sim_k = 1:nIter
    waitbar(sim_k /nIter)
    
    % Gradient Controller
     GC.tauC = 200;
    GC.Kp = 250/(23.571*(GC.tauC+300)); %3.996 % theta = 300?
    GC.Ki = GC.Kp/(min(250,(4*GC.tauC+300)));
    
    GC.err = (0 - Ju_hat);
    GC.u = GC.u + GC.Kp*GC.err + GC.Ki*GC.err - GC.Kp*GC.err0;
    GC.err0 = GC.err;
    
     if sim_k>2*3600
        d_val = 1;
    end
    % Plant simulator
    u_in = vertcat(GC.u);
    Fk = F('x0',xf,'p',vertcat(d_val,u_in));
    xf = full(Fk.xf);
    
    % Gradient Estimator
    [xk_hat,Pk] = EKF_estimation(EKF,xf,xk_hat,u_in,Pk,Qk,Rk,d_val);
    Ju_hat = full(Ju(xf,u_in,d_val));
    
    sim.x(:,sim_k) = xf;
    sim.u(:,sim_k) = u_in;
    sim.Ju(:,sim_k) = Ju_hat;
    sim.c(sim_k) = Ju_hat;
    
end
close(h)
%%

figure(12)
subplot(221)
hold all
plot((sim.x(6,:)))
ylabel('x_G')
grid on
subplot(222)
hold all
% plot(sim.Ju)
plot(sim.c)
ylabel('J_u')
grid on
subplot(224)
hold all
plot((sim.u(1,:)))
ylabel('F_B')
grid on
subplot(223)
hold all
plot((sim.x(7,:)))
ylabel('T_r')
grid on
save('sim_TSA2','sim')

function [Ju,Gu] = Jacobian(sys)
% Function to evaluate the Jacobian of the cost and constraint
% Adetola & Guay (2007)
% William Otto example
%
% Written by: Dinesh Krishnamoorthy, Sep 2019

import casadi.*

J = sys.L; % Economic objective
G = {vertcat(sys.diff)};

rf_function = Function('rf_function',{sys.x,sys.u,sys.d},...
                    {vertcat(G{:}),J,vertcat(sys.x(1),sys.x(6))},{'x','u','d'},{'g','J','c'});
rf  = rootfinder('rf','newton',rf_function);

Ju = rf.factory('Ju',{'x','u','d'},{'jac:J:u'});
Gu = rf.factory('gu',{'x','u','d'},{'jac:c:u'});

end



clc
clear
import casadi.*

Ts = 1; % Sample time

% Load William-Otto reactor model
[sys,par,F] = WilliamOtto(Ts);
d_val = 1.3; % Disturbance Fa

% Compute nominal optimum
[x_opt,u_opt,exitflag] = SSOpt(sys,par,d_val);

% Compute the Jacobians and Hessians for variational calculus
[Juu,Jud,Gu,Gd] = WO_NE(sys,x_opt,u_opt,d_val);

% Initialize WO reactor
u_in = 3;
[xf,exitflag] = solveODE(sys,par,d_val,u_in);

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
    GC.Kp = 520/(32*(GC.tauC+300)); %3.996 % theta = 300?
    GC.Ki = GC.Kp/(min(520,(4*GC.tauC+300)));
    
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
    Ju_hat = Jud*pinv(Gd)*(xf-x_opt) + (Juu - Jud*pinv(Gd)*Gu)*(u_in-u_opt);
    
    sim.x(:,sim_k) = xf;
    sim.u(:,sim_k) = u_in;
    sim.Ju(:,sim_k) = Ju_hat;
    sim.c(sim_k) = Ju_hat;
    
end
close(h)
%%

figure(12)
hold all
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

save('sim_NEC','sim')


function [Juu,Jud,Gu,Gd] = WO_NE(sys,x_opt,u_opt,d_val)
% Function to evaluate the gradient using neighbhouring extremal control 
% Gros et al. 2009
% Written by: Dinesh Krishnamoorthy, Sep 2019

import casadi.*

G = {vertcat(sys.diff)};
J = sys.L;

rf_function = Function('rf_function',{sys.x,sys.u,sys.d},...
                    {vertcat(G{:}),J,vertcat(sys.x(1),sys.x(6))},{'x','u','d'},{'g','J','c'});
rf  = rootfinder('rf','newton',rf_function);

Ju = rf.factory('Ju',{'x','u','d'},{'jac:J:u'});
Ju1 = Function('Ju1',{sys.x,sys.u,sys.d},{Ju(sys.x,sys.u,sys.d)},{'x','u','d'},{'Ju'});

Gu = rf.factory('gu',{'x','u','d'},{'jac:g:u'});
Gd = rf.factory('gd',{'x','u','d'},{'jac:g:d'});
Juu = rf.factory('Juu',{'x','u','d'},{'hess:J:u:u'});
Jud = Ju1.factory('Jud',{'x','u','d'},{'jac:Ju:d'});

Juu = full(Juu(x_opt,u_opt,d_val));
Jud = full(Jud(x_opt,u_opt,d_val));
Gu = full(Gu(x_opt,u_opt,d_val));
Gd = full(Gd(x_opt,u_opt,d_val));

end
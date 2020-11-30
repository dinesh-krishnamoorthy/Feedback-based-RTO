clc
clear
import casadi.*

Ts = 300; % Sample time

% Load William-Otto reactor model
[sys,par,F] = WilliamOtto(Ts);
d_val = 1.3; % Disturbance Fa

% Initialize WO reactor
u_in = 3;
[xf,exitflag] = solveODE(sys,par,d_val,u_in);

%% configure

T = 150;  
f = 1/T;
nEst = 1800;
k = 0.00000005;
%% Simulation

nIter = 8*3600;
h = waitbar(0,'Simulation in Progress...');

GC.u = u_in(1);
GC.u0 = u_in(1);
GC.err = 0;
GC.err0 = 0;
GC.uf = GC.u;
 u_in = [3];
for sim_k = 1:nIter
    waitbar(sim_k /nIter)

    sim.u(sim_k) = u_in + 0.01*sin(2*pi*f*(sim_k));
    Fk = F('x0',xf,'p',vertcat(d_val, sim.u(sim_k) ));
    xf = full(Fk.xf);
    
    % Gradient Estimator
        
    J(sim_k) = full(Fk.qf); % Cost function
    
    if sim_k>3*nEst
        Y = J(sim_k-nEst:sim_k-1)';
        U = [sim.u(sim_k-nEst:sim_k-1);ones(1,nEst)]';
        theta = (U'*U)\(U'*Y);  
        Ju_hat = theta(1);
    else
        
        Ju_hat = 0;
    end
    
    u_in = u_in - k*Ju_hat;
    
    
    sim.x(:,sim_k) = xf;
    
    sim.Ju(:,sim_k) = Ju_hat;
    
end
close(h)
%%

figure(12)
clf
subplot(221)
plot((sim.x(6,:)))
ylabel('x_G')
grid on
subplot(222)
hold all
plot(sim.Ju./300)
ylabel('J_u')
grid on
subplot(224)
plot((sim.u(1,:)))
ylabel('F_B')
grid on
subplot(223)
plot((sim.x(7,:)))
ylabel('T_r')
grid on

% figure(15)
% plot((sim.u(1,600:end)-u_opt(1)),(sim.x(7,600:end)-x_opt(7)),'.')

save('sim_LSE','sim')



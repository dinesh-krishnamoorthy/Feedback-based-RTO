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
k = 0.0009;
%% Simulation

nIter = 8*3600;
h = waitbar(0,'Simulation in Progress...');

u_in = [3];

A = eye(2);
JQ = 1e1*eye(2);
JR = 1e-6*eye(2);
x_hat = [-30;10];
JPk = 10*eye(2);
Ju_hat = 0;
for sim_k = 1:nIter
    waitbar(sim_k /nIter)
   
    sim.u(sim_k) = u_in + 0.01*sin(2*pi*f*(sim_k));
    Fk = F('x0',xf,'p',vertcat(d_val, sim.u(sim_k) ));
    xf = full(Fk.xf);
    sim.J(sim_k) = full(Fk.qf);
    
    % Gradient Estimator
    
    if sim_k>5 && rem(sim_k,5)==0
        Yk = [sim.J(sim_k);sim.J(sim_k-1)];
        Ck = [sim.u(sim_k),1;sim.u(sim_k-1),1];
        
        % Update
        xk_1 = A*x_hat;
        JPk_1 = A*JPk*A' + JQ;
        % Predict
        ek = Yk - Ck*xk_1;
        Kk = JPk_1*Ck'*(Ck*JPk_1*Ck' + JR)^-1;
        x_hat = xk_1 + Kk*ek;
        JPk = (eye(2)-Kk*Ck)*JPk_1;
        
        Ju_hat = x_hat(1)*1e-4;
    
    end
    
    u_in = u_in - 0.0009*Ju_hat;
   
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
plot(sim.Ju*10)
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

save('sim_KF','sim')

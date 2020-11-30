clc
clear
import casadi.*

Ts = 300; % Sample time

% Load William-Otto reactor simulator
[sys,par,F] = WilliamOtto(Ts);
d_val = 1.3; % Disturbance Fa

% Initialize WO reactor
u_in = 3;
[xf,exitflag] = solveODE(sys,par,d_val,u_in);

%%
% Configure Extremum seeking control
a = 0.1;
k = 0.007;
T = 150;         % dither time period
Th = 1200*150;  % High Pass Filter
Tl = 300*150;       % Low Pass filter
f = 1/T;

al = Ts/(Ts + Tl);
ah = Th/(Ts + Th);

Ju_hat = 0; 

%% Simulation

nIter = 8*3600;
h = waitbar(0,'Simulation in Progress...');

u_in = [3];


for sim_k = 1:nIter
    waitbar(sim_k /nIter)
    
    sim.u(sim_k) = u_in + 0.01*sin(2*pi*f*(sim_k));
    Fk = F('x0',xf,'p',vertcat(d_val, sim.u(sim_k) ));
    xf = full(Fk.xf);
    
    % Gradient Estimator
    
    J(sim_k) = full(Fk.qf); % Cost function
    
    % Discrete ESC algorithm
    if sim_k>150
        z(sim_k) = ah*z(sim_k-1) + ah*J(sim_k) - ah*J(sim_k-1);                 % High pass filter
        x(sim_k) = ((1-al)*x(sim_k-1) + al*z(sim_k)*sin(2*pi*f*(sim_k)));  % correlation and Low pass filter
        Ju_hat = x(sim_k)/Ts;
       
    else
        u_hat(sim_k) = u_in;
        z(sim_k) = 0;
        x(sim_k) = 0;
        u(sim_k) = 0;
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
plot(2*sim.Ju./0.01)
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


save('sim_ESC','sim')


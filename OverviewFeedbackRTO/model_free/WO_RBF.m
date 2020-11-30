clc
clear
import casadi.*

Ts = 1; % Sample time

% Load William-Otto reactor simulator
[plant,par,F] = WilliamOtto(Ts);
d_val = 1.3; % Disturbance Fa

% Initialize WO reactor
u_in = 3;
[xf,exitflag] = solveODE(plant,par,d_val,u_in);

% Load William-Otto reactor model
[model,par] = WilliamOtto2reaction;

[Ju,Gu] = Jacobian(model);
% Get RBF weights offline
nm = 100;
nn = 8;
b = 1;
[W,nn,mu] = WO_TrainRBF(Ts,nn,b,0);

%% Simulation

nIter = 4*3600;
h = waitbar(0,'Simulation in Progress...');

GC.u = u_in(1);
GC.u0 = u_in(1);
GC.err = 0;
GC.err0 = 0;
GC.uf = GC.u;

Ju_hat = 0; 

for sim_k = 1:nIter
    waitbar(sim_k /nIter)
    
    % Gradient Controller
    GC.tauC = 200;
    GC.Kp = 300/(37*(GC.tauC+300)); %3.996 % theta = 300?
    GC.Ki = 1.5*GC.Kp/(max(300,4*GC.tauC));
    
    GC.err = (0 - Ju_hat);
    GC.u = GC.u + GC.Kp*GC.err + GC.Ki*GC.err - GC.Kp*GC.err0;
    GC.uf = GC.u;
    GC.err0 = GC.err;
    
    
    if sim_k>600
        u_in = vertcat(GC.uf);
    end
    
    Fk = F('x0',xf,'p',vertcat(d_val,u_in));
    xf = full(Fk.xf);
    xmodel = xf([1,2,4,5,6,7]);
    sim.J(sim_k) = full(Fk.qf);
    sim.u(:,sim_k) = u_in;
    
    Ju_hat_model = full(Ju(xmodel,u_in,d_val));
    
    %% estimate error gradient
    
    if sim_k >nm

        % Compute gradient of the RBF network
        Wp = W(2:end)';
        rSum = 0;
        for i = 1:nn
            rSum = rSum + exp(-b*(u_in-mu(i))^2);
        end
        
        for i = 1:nn
            r(i) = exp(-b*(u_in-mu(i))^2)/rSum;
        end
        
        rK = 0;
        for i = 1:nn
            rK = rK + (u_in - mu(i))*r(i);
        end
        for i = 1:nn
            dr(i) = -2*b*r(i)*((u_in-mu(i)) - rK);
        end
        
        derr = Wp*dr';
        
        % Compute plant grad. = model grad. + RBF grad.
        Ju_hat = Ju_hat_model + derr ;
        
    else
        Ju_hat = Ju_hat_model;
    end
    
    sim.x(:,sim_k) = xmodel;
    sim.c(sim_k) = Ju_hat;
    
    
end
close(h)
%%

figure(12)
clf
subplot(221)
hold all
plot((sim.x(5,:)))
ylabel('x_G')
grid on
subplot(222)
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
plot((sim.x(6,:)))
ylabel('T_r')
grid on




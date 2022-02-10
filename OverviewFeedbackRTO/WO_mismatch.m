clc
clear
import casadi.*

FileName = mfilename('fullpath');
[directory,~,~] = fileparts(FileName);
[parent,~,~] = fileparts(directory);
addpath([directory '/data'])
addpath([directory '/models'])
addpath([directory '/functions'])


Ts = 1; % Sample time

% Load William-Otto reactor simulator  - plant
[plant,par,F] = WilliamOtto(Ts);
d_val = 1.3; % Disturbance Fa

[xopt,uopt,sol] = SSOpt(plant,par,d_val);

% Initialize WO reactor
u_in = 3;
[xf,exitflag] = solveODE(plant,par,d_val,u_in);

% Load William-Otto reactor model
[model,par] = WilliamOtto2reaction;
[xopt1,uopt1,sol] = SSOpt(model,par,d_val);
[Ju,Gu] = Jacobian(model);
[Jup,Gu] = Jacobian(plant);

%% Configure EKF

EKF = prepareEKF(model,Ts);

ny = numel(model.y);
nxEKF = numel(model.x);
xk_hat = xf([1,2,4,5,6,7]);
uEKF = u_in;
Pk = 1e3.*eye(nxEKF);
Qk = 1e3.*eye(nxEKF);
Rk = 1e0.*eye(ny);
K = 0;


%% Simulation

nIter = 4*3600;
h = waitbar(0,'Simulation in Progress...');

GC.u = u_in(1);
GC.u0 = u_in(1);
GC.err = 0;
GC.err0 = 0;
GC.uf = GC.u;

Ju_hat = 0;
Ju_sp = 0;


correction = 1;
A = eye(2);
JQ = 1e1*eye(2);
JR = 1e-6*eye(2);
x_hat = [-30;10];
JPk = 10*eye(2);


Ju_hat_plant = 0;

for sim_k = 1:nIter
    waitbar(sim_k /nIter)
    
    % Gradient Controller
    GC.tauC = 200;
    GC.Kp = 300/(37*(GC.tauC+300)); %3.996 % theta = 300?
    GC.Ki = 1.5*GC.Kp/(max(300,4*GC.tauC));
    
    
    GC.err = (Ju_sp - Ju_hat);
    GC.u = GC.u + GC.Kp*GC.err + GC.Ki*GC.err - GC.Kp*GC.err0;
    GC.uf = GC.u;
    GC.err0 = GC.err;
    
    
    if sim_k>600
        u_in = max(0,GC.uf);
    end
    
    Fk = F('x0',xf,'p',vertcat(d_val,u_in));
    xf = full(Fk.xf);
    xmodel = xf([1,2,4,5,6,7]);
    sim.J(sim_k) = full(Fk.qf);
    sim.u(:,sim_k) = u_in;
    
    Ju_hat_model = full(Ju(xmodel,u_in,d_val));
    % Gradient Estimator
    [xk_hat,Pk] = EKF_estimation(EKF,xmodel,xk_hat,u_in,Pk,Qk,Rk,d_val);
    Ju_hat = EstJu(EKF,xk_hat,u_in,d_val);
    
    sim.x(:,sim_k) = xf;
    sim.c(sim_k) = Ju_hat;
    
    
    %%
    if correction
        %         Ju_hat_plant = full(Jup(xf,u_in,d_val));
        if sim_k>2
            du = sim.u(:,sim_k-1) - sim.u(:,sim_k);
            if abs(du)>0
                Ju_hat_plant = (sim.J(sim_k-1) - sim.J(sim_k)/du)
            end
        end
        
        sim.Jup(sim_k) = Ju_hat_plant;
        Ju_sp = -Ju_hat_plant + Ju_hat;
        
    end
    
end
close(h)
%%

figure(12)
clf
subplot(221)
hold all
plot((sim.x(6,:)))
ylabel('x_G')
grid on
subplot(222)
hold all
plot(sim.c)

ylabel('J_u')
grid on
subplot(224)
hold all
plot((sim.u(1,:)))
plot(uopt1.*ones(size(sim.u(1,:))),':')
plot(uopt.*ones(size(sim.u(1,:))),'--')

ylabel('F_B')
grid on
subplot(223)
hold all
plot((sim.x(7,:)))
plot(xopt1(6).*ones(size(sim.u(1,:))),':')
plot(xopt(7).*ones(size(sim.u(1,:))),'--')
ylabel('T_r')
grid on
%%

% save('sim_mismatch','sim')

function Ju_hat = EstJu(EKF,x_hat,uEKF,d_hat)

A = full(EKF.JacAx(x_hat,uEKF,d_hat));
B = full(EKF.JacBu(x_hat,uEKF,d_hat));
C = full(EKF.JacJx(x_hat,uEKF,d_hat));
D = full(EKF.JacJu(x_hat,uEKF,d_hat));
Ju_hat = -C*(A\B) + D;

end

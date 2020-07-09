import casadi.*

clear
clc

Ts = 1; % sampling time in min

global pF pB pD
pF = 10;             % feed price [$/mol]
pB = 8;              % Bottom product price [$/mol]
pD = 86;             % distillate price [$/mol]

[sys,F_integrator,f,par] = cola_lv(Ts);

% ---------- initialization -----------

par.F = 1.42;   % Feedrate
par.pV = 0.015; % Energy price

L = 3.1;  
V = 3.9;

dist = vertcat(par.F,par.pV);
u_in = vertcat(L,V);
[xf,exitflag] = solveODE(sys,dist,u_in);

% ---------- Setup PI Controllers ------ 

% Concentration controller for xD tuned using SIMC rule
CCd.k = 0.456; CCd.tau1 = 36; CCd.tauC = 10;
CCd.Kp = CCd.tau1/(CCd.k*CCd.tauC);
CCd.Ki = CCd.Kp/CCd.tau1;
CCd.err = 0.95-xf(41);    % initialize error
CCd.err0 = CCd.err;

% Concentration controller for xB tuned using SIMC rule
CCb.k = 0.813; CCb.tau1 = 18; CCb.tauC = 10;
CCb.Kp = CCb.tau1/(CCb.k*CCb.tauC);
CCb.Ki = CCb.Kp/CCb.tau1;
CCb.err = max(0.99,0.99127+0.41077*par.pV-30*par.pV^2)-(1-xf(1));
CCb.err0 = CCb.err;
CCb.V = V;

% -------------- Simulation -------------

for sim_k = 1:10*24*60
    
    % ----- Disturbance in pV ------
        if sim_k>3*24*60 && sim_k<=7*24*60
            par.pV = 0.008;
        else if sim_k>7*24*60
                par.pV = 0.02;
            end
        end
    
    % -------- controllers ---------
        % Concentration controller for xD (always active)
        CCd.err = 0.95-xf(41);
        L = L + CCd.Kp*CCd.err + CCd.Ki*CCd.err -CCd.Kp*CCd.err0;
        CCd.err0 = CCd.err;
    
    
        % Conceentration controller for xB
        CCb.err = max(0.99,0.99127+0.41077*par.pV-30*par.pV^2)-(1-xf(1)); % Max selector
        CCb.V = CCb.V + CCb.Kp*CCb.err + CCb.Ki*CCb.err -CCb.Kp*CCb.err0;
        V = min(4.008,CCb.V);  % Min selector
        CCb.V = V;  % Output from the max-min structure
        CCb.err0 = CCb.err;
    
    % --------- Simulator ----------
        par.d = vertcat(par.F,par.pV);
        Fk = F_integrator('x0',xf,'p',vertcat(L,V,par.d));
        % set new initial values for the next iteration
        xf =  full(Fk.xf); % states
        qf = full(Fk.qf); % cost function
    
    % ----- Extract and store data -----
        sim.L(sim_k)    = L;
        sim.V(sim_k)    = V;
        sim.B(sim_k)    = par.Bs + (xf(par.NT+1) - par.MBs)*par.KcB;
        sim.D(sim_k)    = par.Ds + (xf(2*par.NT) - par.MDs)*par.KcD; 
        sim.xD(sim_k)   = xf(41);
        sim.xB(sim_k)   = 1-xf(1);
        sim.q(sim_k)    = qf;
        sim.pV(sim_k)   = par.pV;
    
end

sim.t_days = (1:10*24*60)/(24*60);
sim.t_hrs = (1:10*24*60)/(60);

% save('sim','sim')
%%

figure(12)
clf
subplot(322)
hold all
plot(sim.t_days,0.95.*ones(size(sim.xD)),'k:','linewidth',1.5)
plot(sim.t_days,sim.xD,'r','linewidth',2)
ylabel('$CV_1: x_D$','Interpreter','Latex')
ylim([0.948,0.952])
box on
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
grid on
xticks([0:2:10])

subplot(324)
hold all
plot(sim.t_days,0.99.*ones(size(sim.xD)),'k:','linewidth',1.5)
plot(sim.t_days,0.99127+0.41077*sim.pV-30*sim.pV.^2,'k--','linewidth',0.75)
plot(sim.t_days,sim.xB,'r','linewidth',2)
ylabel('$CV_2: x_B$','Interpreter','Latex')
ylim([0.987,0.993])
box on
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
xticks([0:2:10])
grid on

subplot(321)
plot(sim.t_days,sim.L,'b','linewidth',2)
ylabel('$MV_1: L$','Interpreter','Latex')
ylim([3,3.3])
box on
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
xticks([0:2:10])
grid on

subplot(323)
hold on
plot(sim.t_days,4.008.*ones(size(sim.xD)),'k:','linewidth',1.5)
plot(sim.t_days,sim.V,'b','linewidth',2)
ylabel('$MV_2: V$','Interpreter','Latex')
ylim([3.85,4.1])
box on
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
xticks([0:2:10])
grid on

subplot(325)
plot(sim.t_days,par.F.*ones(size(sim.D)),'k','linewidth',2)
ylabel('$F$','Interpreter','Latex')
xlabel('time [days]','Interpreter','Latex')
box on
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
xticks([0:2:10])
grid on

subplot(326)
plot(sim.t_days,sim.pV,'k','linewidth',2)
ylabel('$d: pV$','Interpreter','Latex')
xlabel('time [days]','Interpreter','Latex')
box on
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
xticks([0:2:10])
grid on

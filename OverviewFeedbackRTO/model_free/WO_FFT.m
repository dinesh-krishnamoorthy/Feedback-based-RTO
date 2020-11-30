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

%% configure

 nL = 4;
L = 4096/nL;
Fs = 4096;
T = 1/Fs;
k=0.000005;
%% Simulation

nIter = 8*3600;
h = waitbar(0,'Simulation in Progress...');
 
for sim_k = 1:nIter
    waitbar(sim_k /nIter)

    sim.sine(sim_k) = 0.01*sin(2*pi*20*T*sim_k);
    sim.u(sim_k) = u_in + sim.sine(sim_k);
    Fk = F('x0',xf,'p',vertcat(d_val, sim.u(sim_k) ));
    xf = full(Fk.xf);
    
   
    % Gradient Estimator
        
    J(sim_k) = full(Fk.qf); % Cost function
    
    
    if sim_k>nL*L
        yL = J(sim_k-nL*L:sim_k-1)';
        u1L = sim.sine(sim_k-nL*L:sim_k-1);
        
        y = yL - movingavg(yL,L/2,0);
        u10 = u1L - movingavg(u1L,L/2,0);
        
        NFFT = 2^nextpow2(nL*L); % Next power of 2 from length of y
        f1 = Fs/2*linspace(0,1,NFFT/2+1);
        i1 = find(f1==20);
        
        Y = fft(y,NFFT)/NFFT;
        U1 = fft(u10,NFFT)/NFFT;
        
        Ymag = 2*abs(Y(1:NFFT/2+1));
        Yph = angle(Y)*180/pi; %phase information
        
        U1mag = 2*abs(U1(1:NFFT/2+1));
        U1ph = angle(U1)*180/pi; %phase information
        
        sim.yph(sim_k) = Yph(i1);
        sim.U1ph(sim_k) = U1ph(i1);
        
        Ju_hat = sign(sim.yph(sim_k)).*sign(sim.U1ph(sim_k))*Ymag(i1);
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
plot(sim.Ju)
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

save('sim_FFT','sim')

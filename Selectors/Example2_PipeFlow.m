clear
clc

% Define constraints
F_max = 10;
p1_min = 1.5e5;
p1_max = 2.5e5;
z1_max = 1;     % MV constraint

%% Tune PI controllers using SIMC rules

% Flow control for F_max
k_F = 8.6429; tau1 = 10;
tauC = 5; 
FC.K = 1/(k_F*tauC);
FC.Kp = tau1/(k_F*tauC);    % proportional gain
FC.Ki = FC.Kp/tau1;         % integral gain 
FC.Kaw =   FC.Ki/FC.Kp;     % Antiwindup gain

% pressure control for p1_max and p1_min
k_p1 = 1.8033e5; tau1 = 10;
tauC = 5; 
PC.k = 1/(k_p1*tauC);       
PC.Kp = tau1/(k_p1*tauC);   % proportional gain
PC.Ki = PC.Kp/tau1;         % integral gain
PC.Kaw =   PC.Ki/PC.Kp;     % Antiwindup gain


% Initialize simulation
F(1:2) = 10; FC.z(1:2) = 0.5;  FC.e0 = 0;
p1(1:2) = 2e5; PC.z(1:2) = 0.5; PC1.e0 = 0;
PC2.z(1:2) = 0.5; PC2.e0 = 0;
z(1:2) = 0.5;
y(1,1:2) = F(1:2); y(2,1:2) = p1(1:2);


% Disturbance profile
p2 = [1e5.*ones(300,1);...
    0.2e5.*ones(300,1);...
    1.25e5.*ones(300,1);...
    1.5e5.*ones(300,1);...
    1.75e5.*ones(300,1);...
    2e5.*ones(300,1);...
    2e5.*ones(300,1);...
    2e5.*ones(300,1)...
    ];
p0 = [3e5*ones(1800,1);...
    2.5e5*ones(300,1);...
    2.25e5*ones(300,1)];


for i = 2:length(p2)
    
    F(i)= y(1,i);
    p1(i) = y(2,i);
    
    FC.z(i+1) = FC.z(i)+ FC.Kp*((F_max-F(i))-FC.e0)+ FC.Ki*(F_max-F(i)) + FC.Kaw*(z(i) - FC.z(i));
    PC.z(i+1) = PC.z(i) + PC.Kp*(p1_max-p1(i) - PC1.e0) + PC.Ki*(p1_max-p1(i)) + PC.Kaw*(z(i) - PC.z(i));
    PC2.z(i+1) = max(0,PC2.z(i) + PC.Kp*((p1_min-p1(i))-PC2.e0) + PC.Ki*(p1_min-p1(i)) + PC.Kaw*(z(i) - PC2.z(i)));
    z(i+1) = max(0,min(1,max(PC2.z(i+1),min(FC.z(i+1),PC.z(i+1)))));
  
    FC.e0 = (F_max-F(i));
    PC1.e0 = (p1_max-p1(i));
    PC2.e0 = (p1_min-p1(i));
     
    [t,x] = ode45(@(t,y) ODE(t,y,z(i),p2(i),p0(i)),[0 1],y(:,i));
    y(:,i+1) = x(end,:) + 0.*[0.1,1000].*randn(1,2);
end


%%

figure(1)
clf
subplot(411)
set(gca,'FontSize',14)
hold all
plot(F,'k','linewidth',2.5)
plot(F_max*ones(size(F)),'k:','linewidth',1.5)
ylabel('$CV_1: F$ [kg/s]','Interpreter','latex')
ylim([2,12])
xlim([0,i])
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'Latex';

subplot(412)
set(gca,'FontSize',14)
hold all
plot(p1.*1e-5,'k','linewidth',2.5)
plot(1.5*ones(size(p1)),'k:','linewidth',1.5)
plot(2.5*ones(size(p1)),'k:','linewidth',1.5)
ylabel('$CV_2: p_1$ [bar]','Interpreter','latex')
ylim([1,3])
xlim([0,i])
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'Latex';

subplot(413)
set(gca,'FontSize',14)
hold all
plot(z,'k','linewidth',2.5)
plot( PC.z(1:end),'b--','linewidth',1.5)
plot( PC2.z(1:end),'r--','linewidth',1)
plot( FC.z(1:end),'--','linewidth',1.5,'color',[0.1,0.6,1])
plot(1*ones(size(p1)),'k:','linewidth',1.5)
legend('$u$','$u_1$','$u_2$','$u_3$','Interpreter','latex',...
    'box','off','location','best','orientation','horizontal')
ylabel('$MV: z_1$ ','Interpreter','latex')
ylim([0,2])
xlim([0,i])
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'Latex';

subplot(414)
set(gca,'FontSize',14)
hold all
plot(p0*1e-5,'linewidth',2.5,'color',[0.3,0.3,0.3])
plot(p2*1e-5,'linewidth',2.5,'color',[0.6,0.6,0.6])
legend('$p_0$','$p_2$','Interpreter','latex','box','off','location','best')
ylabel('$d: [p_0,p_2]$ [bar] ','Interpreter','latex')
ylim([0,3.5])
xlim([0,i])
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'Latex';
xlabel('time unit','Interpreter','latex')



function dydt = ODE(t,y,z,p2,p0,z2)

if nargin<5
    p0 = 3e5;
    z2 = 1;
end

if nargin<6
    z2 = 1;
end

Cv1 = 2e-3;
Cv2 = 1e-3;
rho = 1000;


a = (Cv1*z/(Cv2*z2))^2;
p1 = ((a*p0 + p2)/(a+1));
F1 = Cv2*1*sqrt((p1-p2)/rho);
F =  F1*rho;

tau = 10;
dydt = [(F-y(1))/tau; (p1-y(2))/tau];

end


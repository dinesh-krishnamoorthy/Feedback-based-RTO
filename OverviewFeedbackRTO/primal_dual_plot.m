sim1 = load('PrimalDual11.mat');
sim2 = load('selector11.mat');
time = (1:9*3600)/3600;

for i = 1:9*3600
    if i<=3*3600
        sim.d(i) = 1.4;
    end
    if i>3*3600
        sim.d(i) = 1.9;
    end
    if i>6*3600
        sim.d(i) = 1.5;
    end
end
% 
% for i = 1:4*3600
%    
%         sim.d(i) = 1.4;
%         
% %         if i >=2*3600
% %             sim1.sim.u(2,i) = 3;
% %         end
%     
% end


sim1.sim.J = -(1043.38.*sim1.sim.x(4,:).*(sim.d+sim1.sim.u(1,:))+...
    20.92*sim1.sim.x(5,:).*(sim.d+sim1.sim.u(1,:)) ...
    - 79.23*sim.d ...
    - 118.34*sim1.sim.u(1,:));

sim2.sim.J = -(1043.38*sim2.sim.x(4,:).*(sim.d+sim2.sim.u(1,:))+...
    20.92*sim2.sim.x(5,:).*(sim.d+sim2.sim.u(1,:)) ...
    - 79.23*sim.d ...
    - 118.34*sim2.sim.u(1,:));

figure(12)
clf
hold all
subplot(421)
hold all
plot(time,0.08.*ones(size(time)),':','linewidth',1.5)
plot(time,(sim1.sim.x(6,:)),'k','linewidth',1.5)
plot(time,(sim2.sim.x(6,:)),'r--','linewidth',1.5)
ylabel('$x_G$ [kg/kg]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
% ylim([0.075,0.085])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on


subplot(422)
hold all
plot(time,0.12.*ones(size(time)),':','linewidth',1.5)
plot(time,(sim1.sim.x(1,:)),'k','linewidth',1.5)
plot(time,(sim2.sim.x(1,:)),'r--','linewidth',1.5)
ylabel('$x_A$ [kg/kg]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
% ylim([0.11,0.14])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on


subplot(424)
hold all
plot(time,(sim1.sim.u(1,:)),'k','linewidth',1.5)
plot(time,(sim2.sim.u(1,:)),'r--','linewidth',1.5)
ylabel('Primal $MV_2:F_B$ [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
legend('Primal-Dual control','Active constraint control with selectors','Interpreter','latex','orientation','horizontal')
xlim([0,9])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on

subplot(423)
hold all
plot(time,(sim1.sim.u(2,:))-273,'k','linewidth',1.5)
plot(time,(sim2.sim.x(7,:))-273,'r--','linewidth',1.5)
ylabel('Primal $MV_1:T_r$ [$^\circ$ C]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on


subplot(426)
hold all
plot(time,(sim1.sim.lam(1,:)),'k','linewidth',1.5)
ylabel(' Dual $MV_1: \lambda_{x_A}$','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on

subplot(425)
hold all
plot(time,(sim1.sim.lam(2,:)),'k','linewidth',1.5)
ylabel(' Dual $MV_2:\lambda_{x_G}$','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on

subplot(427)
hold all
plot(time,(sim1.sim.J),'k','linewidth',1.5)
plot(time,(sim2.sim.J),'r--','linewidth',1.5)
ylabel(' Cost $J$ [\$/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on

subplot(428)
hold all
plot(time,(sim.d),'linewidth',1.5)
ylabel(' $F_a$ [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,9])
xticks(0:1:9)
xlim([0,9])
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on

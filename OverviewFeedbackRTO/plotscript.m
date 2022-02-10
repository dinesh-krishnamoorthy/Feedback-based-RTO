SOC = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_SOC2.mat');
NEC = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_NEC.mat');
TSA = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_TSA2.mat');
FRTO = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_FRTO2.mat');


time = (1:4*3600)/3600;

for i = 1:4*3600
    if i<=2*3600
        sim.d(i) = 1.3;
        sim.J(i) = -76.7095;
        sim.u(i) = 3.1943;
        sim.T(i) = 349.7717;
    end
    if i>2*3600
        sim.d(i) = 1;
        sim.J(i) = -73.0242;
        sim.u(i) = 2.5506;
        sim.T(i) = 347.5048;
    end
end


SOC.sim.J = -(1043.38.*SOC.sim.x(4,:).*(sim.d+SOC.sim.u(1,:))+...
    20.92*SOC.sim.x(5,:).*(sim.d+SOC.sim.u(1,:)) ...
    - 79.23*sim.d ...
    - 118.34*SOC.sim.u(1,:));

NEC.sim.J = -(1043.38*NEC.sim.x(4,:).*(sim.d+NEC.sim.u(1,:))+...
    20.92*NEC.sim.x(5,:).*(sim.d+NEC.sim.u(1,:)) ...
    - 79.23*sim.d ...
    - 118.34*NEC.sim.u(1,:));

TSA.sim.J = -(1043.38.*TSA.sim.x(4,:).*(sim.d+TSA.sim.u(1,:))+...
    20.92*TSA.sim.x(5,:).*(sim.d+TSA.sim.u(1,:)) ...
    - 79.23*sim.d ...
    - 118.34*TSA.sim.u(1,:));

FRTO.sim.J = -(1043.38*FRTO.sim.x(4,:).*(sim.d+FRTO.sim.u(1,:))+...
    20.92*FRTO.sim.x(5,:).*(sim.d+FRTO.sim.u(1,:)) ...
    - 79.23*sim.d ...
    - 118.34*FRTO.sim.u(1,:));



figure(12)
clf
hold all
subplot(221)
hold all
plot(time,0.08.*ones(size(time)),'k:','linewidth',1)
plot(time,(NEC.sim.x(6,:)),'r--','linewidth',1.5)
plot(time,(SOC.sim.x(6,:)),'r','linewidth',2)
plot(time,(TSA.sim.x(6,:)),'k--','linewidth',1.5)
plot(time,(FRTO.sim.x(6,:)),'k','linewidth',1.5)
ylabel('$x_G$ [kg/kg]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,4])
xticks(0:1:4)
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on


subplot(222)
hold all
plot(time,sim.J,'k:','linewidth',1)
plot(time,(NEC.sim.J),'r--','linewidth',1.5)
plot(time,(SOC.sim.J),'r','linewidth',2)
plot(time,(TSA.sim.J),'k--','linewidth',1.5)
plot(time,(FRTO.sim.J),'k','linewidth',1.5)
ylabel('Cost $J$ [\$/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
legend('Ideal','NEC','Nullspace','Two-step','FRTO','Interpreter','latex','location','best','box','off')

xlim([0,4])
xticks(0:1:4)
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on


subplot(224)
hold all
plot(time,sim.u,'k:','linewidth',1)
plot(time,(NEC.sim.u(1,:)),'r--','linewidth',1.5)
plot(time,(SOC.sim.u(1,:)),'r','linewidth',2)
plot(time,(TSA.sim.u(1,:)),'k--','linewidth',1.5)
plot(time,(FRTO.sim.u(1,:)),'k','linewidth',1.5)
ylabel('Primal $MV_2:F_B$ [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
legend('Ideal','NEC','Nullspace','Two-step','FRTO','Interpreter','latex','location','best','box','off')
xlim([0,4])
xticks(0:1:4)
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on

subplot(223)
hold all
plot(time,sim.T,'k:','linewidth',1)
plot(time,(NEC.sim.x(7,:)),'r--','linewidth',1.5)
plot(time,(SOC.sim.x(7,:)),'r','linewidth',2)
plot(time,(TSA.sim.x(7,:)),'k--','linewidth',1.5)
plot(time,(FRTO.sim.x(7,:)),'k','linewidth',1.5)
ylabel('Primal $MV_1:T_r$ [$^\circ$ C]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,4])
xticks(0:1:4)
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on



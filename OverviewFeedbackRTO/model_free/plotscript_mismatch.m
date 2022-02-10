mismatch = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_mismatch.mat');
RBF = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_RBF.mat');
FMA = load('/Users/dineshk/Dropbox/Matlab/Feedback-based-RTO/OverviewFeedbackRTO/data/sim_mismatch2.mat');


time = (1:4*3600)/3600;

for i = 1:4*3600
        sim.d(i) = 1.3;
        
        sim.J(i) = -76.7095;
        sim.u(i) = 3.1943;
        sim.T(i) = 349.7717;
        
        sim.Jm(i) = -44.4596;
        sim.um(i) = 3.3130;
        sim.Tm(i) = 352.4098;
end



figure(12)
clf
hold all
subplot(221)
hold all
plot(time,0.08.*ones(size(time)),'k:','linewidth',1)
plot(time,(mismatch.sim.x(6,:)),'r','linewidth',1)
plot(time,(FMA.sim.x(6,:)),'k-','linewidth',1.5)
plot(time,(RBF.sim.x(6,:)),'r--','linewidth',1.5)
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
plot(time,(mismatch.sim.J),'r','linewidth',1)
plot(time,(FMA.sim.J),'k-','linewidth',1.5)
plot(time,(RBF.sim.J),'r--','linewidth',1.5)
ylabel('Cost $J$ [\$/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
legend('Ideal','mismatch','FMA','RBF','Interpreter','latex','location','best','box','off')

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
plot(time,(mismatch.sim.u(1,:)),'r','linewidth',1)
plot(time,(FMA.sim.u(1,:)),'k-','linewidth',1.5)
plot(time,(RBF.sim.u(1,:)),'r--','linewidth',1.5)
ylabel('Primal $MV_2:F_B$ [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
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
plot(time,(mismatch.sim.x(7,:)),'r','linewidth',1)
plot(time,(FMA.sim.x(7,:)),'k-','linewidth',1.5)
plot(time,(RBF.sim.x(7,:)),'r--','linewidth',1.5)
ylabel('Primal $MV_1:T_r$ [$^\circ$ C]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,4])
xticks(0:1:4)
axs = gca;
axs.TickLabelInterpreter = 'latex';
axs.FontSize = 14;
grid on
box on



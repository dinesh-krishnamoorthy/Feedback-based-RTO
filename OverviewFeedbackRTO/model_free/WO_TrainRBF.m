function [W,nn,mu] = WO_TrainRBF(Ts,nn,b,Validation)
import casadi.*

if nargin <5
    Validation  = 0;
end
% Load William-Otto reactor plant
[plant,plant_par,F] = WilliamOtto(Ts);
d_val = 1.3;

% Load William-Otto reactor plant
[model,model_par] = WilliamOtto2reaction;

h = waitbar(0,'Training RBF weights...');

for i = 1:200
     waitbar(i/200)
    u_in = (7-2)*rand(1,1)+2;
    % Compute plant steady-state
    [x_optP] = solveODE(plant,plant_par,d_val,u_in);
    J_plant = cost(x_optP(4),x_optP(5),d_val,u_in);
    % Compute model steady-state
    [x_opt] = solveODE(model,model_par,d_val,u_in);
    J_model = cost(x_opt(3),x_opt(4),d_val,u_in);
    
    
    sim.u(i) = u_in;
    sim.yP(:,i) = x_optP([1,2,4,5,6,7]);
    sim.yM(:,i) = x_opt;
    sim.err(:,i) = J_plant - J_model; %sim.yP(:,i)-sim.yM(:,i);
end
close(h)


nm = 200;
u = sim.u;

[idx,mu] = kmeans(u',nn);

for k = 1:nm
    rSum = 0;
    for i = 1:nn
        rSum = rSum + exp(-b*(u(k)-mu(i))^2);
    end
    for i =1:nn
        r(k,i) = exp(-b*(u(k)-mu(i))^2)/rSum ;
    end
end

R = [ones(nm,1),r];
Y = sim.err(1:nm)';

W = (R'*R)\(R'*Y);

%% Validate model

if Validation
    for k = 1:200
        rSum = 0;
        for i = 1:nn
            rSum = rSum + exp(-b*(u(k)-mu(i))^2);
        end
        for i =1:nn
            r(k,i) = exp(-b*(u(k)-mu(i))^2)/rSum ;
        end
    end
    
    R = [ones(200,1),r];
    Y_hat = R*W;
    
    
    figure(162)
    clf
    subplot(211)
    hold all
    plot(sim.err')
    plot(Y_hat,'--')
    subplot(212)
    plot((Y_hat-sim.err'))
end

function J = cost(xp,xe,Fa,Fb) 
J = -(1043.38*xp*(Fa+Fb)+20.92*xe*(Fa+Fb) - 79.23*Fa - 118.34*Fb);
end

end


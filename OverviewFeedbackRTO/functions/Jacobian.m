function [Ju,Gu] = Jacobian(sys)
% Function to evaluate the Jacobian of the cost and constraint
% Adetola & Guay (2007)
% William Otto example
%
% Written by: Dinesh Krishnamoorthy, Sep 2019

import casadi.*

J = sys.L; % Economic objective
G = {vertcat(sys.diff)};

rf_function = Function('rf_function',{sys.x,sys.u,sys.d},...
                    {vertcat(G{:}),J,sys.nlcon},{'x','u','d'},{'g','J','c'});
rf  = rootfinder('rf','newton',rf_function);

Ju = rf.factory('Ju',{'x','u','d'},{'jac:J:u'});
Gu = rf.factory('gu',{'x','u','d'},{'jac:c:u'});

end

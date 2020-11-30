function [xf,exitflag,sol] = solveODE(sys,par,d_val,u_in,opts)

% Function that computes the steady-state optimum
% Written by Dinesh Krishnamoorthy, Jul 2019, NTNU

import casadi.*

if nargin<8
    % Default Options
    opts = struct('warn_initial_bounds',false, ...
        'print_time',false, ...
        'ipopt',struct('print_level',1)...
        );
end

lbx = par.lbx;
ubx = par.ubx;
dx0 = par.dx0;

assert(numel(sys.u)==numel(u_in),'Error: Input dimension mismatch !')
assert(numel(sys.d)==numel(d_val),'Error: Disturbance dimension mismatch !')

w = {};
w0 = [];
lbw = [];
ubw = [];

g = {};
lbg = [];
ubg = [];

w = {w{:},sys.x};
lbw = [lbw;lbx];
ubw = [ubw;ubx];
w0 = [w0;dx0];

g = {g{:},vertcat(sys.diff)};
lbg = [lbg;zeros(numel(sys.diff),1)];
ubg = [ubg;zeros(numel(sys.diff),1)];

nlp = struct('x',vertcat(w{:}),'p',vertcat(sys.d,sys.u),'f',0,'g',vertcat(g{:}));
solver = nlpsol('solver','ipopt',nlp,opts);

sol = solver('x0',w0,'p',vertcat(d_val,u_in),'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
wf = full(sol.x);
xf = wf(1:numel(sys.x));


flag = solver.stats();
exitflag =  flag.return_status;

assert(flag.success == 1,'Error: Steady-state ODE solver failed !')


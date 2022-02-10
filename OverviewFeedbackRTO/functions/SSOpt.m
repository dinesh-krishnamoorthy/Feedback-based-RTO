function [xopt,uopt,exitflag,sol] = SSOpt(sys,par,d_val,opts)

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

lbu = par.lbu;
ubu = par.ubu;
u0 = par.u0;

assert(numel(sys.d)==numel(d_val),'Error: Disturbance dimension mismatch !')

w = {};
w0 = [];
lbw = [];
ubw = [];

g = {};
lbg = [];
ubg = [];

w = {w{:},sys.x,sys.u};
lbw = [lbw;lbx;lbu];
ubw = [ubw;ubx;ubu];
w0 = [w0;dx0;u0];

g = {g{:},vertcat(sys.diff)};
lbg = [lbg;zeros(numel(sys.diff),1)];
ubg = [ubg;zeros(numel(sys.diff),1)];

if ~isempty(sys.nlcon)
    assert(numel(sys.nlcon)==numel(sys.lb))
    assert(numel(sys.nlcon)==numel(sys.ub))
    
    g = {g{:},sys.nlcon};
    lbg = [lbg;sys.lb];
    ubg = [ubg;sys.ub];
end


nlp = struct('x',vertcat(w{:}),'p',sys.d,'f',sys.L,'g',vertcat(g{:}));
solver = nlpsol('solver','ipopt',nlp,opts);

sol = solver('x0',w0,'p',d_val,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
wf = full(sol.x);
xopt = wf(1:numel(sys.x));
uopt = wf(numel(sys.x)+1);

flag = solver.stats();
exitflag =  flag.return_status;

assert(flag.success == 1,'Error: Steady-state ODE solver failed !')


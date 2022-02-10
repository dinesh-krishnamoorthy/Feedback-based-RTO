function [EKF] = prepareEKF(sys,Ts)

import casadi.*

nx = numel(sys.x);

diff_EKF= [];


    for i = 1:nx
        diff_EKF = [diff_EKF; sys.x(i) + Ts*sys.diff(i)];
    end
    
    EKF.f_EKF = Function('f_EKF',{sys.x,sys.u,sys.d},{diff_EKF},{'x','u','d'},{'xdot'});
    EKF.JacFx = Function('JacFx',{sys.x,sys.u,sys.d},{jacobian(diff_EKF,sys.x)});
    
    EKF.h_EKF = Function('h_EKF',{sys.x,sys.u,sys.d},{sys.y});
    EKF.JacHx = Function('JacHx',{sys.x,sys.u,sys.d},{jacobian(sys.y,sys.x)});
    
    
    EKF.JacAx = Function('JacAx',{sys.x,sys.u,sys.d},{jacobian(sys.diff,sys.x)});
    EKF.JacBu = Function('JacBu',{sys.x,sys.u,sys.d},{jacobian(sys.diff,sys.u)});
    EKF.JacJx = Function('JacJx',{sys.x,sys.u,sys.d},{jacobian(sys.L,sys.x)});
    EKF.JacJu = Function('JacJu',{sys.x,sys.u,sys.d},{jacobian(sys.L,sys.u)});
    EKF.JacCx = Function('JacCx',{sys.x,sys.u,sys.d},{jacobian(sys.nlcon,sys.x)});
    EKF.JacCu = Function('JacCu',{sys.x,sys.u,sys.d},{jacobian(sys.nlcon,sys.u)});
    



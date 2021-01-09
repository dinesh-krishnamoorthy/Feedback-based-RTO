function [xk_hat,Pk] = EKF_estimation(EKF,ymeas,xk_hat,uEKF,Pk,Qk,Rk,d_val)

nxEKF = numel(xk_hat);

Fj = full(EKF.JacFx(xk_hat,uEKF,d_val));

xk_hat_1 =  full(EKF.f_EKF(xk_hat,uEKF,d_val));
Pk_1 = Fj*Pk*Fj' + Qk;

Hj = full(EKF.JacHx(xk_hat_1,uEKF,d_val));
ek = full(ymeas - EKF.h_EKF(xk_hat_1,uEKF,d_val));
Sk = Hj*Pk_1*Hj' + Rk;
Kk = (Pk_1*Hj')/(Sk);
xk_hat = xk_hat_1 + Kk*ek;
Pk = (eye(nxEKF) - Kk*Hj)*Pk_1;


end
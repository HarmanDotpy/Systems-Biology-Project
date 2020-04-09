function v = dynamic_LacI_AraC(t,y,Z,Parameter_Vec)
%% Parameters
%[alpha_r; alpha_a; gamma_r; gamma_a; f; Ca; Cr_T; N; beta; Ro; A_T] = Parameter_Vec;
alpha_r = Parameter_Vec(1);
alpha_a = Parameter_Vec(2);
gamma_r = Parameter_Vec(3);
gamma_a = Parameter_Vec(4);
f = Parameter_Vec(5);
Ca = Parameter_Vec(6);
Cr_T = Parameter_Vec(7);
N = Parameter_Vec(8);
beta = Parameter_Vec(9);
Ro = Parameter_Vec(10);
A_T = Parameter_Vec(11);
%%
lag_a = Z(:,1);
lag_r = Z(:,2);
%%
v = zeros(2,1);
%% LacI
% dr/dt = v(1), y(1) = r(t), lag_r(2) = a(t - lag_r), lag_r(1) = r(t - lag_r)
v(1) = (alpha_r*(inv(f) + lag_r(2)/Ca)/((1 + lag_r(2)/Ca)*(1+ lag_r(1)/Cr_T)^N)...
    -beta*y(1) - gamma_r*y(1)/(Ro + y(1) + y(2)))/A_T;
%% AraC
% da/dt =  v(2), y(2) = a(t), lag_a(2) = a(t - lag_a), lag_a(1) = r(t - lag_a)
v(2) = (alpha_a*(inv(f) + lag_a(2)/Ca)/((1 + lag_a(2)/Ca)*(1 + lag_a(1)/Cr_T)^N)...
    -beta*y(2) - gamma_a*y(2)/(Ro + y(1) + y(2)))/A_T;

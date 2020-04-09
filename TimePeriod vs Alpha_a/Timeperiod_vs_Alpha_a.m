%% Parameter Values
%maximal transcription rates
alpha_r =  265;
alpha_a = 92.75;
%transcriptional delay times
tau_r = 13.5;
tau_a = 15;
%dilution rate constant
beta = 0.0275;
%Michaelis–Menten constants
gamma_r = 76;
gamma_a = 76;
Ro = 1.8;
%Hill coefficient for LacI repression
N = 4;
%measure of the strength of the positive feedback loop
f = 2;
%concentration of AraC
Ca = 5;
%min and max values of Cr(T)
Cr_max = 830;
Cr_min = 50;
%reference temp
Tref = 37; %degree
T = 90;
%To is the temperature at which Cr(T) is half maximal
To = 38; %
%Hill Coefficient
b = 20;
%Arrhenius constant
theta = 4500; 

%% Arrhenius scaling term
A_T = exp(theta*(inv(T+273) - inv(Tref + 273))); 

%% delay values that will be used in the delay block
delay_a = A_T*tau_a;
delay_r = A_T*tau_r;

%% temp. dependent concentration of LacI
TTo = (T/To)^b;
Cr_T = (Cr_max - Cr_min)*(TTo/(1 + TTo)) + Cr_min;

%history values
history_r = 1300; 
history_a = 1300;
Parameter_Vec = [alpha_r; alpha_a; gamma_r; gamma_a; f; Ca; Cr_T; N; beta; Ro; A_T];

%% collecting timeperiod data as function of alpha_a by solving the dde multiple times using the dde3 function
% change num_iter to take a range of values of alpha_a, not that the
% optimal value of alpha a as taken by the author of temperature
% Compensation paper  will be kept as the middle value. For custom ranges
% of Alpha , change range alphaA to for example 90 : 1: 300 for taking
% alpha_A from 90 to 300
num_iter = 10;
avgperiod_laci = zeros(1,num_iter);
avgperiod_arac = zeros(1,num_iter);
range_alphaA = alpha_a-(num_iter/2 - 1): 1 : alpha_a+(num_iter/2);
for j = 1:num_iter
    Parameter_Vec = [alpha_r; range_alphaA(j); gamma_r; gamma_a; f; Ca; Cr_T; N; beta; Ro; A_T];
    sol = dde23(@(t,y,Z)dynamic_LacI_AraC(t,y,Z,Parameter_Vec),[delay_a, delay_r],[history_r, history_a],[0, 50]);
    [peaks, locations] = findpeaks(sol.y(1,:), sol.x); % Peaks of LacI
    [peaks2, locations2] = findpeaks(sol.y(2,:), sol.x); % Peaks of AraC
    [m,n] = size(locations);
    timeperiod_lacI = zeros(1,n);
    for i = 2:1:n
        timeperiod_lacI(1,i) = locations(1,i) - locations(1,i-1);
    end
avgperiod_laci(1,j) = sum(timeperiod_lacI)/(n-1);

    [m2,n2] = size(locations2);
    timeperiod_arac = zeros(1,n2);
    for i = 2:1:n2
        timeperiod_arac(1,i) = locations2(1,i) - locations2(1,i-1);
    end
avgperiod_arac(1,j) = sum(timeperiod_arac)/(n2-1);
end 
    
plot(range_alphaA,avgperiod_laci,'LineWidth', 2,'DisplayName', 'TimePeriod lacI');
hold on
plot(range_alphaA,avgperiod_arac,'LineWidth', 2,'DisplayName', 'TimePeriod AraC');
hold off
legend;
xlabel('Different values of Alpha a');
ylabel('Time period (min) ');
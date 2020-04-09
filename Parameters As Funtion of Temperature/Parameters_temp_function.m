

%% for multiple temperatures
num_iter = 11;
avgperiod_laci = zeros(1,num_iter);
avgperiod_arac = zeros(1,num_iter);
Tarray = 36:0.5:41;
y = 1;
for j = 36:0.5:41
    %% Parameter Values
    %maximal transcription rates
    T = j;
    multfactor = (2.5/10)*(j-36);
    %inc_duetotemp = 2.5 * (j-)/()
    alpha_r =  265 * multfactor;
    alpha_a = 92.75 * multfactor ;
    %transcriptional delay times
    tau_r = 13.5;
    tau_a = 15;
    %dilution rate constant
    beta = 0.0275 * multfactor;
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

    Parameter_Vec = [alpha_r; alpha_a; gamma_r; gamma_a; f; Ca; Cr_T; N; beta; Ro; A_T];
    %% dde23 solver
    % constant history of a(t) and r(t)
    history_r = 300; 
    history_a = 300;
    Tspan = [0,5]; %[initial time, final time]
    sol = dde23(@(t,y,Z)dynamic_LacI_AraC(t,y,Z,Parameter_Vec),[delay_a, delay_r],[history_r, history_a],[0, 100]);
    
    %% find time period of oscillations
    [peaks, locations] = findpeaks(sol.y(1,:), sol.x);
    [m,n] = size(locations);
    timeperiod_lacI = zeros(1,n);
    for i = 2:1:n
        timeperiod_lacI(1,i) = locations(1,i) - locations(1,i-1);
    end
    avgperiod_laci(y) = sum(timeperiod_lacI)/(n-1);
    
    [peaks2, locations2] = findpeaks(sol.y(2,:), sol.x);
    [m2,n2] = size(locations2);
    timeperiod_arac = zeros(1,n2);
    for i = 2:1:n2
        timeperiod_arac(1,i) = locations2(1,i) - locations2(1,i-1);
    end
    avgperiod_arac(y) = sum(timeperiod_arac)/(n2-1);
    
    y = y+1;

end



 plot(Tarray,avgperiod_laci,'LineWidth', 2,'DisplayName', 'TimePeriod LacI');
 hold on
 plot(Tarray,avgperiod_arac,'LineWidth', 2,'DisplayName', 'TimePeriod AraC');
 hold off
 legend;
 xlabel('Temperature T in Degrees');
 ylabel('Time Period of LacI and AraC');
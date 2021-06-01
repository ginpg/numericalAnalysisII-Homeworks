pkg load optim
clearvars;close all;clc;

%Giselt Parra, 26609640

function [A] = getA(alpha, gamma, pipi, kappa, eta, delta)
    tgamma1 = gamma(1) - kappa;
    A = zeros(8);                                     % 1: N (No susceptibles)
    A(1,2) = eta;                                     % 2: S (Susceptibles)
    A(2,2) = -eta; A(2,7) = delta;                    % 3: L (Latentes)
    A(3,3) = -alpha(1);                               % 4: P (Presintomáticos)
    A(4,3:4) = [alpha(1)*pipi, -alpha(2)];            % 5: T (Sintomáticos) 
    A(5,4:5) = [alpha(2), -gamma(1)];                 % 6: A (Asintomáticos)
    A(6,3:6) = [alpha(1)*(1-pipi), 0, 0, -gamma(2)];  % 7: R (Recuperados)
    A(7,5:7) = [tgamma1, gamma(2), -delta];           % 8: D (Fallecidos)
    A(8,5) = kappa;

end

function [Y] = RK4(Fun,Y,A,F,dt)                      % Runge-Kutta of order 4
    k_1 = Fun(Y,A,F);
    k_2 = Fun(Y+0.5*dt*k_1,A,F);
    k_3 = Fun(Y+0.5*dt*k_2,A,F);
    k_4 = Fun(Y+k_3*dt,A,F);
    Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
end
    
    
function [NS, S, L, P, T, A, R, D, VITOT]  = COVID19E(p,t,h,y)
  
    
    alpha = [p(1) p(2)];
    beta = [p(3)  p(4)  p(5)];
    gamma = [p(6)  p(7)];
    pipi = p(8);
    kappa = p(9);
    eta = p(10);
    delta = p(11);

    Npop  = p(12);
    NS0   = p(13);                      % Número inicial de no susceptibles
    L0    = p(14);                      % Número inicial de latentes
    P0    = p(15);                      % Número inicial de presintomáticos 
    T0    = p(16);                      % Número inicial de sintomáticos
    A0    = p(17);                      % Número inicial de asintomáticos
    R0    = p(18);                      % Número inicial de recuperados
    D0    = p(19);
    S0    = p(20);
    
    N = numel(t);
    Y = zeros(8,N);


    ITOT = L0+P0+T0+A0+R0+D0;                % Cantidad de infecciones que se acumulan durante el período
    VITOT = [ITOT];                          % Vector con el histórico de infecciones acumuladas
    Y(1,1) = NS0;                            % No Susceptibles
    Y(2,1) = Npop-NS0-L0-P0-T0-A0-R0-D0;     % Suceptibles iniciales
    Y(3,1) = L0;                             % Latentes
    Y(4,1) = P0;                             % Pre-sintomáticos (Infecciosos)
    Y(5,1) = T0;                             % Sintomáticos
    Y(6,1) = A0;                             % Asintomáticos
    Y(7,1) = R0;                             % Recuperados
    Y(8,1) = D0;                             % Fallecidos


    if round(sum(Y(:,1))-Npop)~=0
        error('La suma debe ser cero porque la población (incluyendo los fallecidos) se asume constante');
    end


    betaP = beta(1);    % Tasa de infección presintomáticos
    betaT = beta(2);    % Tasa de infección sintomáticos
    betaA = beta(3);    % Tasa de infección asintomáticos


    modelFun = @(Y,A,F) A*Y + F;
    dt = median(diff(t));


    for ii=1:N-1
       % N-1, S-2, L-3, P-4, T-5, A-6, R-7, D-8
       A = getA(alpha, gamma, pipi, kappa, eta, delta);     % Matriz A (Parte lineal)
       FI = betaP*Y(4,ii) + betaT*Y(5,ii) + betaA*Y(6,ii);  % Fuerza de infección
       SF = FI*Y(2,ii);                                     % Fuerza de infección sintomática aplicada
       F = zeros(8,1);                                      % Vector parte no lineal
       F(2:3,1) = SF*[-1/Npop; 1/Npop];
       Nuevos = F(3,1);
       ITOT = ITOT + F(3,1);
       VITOT = [VITOT, ITOT];
       Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);             % Resolución ODE
    end


    NS = Y(1,1:N);
    S = Y(2,1:N);
    L = Y(3,1:N);
    P = Y(4,1:N);
    T = Y(5,1:N);
    A = Y(6,1:N);
    R = Y(7,1:N);
    D = Y(8,1:N);
    


end

function VITOT  = MODEL(p,t,h,y)
    
    alpha = [p(1) p(2)];
    beta = [p(3)  p(4)  p(5)];
    gamma = [p(6)  p(7)];
    pipi = p(8);
    kappa = p(9);
    eta = p(10);
    delta = p(11);

    Npop  = p(12);
    NS0   = p(13);                        % Número inicial de no susceptibles
    L0    = p(14);                      % Número inicial de latentes
    P0    = p(15);                      % Número inicial de presintomáticos 
    T0    = p(16);                      % Número inicial de sintomáticos
    A0    = p(17);                      % Número inicial de asintomáticos
    R0    = p(18);                      % Número inicial de recuperados
    D0    = p(19);
    S0    = p(20);
    
    N = numel(t);
    Y = zeros(8,N);


    ITOT   = L0+P0+T0+A0+R0+D0;                % Cantidad de infecciones que se acumulan durante el período
    VITOT  = [ITOT];                          % Vector con el histórico de infecciones acumuladas
    Y(1,1) = NS0;                            % No Susceptibles
    Y(2,1) = Npop-NS0-L0-P0-T0-A0-R0-D0;     % Suceptibles iniciales
    Y(3,1) = L0;                             % Latentes
    Y(4,1) = P0;                             % Pre-sintomáticos (Infecciosos)
    Y(5,1) = T0;                             % Sintomáticos
    Y(6,1) = A0;                             % Asintomáticos
    Y(7,1) = R0;                             % Recuperados
    Y(8,1) = D0;                             % Fallecidos


    if round(sum(Y(:,1))-Npop)~=0
        error('La suma debe ser cero porque la población (incluyendo los fallecidos) se asume constante');
    end


    betaP = beta(1);    % Tasa de infección presintomáticos
    betaT = beta(2);    % Tasa de infección sintomáticos
    betaA = beta(3);    % Tasa de infección asintomáticos


    modelFun = @(Y,A,F) A*Y + F;
    dt = median(diff(t));


    for ii=1:N-1
       % N-1, S-2, L-3, P-4, T-5, A-6, R-7, D-8
       A = getA(alpha, gamma, pipi, kappa, eta, delta);     % Matriz A (Parte lineal)
       FI = betaP*Y(4,ii) + betaT*Y(5,ii) + betaA*Y(6,ii);  % Fuerza de infección
       SF = FI*Y(2,ii);                                     % Fuerza de infección sintomática aplicada
       F = zeros(8,1);                                      % Vector parte no lineal
       F(2:3,1) = SF*[-1/Npop; 1/Npop];
       Nuevos = F(3,1);
       ITOT = ITOT + F(3,1);
       VITOT = [VITOT, ITOT];
       Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);             % Resolución ODE
    end


    NS = Y(1,1:N);
    S = Y(2,1:N);
    L = Y(3,1:N);
    P = Y(4,1:N);
    T = Y(5,1:N);
    A = Y(6,1:N);
    R = Y(7,1:N);
    D = Y(8,1:N);
    

    set(h,'ydata',VITOT);
    legend('Nro casos acumulados en el historico', 'Nro de casos acumulados del modelo');
    drawnow;
    pause(.04);


end

%________________________________________________________________________________________________________________________________________________


ND = 156;                                     % Número de días de simulación
dt = 1;
time1=[1:ND];
N = numel(time1);
t = [0:N-1].*dt;
Npop =  30000000;
datay = [2    2    17    17    33    36    42    65    70    70    84    84    106    107    113    119    119    135    135    144    146    153    155    159    165    166    167    171    175    175    181    189    189    197    204    204    227    256    285    288    298    298    318    323    325    329    329    331    333    335    345    357    357    361    379    381    388    402    414    422    423    440    455    459    504    541    618    749    824    882    944    1010    1121    1177    1211    1245    1327    1370    1459    1510    1662    1819    1952    2087    2145    2316    2377    2473    2632    2738    2814    2879    2904    2978    3062    3150    3386    3483    3591    3789    3917    4048    4187    4366    4563    4779    5130    5297    5530    5832    6062    6273    6537    6750    7169    7411    7693    8008    8803    9178    9465    9707    10010    10428    10854    11191    11483    11891    12334    12774    13164    13613    14263    14929    15463    15988    16571    17158    17859    18574    19443    20206    20754    21438    22299    23280    24166    24961    25805    26800    27938    29088    30369    31381    32607    33755];
datax = 0:155; 


figure();
plot(datax,datay,'ko');
hold on;
h = plot(datax,datax,'b-','LineWidth',1.25);
hold off;

NS0 = 0;
S0 = Npop;
L0 =  29;
P0 =  1;
T0 =  1;
A0 =  1;
R0 =  1;
D0 =  1;

NS0l = 0;                     % Número inicial de no susceptibles
L0l = 1;                      % Número inicial de latentes
P0l = 1;                      % Número inicial de presintomáticos 
T0l = 1;                      % Número inicial de sintomáticos
A0l = 1;                      % Número inicial de asintomáticos
R0l = 1;                      % Número inicial de recuperados
D0l = 1;                      % Número inicial de fallecidos

NS0u = 0;                       % Número inicial de no susceptibles
L0u = 196;                      % Número inicial de latentes
P0u = 122;                      % Número inicial de presintomáticos 
T0u = 120;                      % Número inicial de sintomáticos
A0u = 165;                      % Número inicial de asintomáticos
R0u = 10;                        % Número inicial de recuperados
D0u = 10;                        % Número inicial de fallecidos


alphal   =  [0.0001   0.0001];                 
betal    =  [0.0001   0.0001   0.0001];        
gammal   =  [0.0001   0.0001];                  
pipil    =  0.0001;                            
kappal   =  0.0001;                           
etal     =  0.0001;                            
deltal   =  0;

alphau   =  [0.9   0.9];                 
betau    =  [0.9   0.9   0.9];        
gammau   =  [0.9   0.9];                  
pipiu    =  0.9;                            
kappau   =  0.9;                           
etau     =  0.9;                   
deltau   =  0;


  alpha   =   [0.03462452   0.03462452];
  beta    =   [0.24805   0.76061   0.13366];
  gamma   =   [3.9015e-8   1.1685e-01];              
  pipi    =   0.25221;
  kappa   =   0.0135;
  eta     =   0.002253;
  delta   =   0;



xi = [alpha beta gamma pipi kappa eta delta Npop NS0 L0 P0 T0 A0 R0 D0 S0];
lb = [alphal betal gammal pipil kappal etal deltal Npop NS0l L0l P0l T0l A0l R0l D0l S0];
ub = [alphau betau gammau pipiu kappau etau deltau Npop NS0u L0u P0u T0u A0u R0u D0u S0];

[p,resnorm,residual,exitflag,output] = lsqcurvefit(@(xi) MODEL(xi,t,h,datay), xi, time1, datay,lb,ub);
resnorm

%Parametros estimados de salida
alpha = [p(1) p(2)]
beta = [p(3)  p(4)  p(5)]
gamma = [p(6)  p(7)]
pipi = p(8)
kappa = p(9)
eta = p(10)
delta = p(11)

Npop  = p(12)
NS0   = p(13)
L0    = p(14)
P0    = p(15)
T0    = p(16)
A0    = p(17)
R0    = p(18)
D0    = p(19)
xi = [alpha beta gamma pipi kappa eta delta Npop NS0 L0 P0 T0 A0 R0 D0 S0];

figure();
plot(datax,datay,'ko',datax,MODEL(xi,t,h,datay))';
legend('Nro casos acumulados en el historico', 'Nro de casos acumulados del modelo');
xlabel('Número de días');
ylabel('Nro casos acumulados');
title('Estimacion del historico: casos acumulados');
saveas(gcf,'Graficos1/Funcion.png');
MODEL(xi,t,h,datay)'

[NS, S, L, P, T, A, R, D, VITOT] = COVID19E(xi,t,h,datay);


PerNS = (NS/Npop)*1000000;
PerS = (S/Npop)*1000000;
PerL = (L/Npop)*1000000;
PerP = (P/Npop)*1000000;
PerT = (T/Npop)*1000000;
PerA = (A/Npop)*1000000;
PerR = (R/Npop)*1000000;
PerD = (D/Npop)*1000000;
PerAcum = (VITOT/Npop)*1000000;

PerINAis = PerP + PerT + PerA;  % Proporción de Infeciosos No Aislados
PerIN = PerP + PerT + PerA;     %Proporción de Infeciosos
PerLate = PerL;                 %Proporción de Latentes
PerAI = (PerA./PerINAis);       %Proporción de Asintomáticos sobre Infecciosos No Aislados
PerTI = (PerT./PerINAis);       %Proporción de Sintomáticos sobre Infecciosos No Aislados
PerPI = (PerP./PerINAis);       %Proporción de Presintomáticos sobre Infecciosos No Aislados

figure
clf; close all;

figure
clf; close all;


plot(time1, PerIN, 'r', time1, PerLate, 'b','linewidth',1.5);
legend('Infecciosos', 'Latentes');
xlabel('Número de días');
ylabel('Proporción de la población (por millón)');
title('Proporción de la población Latente e Infecciosa');
saveas(gcf,'Graficos1/ProporcionLatentesInfecciosos.png')



% Despliegue de gráfico de estados Infecciosos
plot(time1, PerT, 'r', time1, PerA, 'b', time1, PerP, 'm','linewidth',1.5);
legend('Sintomáticos', 'Asintomáticos', 'Pre-sintomáticos');
xlabel('Número de días');
ylabel('Proporción de la población (por millón)');
title('Proporción de la población en Estados Infecciosos');
saveas(gcf,'Graficos1/ProporcionInfecciosos.png');


% Despliegue de gráfico de estados Infecciosos Proporción
plot(time1, PerTI*100, 'r', time1, PerAI*100, 'b', time1, PerPI*100, 'm','linewidth',1.5);
legend('Sintomáticos', 'Asintomáticos', 'Pre-sintomáticos');
xlabel('Número de días');
ylabel('Porcentaje');
title('Porcentaje en estados infecciosos');
saveas(gcf,'Graficos1/PorcentajeRespectoInfecciosos.png')



% Despliegue de gráfico de estados Recuperados
plot(time1, PerR, 'g', 'linewidth',1.5);
legend('Recuperados');
xlabel('Número de días');
ylabel('Proporción de la población (por millón)');
title('Proporción de la población Recuperada');
saveas(gcf,'Graficos1/ProporcionRecuperados.png')



% Despliegue de gráfico de estados Fallecidos
plot(time1, PerD, 'k', 'linewidth',1.5);
legend('Fallecidos');
xlabel('Número de días');
ylabel('Proporción de la población (por millón)');
title('Proporción de la población Fallecida');
saveas(gcf,'Graficos1/ProporcionFallecidos.png')



% Despliegue de gráfico de estados Susceptibles y No Susceptibles
plot(time1, PerS, 'm',time1, PerNS, 'k', 'linewidth',1.5);
legend('Susceptibles', 'No Susceptibles');
xlabel('Número de días');
ylabel('Proporción de la población (por millón)');
title('Proporción de la población en Estados Susceptibles y No Susceptibles ');
saveas(gcf,'Graficos1/ProporcionSusceptiblesNoSusceptibles.png')

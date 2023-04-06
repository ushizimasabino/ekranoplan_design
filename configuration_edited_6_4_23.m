clc;
clear all;
close all;

% Cargo = 5400;                       %ft^3
Mcruise = 0.4;                        % cruise Mach number
% Range = 7500;                       % Nautical Miles
% RangeFt = Range * 6076.12;          % nmi to ft
% Alt = 35000;                        % ft
% Altm = Alt*0.3048;                  % m
Vapproach = 100;                    % knots
Vapproachft = Vapproach*1.68781;    % ft/s
%CargoDensity = 10; % lbs/ft^3

TEU_volume = 1172;                  % ft^3
TEU_gross_mass = 56000;             % Loaded mass per TEU
TEU_mass = 4780;                    % lbs
TEU_n = 60;                        % Number of TEUs
Range = 6000;                       % Nautical Miles - LA - Shanghai
RangeFt = Range * 6076.12;          % nmi to ft
Alt = 100;                          % ft ***

W_TEU = TEU_gross_mass*TEU_n;
W_TEU_empty = TEU_mass*TEU_n;

Ncrew = 6;                         % Number of crew
PassWeight = 220;                   % Weight of human in lbs
Wcrew = Ncrew * PassWeight;
Wpayload = Wcrew + W_TEU;

%% Design Gross Takeoff Weight Estimate

L_Dfactor = 1;                  % Propeller most efficient cruise L/D factor
% L_Dfactor = 0.866;              % Turbojet most efficient cruise L/D factor
L_D = 35;                           % Estimate from Raymer
a35000 = 969.16;                    % Speed of sound at 35000ft
Vcruise = Mcruise * a35000; 
SFCcruise = 0.3./(3600);              % Estimate from Raymer from 1/hr to 1/s

W1_0 = 0.95;                        % Takeoff Weight factor *** better estimate?
W3_2 = exp((-RangeFt*SFCcruise)./(Vcruise*L_D*L_Dfactor)); % Breguet range equation
W5_4 = 0.995;                       % Landing Weight factor

Ws_o = W1_0*W3_2*W5_4;    % Whole mission Weight fraction
Wf_o = 1.06.*(1-Ws_o);              % Fuel Weight fraction

C = -0.05; % Flying boat estimation
A = 1.09;
K = 1; % fixed sweep

f = @(x) (Wpayload)./(1-Wf_o-x);  %Eqn 3.4
g = @(x) exp(1./C.*log(1./(A.*K).*x));  %Table 3.1
h = @(x) f(x) - g(x);                   %Function handle for intersection
x_intersect = fzero(h, 0.4);
y_intersect = f(x_intersect);

% range of empty weight fraction
we1 = linspace(0.4,0.5,1000); % we/w0 Box 3.1 0.436 to 0.432

plot(we1,f(we1),'b-',we1,g(we1),'r-',x_intersect,y_intersect,'go','LineWidth',2);
legend('Eqn3.4','Table3.1')
xlabel('weight fraction: we/w0')
ylabel('takeoff weight w0')

TEU_frac_gross = W_TEU_empty./y_intersect;
TEU_frac_empty = W_TEU_empty./(y_intersect.*x_intersect);

%% Thrust to Weight Ratio and Wing Loading 

SweepQuart = 0;                % fig 4.19
ClMax = 1.7;                    % Double slotted flaps (same as 767) fig 5.3
Cdo = 0.0105*0.777;             % For jet aircraft section 5.3.7
AR = 3;                         % Aspect ratio from chapter 3
TOP4 = 300;                     % Takeoff parameter for 4 engines Fig 5.4
e = 0.85;                        % Oswald efficieny factor section 5.3.7
gft = 32;                       % ft/s^2

Vstall = Vapproachft./1.3;

rhosl = 1.225; 
rho_cruise = 1.225;

% [T,a,P,rho_cruise] = atmosisa(Altm);        % Atmosphere data
rho_cruise_ft = rho_cruise./515.37881852553;% Convert to slugs/ft^3

rho_td = 1.225;                             % rho touchdown kg/m^3
rho_td_ft = rho_td./515.37881852553;        % Convert to slugs/ft^3

sigma = rho_td./rhosl;                      % Density ratio

qcruise = 0.5*rho_cruise_ft*Vcruise.^2; % Dynamic pressure at cruise 

vto = 249.333; % ft/s takeoff speed of 737
dhdt = 10; % ft/s climb rate set by FAR 25


f = @(x) x./(TOP4*sigma*ClMax);       % Eqn 5.9, f = W/S, x = T/W, Takeoff
g = @(x) ((qcruise*Cdo)./x) + x.*(1./(qcruise.*pi.*AR.*e)); % Eq 5.24, g = W/S, Cruise
W_Sstall = Vstall.^2*rho_td_ft*ClMax*0.5; % Relation between W/S and Vstall, Stall

W_S = linspace(0,100,1000);         % W/S from 0 to 200

figure()
plot(W_S,f(W_S),'LineWidth',2)
hold on
plot(W_S,g(W_S),'LineWidth',2)
hold on
xline(W_Sstall,'LineWidth',2)
legend('Takeoff','Cruise','VStall')
xlabel('W/S - Wing Loading (lb/ft^2)')
ylabel('T/W - Thrust to Weight ratio')

ylim([0,0.4])
xlim([0,100])

W_Sgraph = 35; % Selected W_S and T_W from graph
T_Wgraph = 0.08;

Sref = y_intersect/W_Sgraph; % Reference area from graph ft^2 

b = sqrt(Sref*AR); % Wing span 

%% Plane sizing

a = 0.67; % Table 6.3 for general aviation - twin engine 
C = 0.43;

fuselage_length_raymer = a*y_intersect.^C; % ft


Thrust_req = T_Wgraph .* y_intersect; % Total thrust required in lbf
Thrust_per_engine = Thrust_req/2;
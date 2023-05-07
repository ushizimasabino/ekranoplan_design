clear
close all
clc
%% Constants
phi = 0.25;
LHV = 43E6; % kJ/kg
rF = 1.5;
rP = linspace(20,70);
BPR = linspace(5,15);
effFC = 0.9;
effT = 0.95;
gamma = 1.4;
R = 287; % J/kg*K
M = 0.7;
T0 = 288.15; % K
T1 = T0; % K
p0 = 101.325; % kPa
p1 = p0; % K
govergless1 = gamma/(gamma-1);
gless1over2 = (gamma-1)/2;
gless1overg = 1/govergless1;
cp = govergless1*R; % J/kg*K
p03 = zeros(100,100);
T03 = zeros(100,100);
p04 = zeros(100,100);
T04 = zeros(100,100);
T045 = zeros(100,100);
p045 = zeros(100,100);
T05 = zeros(100,100);
p05 = zeros(100,100);
T09 = zeros(100,100);
p09 = zeros(100,100);
Vc = zeros(100,100);
specificnetthrust = zeros(100,100);
tsfc = zeros(100,100);
effProp = zeros(100,100);
%% Constant Calculations
% Diffuser
p01 = p1*(1+gless1over2*M^2)^govergless1; % kPa
p02 = p01; % kPa
T01 = T1*(1+gless1over2*M^2); % K
T02 = T01; % K
% Fan
p013 = rF*p02; % kPa
p023 = p013; % kPa
T013 = T02*(1+(1/effFC)*(rF^gless1overg-1)); % K
T023 = T013; % K
% Bypass
T019 = T013; % K
p019 = p013; % kPa
% Combustor
mfoverma = (12*12.01+26*1.008)/(18.5*2*16+18.5*3.76*2*14.01);
% Velocity
Vb = sqrt(2*govergless1*R.*T013.*(1-(p1./p013).^gless1overg)); % m/s
V = M*sqrt(gamma*R*T0); % m/s
%% Variable Calculations
for i = 1:100
    for j = 1:100
        % Core
        p03(i,j) = rP(i)*p023; % kPa
        T03(i,j) = T023*(1+(1/effFC)*(rP(i).^gless1overg-1)); % K
        p04(i,j) = p03(i,j); % kPa
        T04(i,j) = T03(i,j)+mfoverma*LHV*phi/cp; % K
        T045(i,j) = T04(i,j)-(T023/effFC)*(rP(i).^gless1overg-1); % K
        p045(i,j) = p04(i,j).*(1-((T023)./(effFC*effT*T04(i,j))).*(rP(i).^gless1overg-1)).^govergless1; % kPa
        T05(i,j) = T045(i,j)-(1+BPR(j)).*(T02/effFC).*(rF^gless1overg-1); % K
        p05(i,j) = p045(i,j).*(1-(1+BPR(j)).*((T02)./(effFC*effT*T045(i,j))).*(rF^gless1overg-1)).^govergless1; % kPa
        p09(i,j) = p05(i,j); % kPa
        T09(i,j) = T05(i,j); % K
        % Velocity
        Vc(i,j) = sqrt(2*govergless1*R.*T05(i,j).*(1-(p1./p05(i,j)).^gless1overg)); % m/s
        specificnetthrust(i,j) = Vc(i,j)-V+BPR(j).*(Vb-V);
        tsfc(i,j) = phi*mfoverma./specificnetthrust(i,j);
        effProp(i,j) = (2*specificnetthrust(i,j)*V)/(BPR(j).*(Vb^2-V^2)+(Vc(i,j).^2-V^2));
        if (1-(p1/p05(i,j))^gless1overg) < 0
             specificnetthrust(i,j) = NaN;
             tsfc(i,j) = NaN;
             effProp(i,j) = NaN;
        end
    end
end
%% Plots
% Specific Net Thrust
figure()
ThrustGraph = pcolor(BPR,rP,specificnetthrust);
title('Specific Net Thrust')
xlabel('BPR')
ylabel('Pressure Ratio')
ThrustGraph.FaceColor = 'interp';
% Thrust Specific Fuel Consumption
figure()
tsfcGraph = pcolor(BPR,rP,real(tsfc));
title('Thrust Specific Fuel Consumption')
xlabel('BPR')
ylabel('Pressure Ratio')
tsfcGraph.FaceColor = 'interp';
% Propulsive Efficiency
figure()
propGraph = pcolor(BPR,rP,real(effProp));
title('Propulsive Efficiency')
xlabel('BPR')
ylabel('Pressure Ratio')
propGraph.FaceColor = 'interp';
%% Diameter calculations
maxspecificnetthrustcolumns = max(specificnetthrust);
maxspecificnetthrust = max(maxspecificnetthrustcolumns)
maxpropeffcolumns = max(effProp);
maxpropeff = max(maxpropeffcolumns)
mintsfccolumns = min(tsfc);
mintsfc = min(mintsfccolumns)
Vnew = 293.333;
masscore = 118.33;
rho = 1.225;
coreD = sqrt(masscore*4/(rho*Vnew*pi()));
massbypass = masscore*13.28;
bypassD = sqrt((4*massbypass/(rho*pi()*Vnew))+coreD^2);
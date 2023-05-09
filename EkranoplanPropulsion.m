clear
close all
clc
%% Constants
% Design Parameters
LHV = 43E6; % kJ/kg
rF = 1.3;
rP = linspace(20,60);
BPR = linspace(6,16);
cooling = 0.15; % mass of cooling stream over total mass flow
% Design Limitations
nozzlepressureratio = 0.985;
diffuserpressureratio = 0.97;
combustpressureratio = 0.96;
effCombust = 0.995;
effF = 0.89;
effC = 0.9;
effT = 0.89;
effN = 1;
effD = 1;
TIT = 2000; % K
% Perfect Gas Assumption
gamma = 1.4;
R = 287; % J/kg*K
M = 0.7;
% Ambient Values
T0 = 288.15; % K
p0 = 101.325; % kPa
T1 = T0; % K
p1 = p0; % K
govergless1 = gamma/(gamma-1);
gless1over2 = (gamma-1)/2;
gless1overg = 1/govergless1;
cp = govergless1*R; % J/kg*K
phi = zeros(100,100);
p0_3 = zeros(100,100);
T0_3 = zeros(100,100);
p0_4 = zeros(100,100);
T0_4 = zeros(100,100);
T0_45 = zeros(100,100);
p0_45 = zeros(100,100);
T0_5 = zeros(100,100);
p0_5 = zeros(100,100);
T0_9 = zeros(100,100);
p0_9 = zeros(100,100);
Vc = zeros(100,100);
specificnetthrust = zeros(100,100);
tsfc = zeros(100,100);
effProp = zeros(100,100);
%% Constant Calculations
% Diffuser
p0_1 = p1*(1+gless1over2*M^2)^govergless1; % kPa
p0_2 = p0_1*diffuserpressureratio; % kPa
T0_1 = T1*(1+gless1over2*M^2); % K
T0_2 = T0_1; % K
% Fan
p0_13 = rF*p0_2; % kPa
p0_23 = p0_13; % kPa
T0_13 = T0_2*(1+(1/effF)*(rF^gless1overg-1)); % K
T0_23 = T0_13; % K
% Bypass
T0_19 = T0_13; % K
p0_19 = p0_13; % kPa
% Combustor
mfoverma = (12*12.01+26*1.008)/(18.5*2*16+18.5*3.76*2*14.01); % stoichiometric fuel to air mass ratio
% Velocity
Vb = sqrt(2*govergless1*R.*T0_13.*(1-(p1./p0_13).^gless1overg)); % m/s
V = M*sqrt(gamma*R*T0); % m/s
%% Variable Calculations
for i = 1:100
    for j = 1:100
        % Core
        p0_3(i,j) = rP(i)*p0_23; % kPa
        T0_3(i,j) = T0_23*(1+(1/effC)*(rP(i).^gless1overg-1)); % K
        p0_4(i,j) = p0_3(i,j).*combustpressureratio; % kPa
        phi(i,j) = (cp.*(TIT-T0_3(i,j)))./(mfoverma*LHV); % K
        T0_4(i,j) = T0_3(i,j)+(1-cooling)*effCombust*mfoverma*LHV*phi(i,j)/cp; % K
        T0_45(i,j) = T0_4(i,j)-(T0_23/effC)*(rP(i).^gless1overg-1); % K
        p0_45(i,j) = p0_4(i,j).*(1-((T0_23)./(effC*effT*T0_4(i,j))).*(rP(i).^gless1overg-1)).^govergless1; % kPa
        T0_5(i,j) = T0_45(i,j)-(1+BPR(j)).*(T0_2/effF).*(rF^gless1overg-1); % K
        p0_5(i,j) = p0_45(i,j).*(1-(1+BPR(j)).*((T0_2)./(effF*effT*T0_45(i,j))).*(rF^gless1overg-1)).^govergless1; % kPa
        p0_9(i,j) = p0_5(i,j); % kPa
        T0_9(i,j) = T0_5(i,j); % K
        % Velocity
        Vc(i,j) = sqrt(2*govergless1*R.*T0_5(i,j).*(1-(p1./p0_5(i,j)).^gless1overg)); % m/s
        specificnetthrust(i,j) = Vc(i,j)-V+BPR(j).*(Vb-V);
        tsfc(i,j) = (1-cooling)*phi(i,j)*mfoverma./specificnetthrust(i,j);
        effProp(i,j) = (2*specificnetthrust(i,j)*V)/(BPR(j).*(Vb^2-V^2)+(Vc(i,j).^2-V^2));
        if (1-(p1/p0_5(i,j))^gless1overg) < 0
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
colorbar
% Thrust Specific Fuel Consumption
figure()
tsfcGraph = pcolor(BPR,rP,real(tsfc));
title('Thrust Specific Fuel Consumption')
xlabel('BPR')
ylabel('Pressure Ratio')
tsfcGraph.FaceColor = 'interp';
colorbar
% Propulsive Efficiency
figure()
propGraph = pcolor(BPR,rP,real(effProp));
title('Propulsive Efficiency')
xlabel('BPR')
ylabel('Pressure Ratio')
propGraph.FaceColor = 'interp';
colorbar
%% Diameter calculations
maxspecificnetthrustcolumns = max(specificnetthrust);
maxspecificnetthrust = max(maxspecificnetthrustcolumns);
maxpropeffcolumns = max(effProp);
maxpropeff = max(maxpropeffcolumns);
mintsfccolumns = min(tsfc);
mintsfc = min(mintsfccolumns);
myspecificthrust = 1.1856*10^3;
totalThrust = 1000000*4.44822; % N
numberEngines = [10 8];
engineThrust = totalThrust./numberEngines; % N
masscore = engineThrust./myspecificthrust; % kg/s
rho = 1.225; % kg/m^3
coreD = sqrt(masscore.*4/(rho*V*pi())) % m
massbypass = masscore.*9.8182; % kg/s
bypassD = sqrt((4*massbypass./(rho*pi()*V))+coreD.^2) % m

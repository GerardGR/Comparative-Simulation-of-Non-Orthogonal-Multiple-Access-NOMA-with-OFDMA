%% Scenario 1: Two users located within the cell 
%% -- Rural, subUrban and Urban areas
%% -- Analysis for the following carrier frequencies:
%% -- User 1 (weak user, cell-centered)
%% -- User 2 (strong user, cell-edge)

clear all;

%% ------------------------------------------------------------------------
%% Definition of the Scenario: Initial Assumptions
%% ------------------------------------------------------------------------
%% Pathloss and System Parameters
Conf_.MBS_ant_height= 180;     %% m  antenna height
Conf_.UE_ant_height= 5;        %% m  antenna height

Conf_.NumUEs = 2;                %% Number of users in the cell
Conf_.CellSizes = [0.5 1.5 4 10];    %% Cell size, km
Conf_.distances = [0.1 0.4];        %% User distances, km
Conf_.Flag = 0;                  % 0 --> Small urban areas
                                 % 1 --> Metropolitan urban areas
                                 % 2 --> subUrban areas
                                 % 3 --> Rural areas
Conf_.c = 3e8;

Conf_.MBS_Pmax_dB =  43;         %% Max Power Transmitted by the BS, dBm
Conf_.UE_Pmax_dB =  23;          %% Max Power Transmitted by the UE, dBm
Conf_.Antenna_gain = 15;         %% Antenna Gain, dBi [12-15]

Conf_.BW = 20*10^6;              %% Bandwidth, Hz [i.e 20 MHz]
Conf_.Freq_Carr= 450;           %% Carrier frequency, MHz, [i.e. 450 MHz]
Conf_.No = -174;                 %% Thermal noise density, dBm/Hz 
Conf_.F_UTnoise = 7;             %% Noise Figure, dB

%%
%% Noise
%%
Conf_.Eff_Noise  = Conf_.No + 10*log10(Conf_.BW)+Conf_.F_UTnoise;

%%
%% Path Loss Computation: 
%% 
[Pathloss] = pathloss_OH(Conf_.distances, Conf_.Freq_Carr, Conf_.MBS_ant_height, Conf_.UE_ant_height, Conf_.Flag);

%% Gains Computation:

G1 = 10^(0.1*(Conf_.Antenna_gain - Pathloss(1)));
G2 = 10^(0.1*(Conf_.Antenna_gain - Pathloss(2)));

channel_1 = sqrt(G1/2)*(randn(1,1)+1i*randn(1,1));
channel_2 = sqrt(G2/2)*(randn(1,1)+1i*randn(1,1));
H1 = abs(channel_1)^2;
H2 = abs(channel_2)^2;

G1 = G1*H1;
G2 = G2*H2;

P = 10^(0.1*(Conf_.MBS_Pmax_dB-30));
N = 10^(0.1*(Conf_.Eff_Noise-30));

SNR = 10*log10(P*G1/N);

%% -------------------------------
%% ------ Energy Efficiency ------
%% -------------------------------
%%
run EnergyEff.m;  %%1-> Effective Transmit Power
                  %%2-> Circuit Power
                  %%3-> Power Consumption
                  
%%
%% OMA (OFDMA) Computation: Rates
%%
count = 1;
for alpha = 0:0.01:1 %bandwidth splitting factor

P1 = P/2;
P2 = P/2;
R1(count) = alpha*Conf_.BW*log2(1 + P1*G1/(alpha*N));
R2(count) = (1-alpha)*Conf_.BW*log2(1 + P2*G2/((1-alpha)*N));
count = count + 1;
end

figure(1)
hold on;
plot(R1,R2,'k');
grid on;

%%
%% NOMA Computation: Rates
%%
count = 1;
for gamma = 0:0.01:1 % power splitting factor

P1 = P*gamma;
P2 = P - P1;
R1(count) = Conf_.BW*log2(1 + P1*G1/N);
R2(count) = Conf_.BW*log2(1 + P2*G2/(P1*G2 + N));
count = count + 1;
end

hold on;
plot (R1,R2,'r');
grid on;
xlabel('Rate of user 1 (bps)');
ylabel('Rate of user 2 (bps)');
title('Comparison Rates: NOMA vs OMA');
grid on;
box on;
legend('OFDMA','NOMA')

%%
%% NOMA Computation: EE
%%
count = 1;
for p = 1:1:100 %W
P1 = p*0.1; %allocate less power to UE1
P2 = p - P1;
R1 = Conf_.BW*log2(1 + P1*G1/N);
R2 = Conf_.BW*log2(1 + P2*G2/(P1*G2 + N));
Sum_Rate_NOMA = R1 + R2;
SE(count) = Sum_Rate_NOMA/Conf_.BW; % bit/sec/Hz
EE(count) = Sum_Rate_NOMA/(EnergyEff_.CP + p); % bit/watt.sec
count = count + 1;
end

figure(2)
hold on;
plot(SE,EE,'k');
xlabel('SE (bit/sec/Hz)');
ylabel('EE (bit/joule)');
grid on;

%%
%% OMA (OFDMA) Computation: EE
%%
count = 1;
for p = 1:1:100 %Watt
P1 = p/2;
P2 = p/2;
R1 = (Conf_.BW/2)*log2(1 + P1*G1/(N/2));
R2 = (Conf_.BW/2)*log2(1 + P2*G2/(N/2));
Sum_Rate_OMA = R1 + R2;
SE(count) = Sum_Rate_OMA/Conf_.BW; % bit/sec/Hz
EE(count) = Sum_Rate_OMA/(EnergyEff_.CP + p); % bit/watt.sec
count = count + 1;
end

hold on;
plot(SE,EE,'g-');
xlabel('SE (bit/sec/Hz)');
ylabel('EE (bit/joule)');
grid on;
box on;
legend('NOMA','OFDMA')
title('Comparison EE/SE: NOMA vs OMA');
hold off;

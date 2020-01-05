clear all;
close all;

%% Initial Assumptions
Conf_.CellSize = 0.5;
Conf_.NumUEs = 10;
Conf_.Freq_Carr = 450;
Conf_.MBS_ant_height = 180;
Conf_.UE_ant_height = 5;
Conf_.BW = 20*10^6;
Conf_.MBS_Pmax_dB =  43;         
Conf_.UE_Pmax_dB =  23; 
Conf_.Antenna_gain = 15; 
Conf_.No = -174;                 
Conf_.F_UTnoise = 7;
Conf_.Flag = 0;           % 0 --> Small Urban Areas
                          % 1 --> Metropolitan Urban Areas
                          % 2 --> SubUrban Areas
                          % 3 --> Rural Areas

Conf_.add_flag = 1;       % 0 -> NEXT QUEUED USER
                          % 1 -> MAX RATE USER
sumNOMA = 0;
sumOMA = 0;
Mark = zeros(1,30);
Mark_OMA = zeros(1,30);
xmark = zeros(1,30);
%% Setting up  the environment
% Computation of users' distances wrt the center-cell and returns the
% users that might be potentially functioning as NOMA.

for(seq=1:30)
[Usr_.distances,Usr_.angles, l_noma] = Cell_User_deployment(Conf_.CellSize, Conf_.NumUEs);

Num_usr_OMA = Conf_.NumUEs-length(l_noma);
Num_usr_NOMA = length(l_noma);

Usr_.PL_OH = zeros(1,Conf_.NumUEs);
Usr_.Gain = zeros(1,Conf_.NumUEs);
Usr_.Power = zeros(1,Conf_.NumUEs);
Usr_.Rate = zeros(1,Conf_.NumUEs);

% Bandwidth need to be splitted among users. 
BW_split = Conf_.BW/Conf_.NumUEs;
Conf_.Eff_Noise  = Conf_.No + 10*log10(Conf_.BW)+Conf_.F_UTnoise;
Noise = 10^(0.1*(Conf_.Eff_Noise-30));
% Reference Power BS
Power_ref = 10^(0.1*(Conf_.MBS_Pmax_dB-30));

%% Computation of the general and common attributes of all users
% Pathloss and channel Gain 
l_oma = [];
l_noma_new = [];
for ii = 1:Conf_.NumUEs
    Usr_.PL_OH(ii) = pathloss_OH(Usr_.distances(ii), Conf_.Freq_Carr, Conf_.MBS_ant_height, Conf_.UE_ant_height, Conf_.Flag);
    G = 10^(0.1*(Conf_.Antenna_gain - Usr_.PL_OH(ii)));
    channel = sqrt(G/2)*(randn(1,1)+1i*randn(1,1));
    H = abs(channel)^2;
    Usr_.Gain(ii) = G*H;
    
    % Listing the number of NOMA and OMA users initially
    if (ismember(ii,l_noma)) 
       disp('NOMA Usr');
    else
       disp('OMA Usr');
       l_oma = [l_oma ii];
    end
end

%% Applying the decision making to the NOMA potential users
for ii = 1:2:Num_usr_NOMA
    if Usr_.Gain(l_noma(ii)) > Usr_.Gain(l_noma(ii+1))
        G1 = Usr_.Gain(l_noma(ii));
        G2 = Usr_.Gain(l_noma(ii+1));
        flag = true;
    else
        G2 = Usr_.Gain(l_noma(ii));
        G1 = Usr_.Gain(l_noma(ii+1));
        flag = false;
    end
    [p1, p2] = eqPower_Allocation_NOMA(G1, G2, Power_ref, Noise);
    disp('Power Allocation');
    disp([double(p1),double(p2)]);
   
    if isempty(p1) || isempty(p2)
        % Test Fails: reallocating the BW and updating the number of Users
        % (OMA and NOMA counters)
        test_F = sprintf('\nNOMA test failed. Users %d and %d need to be configured to OMA\n ',l_noma(ii), l_noma(ii+1));
        fprintf(test_F);
        l_oma = [l_oma l_noma(ii) l_noma(ii+1)];
        
        Num_usr_OMA = Num_usr_OMA + 2;
        Num_usr_NOMA = Num_usr_NOMA - 2;
        BW_split = Conf_.BW/(Num_usr_OMA + Num_usr_NOMA*0.5);
                
    else
        % Test is a success: allocating the P1 and P2 to the specific
        % paired NOMA users.
        test_P = sprintf('NOMA test passed. Users %d and %d configured to NOMA\n ',l_noma(ii), l_noma(ii+1));
        fprintf(test_P);
        l_noma_new = [l_noma_new l_noma(ii) l_noma(ii+1)];
        if flag == true
            Usr_.Power(l_noma(ii)) = double(p1);
            Usr_.Power(l_noma(ii+1)) = double(p2);
        else
            Usr_.Power(l_noma(ii)) = double(p2);
            Usr_.Power(l_noma(ii+1)) = double(p1);
        end
    end
end

%% User-Rate computation depending on the radio access technology: OMA & NOMA

for jj = 1:Num_usr_OMA
    Usr_.Power(l_oma(jj)) = Power_ref;
    Usr_.Rate(l_oma(jj)) = (BW_split/Conf_.NumUEs)*log2(1 + Usr_.Power(l_oma(jj))*Usr_.Gain(l_oma(jj))/(Noise/Conf_.NumUEs));
end

for jj = 1:2:Num_usr_NOMA
    Usr_.Rate(l_noma_new(jj)) = BW_split*log2(1 + Usr_.Power(l_noma_new(jj))*Usr_.Gain(l_noma_new(jj))/Noise);
    Usr_.Rate(l_noma_new(jj+1)) = BW_split*log2(1 + Usr_.Power(l_noma_new(jj+1))*Usr_.Gain(l_noma_new(jj+1))/(Usr_.Power(l_noma_new(jj))*Usr_.Gain(l_noma_new(jj+1)) + Noise));
    
    % Update the plot
    hold on
    [x,y] = pol2cart(Usr_.angles(l_noma_new(jj)),Usr_.distances(l_noma_new(jj)));
    scatter(x,y,'co');
    [x2,y2] = pol2cart(Usr_.angles(l_noma_new(jj+1)),Usr_.distances(l_noma_new(jj+1)));
    scatter(x2,y2,'co');
end
hold off

%% Additional User: 
Num_Additional = Conf_.NumUEs - (Num_usr_OMA + Num_usr_NOMA*0.5);
l_addusers = zeros(1,2*Num_Additional); 
% NEXT QUEUED USER
if (Conf_.add_flag == 0)
    n = 1;
    m = 0;
    for kk = 1:Num_Additional
        m = m + 1;
        if (ismember(m,l_noma)) 
            aux_k = find(l_noma_new == m);
            Usr_.Rate(Conf_.NumUEs + n) = Usr_.Rate(l_noma_new(aux_k));
            Usr_.Rate(Conf_.NumUEs + n + 1) = Usr_.Rate(l_noma_new(aux_k+1));
            l_addusers(n) = l_noma_new(aux_k);
            l_addusers(n + 1) = l_noma_new(aux_k + 1);
            l_noma(aux_k) = 0;
            l_noma(aux_k+1) = 0;
            n = n + 1;
        elseif(~ismember(m,l_noma) && ismember(m,l_noma_new))
            m = m + 1;
            if ismember(m,l_noma_new)
                Usr_.Rate(Conf_.NumUEs + n) = Usr_.Rate(m);
                Usr_.Rate(Conf_.NumUEs + n + 1) = Usr_.Rate(l_noma_new(m + 1));
                l_addusers(n) = m;
                l_addusers(n + 1) = l_noma_new(m + 1);
                n = n + 1;
            else
                Usr_.Rate(Conf_.NumUEs + n) = Usr_.Rate(m);
                l_addusers(n) = m;
            end
        else
            Usr_.Rate(Conf_.NumUEs + n) = Usr_.Rate(m);
            l_addusers(n) = m;
        end
        n = n + 1;
    end
else
% MAX RATE USER
    Max_Rates = Usr_.Rate;
    n = 1;
    for kk = 1:Num_Additional
        [Max_rate, index] = max(Max_Rates);
        Max_Rates(index) = 0;
        l_addusers(n) = index;
        if (ismember(index,l_noma_new))
            aux_i = find(l_noma_new == index);
            if (mod(aux_i,2) == 0)
                Usr_.Rate(Conf_.NumUEs + n) = Max_rate;
                Usr_.Rate(Conf_.NumUEs + n + 1) = Usr_.Rate(l_noma_new(aux_i-1));
                Max_Rates(l_noma_new(aux_i-1)) = 0;
                l_addusers(n+1) = l_noma_new(aux_i-1);
            else
                Usr_.Rate(Conf_.NumUEs + n) = Max_rate;
                Usr_.Rate(Conf_.NumUEs + n + 1) = Usr_.Rate(l_noma_new(aux_i+1));
                Max_Rates(l_noma_new(aux_i+1)) = 0;
                l_addusers(n+1) = l_noma_new(aux_i+1);
            end
            n = n + 1;
        else
            Usr_.Rate(Conf_.NumUEs + n) = Max_rate;
        end
        n = n + 1;
    end
end

% SumRate NOMA - OMA scenario
Usr_.SumRate = sum(Usr_.Rate);

%%
%% FULL - OMA 
%%

%% User-Rate computation 

for jj = 1:Conf_.NumUEs
    Usr_.Rate_OMA(jj) = (BW_split/Conf_.NumUEs)*log2(1 + Power_ref*Usr_.Gain(jj)/(Noise/Conf_.NumUEs));
end
% SumRate FULL - OMA scenario
Usr_.SumRate_OMA = sum(Usr_.Rate_OMA);

sumNOMA = sumNOMA + Usr_.SumRate;
sumOMA = sumOMA + Usr_.SumRate_OMA;



Mark(seq) = Usr_.SumRate;
Mark_OMA(seq) = Usr_.SumRate_OMA;
close all;
end
for (jjj=1:30)
xmark(jjj)=jjj;
end

figure(50);
scatter(xmark,Mark,20,'r','filled');
hold on;
scatter(xmark,Mark_OMA,20,'b','filled');
hold on;

aveNOMA = sumNOMA / seq;
aveOMA = sumOMA / seq;

plot([0,30],[aveNOMA,aveNOMA],'r','linewidth',1);
hold on;
plot([0,30],[aveOMA,aveOMA],'b','linewidth',1);
legend('Sum Rate of NOMA','Sum Rate of OMA','Average Sum Rate of NOMA','Average Sum Rate of OMA')
title('Data Rate of NOMA and OMA in Small Urban Areas')
xlabel('Trials');
ylabel('Sum Data Rate (bit/s)');
grid

run EnergyEff.m;

%% NOMA Computation: EE
SE_noma = aveNOMA/Conf_.BW; % bit/sec/Hz
EE_noma = aveNOMA/EnergyEff_.CP; % bit/watt.sec

%% OMA (OFDMA) Computation: EE
SE_oma = aveOMA/Conf_.BW; % bit/sec/Hz
EE_oma = aveOMA/EnergyEff_.CP; % bit/watt.sec


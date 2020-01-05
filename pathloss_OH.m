function [PL_OkumuraHata] = pathloss_OH(d, fc, hb, hm, flag)
%%
%% d: distance in km
%% flags:   0 --> Small urban areas
%%          1 --> Metropolitan urban areas
%%          2 --> subUrban areas
%%          3 --> Rural areas
%% fc: carrier frequency in MHz
%% hm, hb: antennas height in m
%% Output: PathLoss = A + Blog(d) + C
%%

a_hm_ub = (1.1*log10(fc)-0.7)*hm - (1.56*log10(fc)-0.8);
if (flag == 0)
    % Pathloss for urban: small or medium urban area
    a_hm = a_hm_ub;
    C = 0;
elseif (flag == 1)
    % Pathloss for urban: metropolitan (large cities)
    if (fc <= 200)
        a_hm = 8.29*(log10(1.54*hm))^2 - 1.1;
    elseif (fc >= 400)
        a_hm = 3.2*(log10(11.75*hm))^2 - 4.97;
    end
    C = 0;
elseif (flag == 2)
    % Pathloss for suburban environments
    C = -2*(log10(fc/28))^2 -5.4;
    a_hm = a_hm_ub;
else
    % Pathloss for rural areas
    C = -4.78*(log10(fc))^2 + 18.33*log10(fc)-40.98;
    a_hm = a_hm_ub;
end    
    %% PL = A + Blog(d) + C
    A = 69.55 + 26.16*log10(fc) - 13.82*log10(hb)-a_hm;
    B = 44.9 - 6.55*log10(hb);
    PL_OkumuraHata = A + B*log10(d) + C;
    
    
    
    
%% Group Member Names


%% Lab Section: 


%% SETUP
clear; clc;


%%  Load In Noisy_Finale.wav




%% Identify Frequency of Noise




%% Design Filter




%% Remove Noise




%% Name The Song



%% All functions
function H = DTFT(x,w)
%  ===> Describe function here <===    
    H = zeros(length(w),1);    
    for nn = 1:length(x)        
        H = H + x(nn).*exp(-1j*w.'*(nn-1));    
    end
end


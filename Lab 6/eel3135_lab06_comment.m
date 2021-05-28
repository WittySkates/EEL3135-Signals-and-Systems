%% USER-DEFINED VARIABLES
clear
close all
clc

%% DEFINE FILTER
N = 10;
h = (1/N)*ones(N,1);

% <-- Answer Question: What is the impulse response of this filter? 
%     Use d in place of delta.
% 1/10(d[n] + d[n-1] + d[n-2] + d[n-3] + d[n-4] + d[n-5] + d[n-6] + d[n-7]
% + d[n-8] + d[n-9])

% COMPUTE THE DTFT
n = 0:(N-1);
w = -pi:pi/5000:pi;
H = DTFT(h,w);

% PLOT THE IMPULSE RESPONSE AND DTFT
figure
subplot(3,1,1)
stem(n,h)
xlim([-0.5 20.5])
title('Impulse Response of h')
xlabel('Time Index (n)')
ylabel('Amplitude')
subplot(3,1,2)
plot(w,abs(H))
grid on;
title('Magnitude Response of H')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
subplot(3,1,3)
plot(w,angle(H))
grid on;
title('Phase Response of H')
ylabel('Phase [rad]')
xlabel('Normalized Angular Frequency [rad/s]')


function H = DTFT(x,w)
%  ===> Computes the DTFT of the input <===
  
    H = zeros(length(w),1);
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.'*(nn-1));
    end
    
end

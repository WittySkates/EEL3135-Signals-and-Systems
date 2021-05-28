% Connor Dupuis
% Section: 28944
% TA: Noaki Sawahashi
%% QUESTION 1 COMMENTING

% DO NOT REMOVE THE LINE BELOW
% MAKE SURE 'eel3135_lab06_comment.m' IS IN THE SAME DIRECTORY AS THIS FILE
clear; close all; clc;
type('eel3135_lab06_comment.m')

%% QUESTION 2: DTFT OF COMMON FUNCTIONS 
% COMPUTE THE DTFT
N = 20;
w = -pi:pi/5000:pi;

%% 2 (a) PLOT DTFT
% ALSO ANSWER: Is the data predominantly low frequency, high frequency, 
%              or neither?
% Neither predominantly high or low frequency

n = 0:(N-1);
h = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
H = DTFT(h,w);

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

%% 2 (b) PLOT DTFT
% ALSO ANSWER: Is the data predominantly low frequency, high frequency, 
%              or neither?
% Neither predominantly high or low frequency

n = (0:(N-1));
h = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
H = DTFT(h,w);

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



%%
N = 20;
w = -pi:pi/5000:pi;
x = -pi:pi:pi;
n = (0:(N-1));
h = (cos(x)-1).*exp(-1j.*x);
h1 = ((1/2).^n);
H = DTFT(h,w);

figure
subplot(2,1,1)
plot(w,abs(H))
grid on;
title('Magnitude Response of H')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
subplot(2,1,2)
plot(w,angle(H))
grid on;
title('Phase Response of H')
ylabel('Phase [rad]')
xlabel('Normalized Angular Frequency [rad/s]')


%% 2 (c) PLOT DTFT
% ALSO ANSWER: Is the data predominantly low frequency, high frequency, 
%              or neither?
% Predominantly low frequency

n = (0:(N-1));
h = ((1/2).^n);
H = DTFT(h,w);

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


%% 2 (d) PLOT DTFT
% ALSO ANSWER: Is the data predominantly low frequency, high frequency, 
%              or neither?
% Predominantly high frequency

n = (0:(N-1));
h = (-1/2).^n;
H = DTFT(h,w);

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

%% 2 (e) PLOT DTFT
% ALSO ANSWER: Is the data predominantly low frequency, high frequency, 
%              or neither?
% Predominantly low frequency with small bands over some high frequencies

n = (0:(N-1));
h = [1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
H = DTFT(h,w);

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

%% 2 (f) PLOT DTFT
% ALSO ANSWER: Is the data predominantly low frequency, high frequency, 
%              or neither?
% Predominantly low frequency with small bands over some high frequencies

n = (0:(N-1));
h = cos((pi/4).*n);
H = DTFT(h,w);

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

%% QUESTION 3: DTFT PROPERTIES

%% 3(a) PLOT DTFT
% ALSO ANSWER: describe how each system changes the frequency domain 
% Is the original domain
n = (0:(N-1));
xn = (1 - cos((pi/5).*n)) .* [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
h = xn;
H = DTFT(h,w);

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

%% 3(b) PLOT DTFT
% ALSO ANSWER: describe how each system changes the frequency domain
%
% Only the impusle response is shifted
n = (0:(N-1));
xn = (1 - cos((pi/5).*n)) .* [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
h = xn;
H = DTFT(h,w);

figure
subplot(3,1,1)
stem(n-5,h)
xlim([-5.5 14.5])
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
%% 3(c) PLOT DTFT
% ALSO ANSWER: describe how each system changes the frequency domain
%
% There are bands over the frequency +-pi/2

n = (0:(N-1));
xn = (1 - cos((pi/5).*n)) .* [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
h = xn .*(cos((pi/2).*n));
H = DTFT(h,w);

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

%% 3(d) PLOT DTFT
% ALSO ANSWER: describe how each system changes the frequency domain 
%
% There are evev larger bands over the frequency +-pi/2

n = (0:(N-1));
xn = (1 - cos((pi/5).*n)) .* [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
h = xn .* (xn.*(cos((pi/2).*n)));
H = DTFT(h,w);

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

%% 3(e) PLOT DTFT
% ALSO ANSWER: describe how each system changes the frequency domain
%
% The freqeuncy is centered around 0 with a decline towards +-pi with a
% small spike at +-pi/2 

n = (0:(N-1));
xn = (1 - cos((pi/5).*n)) .* [1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
h = xn + (xn.*(cos((pi/2).*n)));
H = DTFT(h,w);

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

%% QUESTION 4: NULLING FILTER

% DO NOT REMOVE THE LINE BELOW
% MAKE SURE 'NoisyWannabe.wav' IS IN THE SAME DIRECTORY AS THIS FILE
[x, fs] = audioread('Noisy.wav');

%% 4(a) EVALUATE DTFT OF INPUT SIGNAL
w = -pi:pi/5000:pi;
H = DTFT(x,w);

figure
plot(w, abs(H))
grid on
title('Magnitude Response of H')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
%% 4(b) IDENTIFY FREQUENCY

% <== ANSWER TO QUESTION ==>
% The normalized angular noise freqeuncy is roughly 0.08922
% The continuous-time cyclic frequency is 626.2144

w1 = 0.08922;
freq1 = (w1*fs)/(2*pi);
%% 4(c) DESIGN FILTER

h = [1, -2*cos(w1), 1];
N = length(h);
n = (0:(N-1));
H = DTFT(h,w);

figure
subplot(3,1,1)
stem(n,h)
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
%% 4(d) APPLY FILTER

% <== ANSWER TO QUESTION ==>
% There are no longer the massive peaks in magnitude in the filtered audio.
% The plot also follows the filter magnitude as if it were overlayed on top.

filter = conv(h,x);
H = DTFT(filter,w);

figure
plot(w,abs(H))
grid on;
title('Magnitude Response of H')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
%% 4(e) LISTEN TO AUDIO

soundsc(filter,fs);
filter_sc = filter/max(abs(filter));
audiowrite('lab6.wav', filter_sc, fs);

%% ALL FUNCTIONS SUPPORTING THIS CODE 

function H = DTFT(x,w)
%  ===> Computes the DTFT of the input <===
  
    H = zeros(length(w),1);
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.'*(nn-1));
    end
    
end


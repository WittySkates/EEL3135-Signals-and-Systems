%% [Connor Dupuis]
%% [Friday 1:55pm] - [28944] - [Naoki Sawahashi]
%% QUESTION 2 
% DO NOT REMOVE THE LINE BELOW
% MAKE SURE 'eel3135_lab11_comment.m' IS IN THE SAME DIRECTORY AS THIS FILE
clear; close all; clc;
type('eel3135_lab11_comment.m')

%% Question 2
% This question will cover the creation of their own DFT function by
% modifying their DTFT function
% Keep in mind that the DFT can be thought of as a sampled version of the
% DTFT

n=0:59;
x = 0.75 + cos(pi*n/20) + cos(pi*n/15) + cos(pi*n + 2*pi/3);

w_DTFT = linspace(0, 2*pi-pi/5000, 10000);
X_DTFT = DTFT(x, w_DTFT);

% NOTE: USE THE FOLLOWING COMMENTED LINE FOR PLOTTING THE DFT ATOP THE DTFT
% (YOU NEED TO DEFINE w_DFT), THIS WILL MAKE THE PLOTS EASIER TO INTERPRET

%% Question 2a
% Take your DTFT function from previous labs, and modify it to output the
% DFT given an input X

w_DFT = ((2*pi)/length(n))*(n);
X_DFT = DFT(x);

%% Question 2b
% Given x and w_DTFT above, plot the magnitude of the DTFT. Be sure to 
% label axes with units and title the plot.

figure
plot(w_DTFT,abs(X_DTFT));
xlabel('Normalized Frequency (rad/s) (0 - 2pi)');
ylabel('Magnitude')
title('DTFT')
%% Question 2c
% Now plot the magnitude of the 60-length DFT of x. Use hold on
% and hold off to plot the DFT on top of the DTFT plot

figure
plot(w_DTFT,abs(X_DTFT)); 
hold on; plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); hold off;
xlabel('Normalized Frequency (rad/s) (0 - 2pi)');
ylabel('Magnitude')
title('Length 60 DFT')
%% Question 2d
% Now plot the magnitude of the 55-length DFT (i.e. remove the
% last 5 values) of x. Use hold on and hold off to plot the DFT on top of
% the DTFT plot

x5 = x(1:55);

X_DTFT = DTFT(x5, w_DTFT);

w_DFT = ((2*pi)/length(x5))*(0:length(x5)-1);
X_DFT = DFT(x5);

figure
plot(w_DTFT,abs(X_DTFT)); 
hold on; plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); hold off;
xlabel('Normalized Frequency (rad/s) (0 - 2pi)');
ylabel('Magnitude')
title('Length 55 DFT')
%% Question 2e
% Now plot the magnitude of the 65-length DFT (i.e. include 5 
% zeros values) of x. Use hold on and hold off to plot the DFT on top of
% the DTFT plot

x0 = [x zeros(1,5)];

X_DTFT = DTFT(x0, w_DTFT);

w_DFT = ((2*pi)/length(x0))*(0:length(x0)-1);
X_DFT = DFT(x0);

figure
plot(w_DTFT,abs(X_DTFT)); 
hold on; plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); hold off;
xlabel('Normalized Frequency (rad/s) (0 - 2pi)');
ylabel('Magnitude')
title('Length 65 DFT')
%% Question 2f
% Now plot the magnitude of the 200-length DFT of x. Use hold 
% on and hold off to plot the DFT on top of the DTFT plot

x200 = [x zeros(1,200)];

X_DTFT = DTFT(x200, w_DTFT);

w_DFT = ((2*pi)/length(x200))*(0:length(x200)-1);
X_DFT = DFT(x200);

figure
plot(w_DTFT,abs(X_DTFT)); 
hold on; plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); hold off;
xlabel('Normalized Frequency (rad/s) (0 - 2pi)');
ylabel('Magnitude')
title('Length 200 DFT')
%% Question 2g
% Answer in your comments:  Based on the last several questions, what is
% the relationship between the DTFT and the DFT? Under what conditions will
% the theoretical DTFT (i.e. when w_hat is continuous) and DFT have the
% same result?

% The DFT is a sampled DTFT. If the length of the sample (N) is inifinite,
% then the DFT would equal the DTFT.

%% Question 3
% This question will cover the creation of your own DFT function by
% modifying your DTFT function. Keep in mind that the DFT can be thought of
% as a sampled version of the DTFT. In this problem. consider the following
% discrete-time signal: x[n] = u[n-10] - u[n-45]

clear all; clc;
x = [zeros(1,10) ones(1,35)];

%% Question 3a
% Modify the DFT question from Question 2 to create an inverse DFT function
% x = IDFT(x), which inputs a DFT-transformed signal X and outputs the
% time-domain signal x. Hint: the inverse DFT is defined in the lab
% document.

%% Question 3b
% Use stem to plot x[n] for n from 0 to 99, inclusive

x99 = [x zeros(1,55)];
figure
stem(0:length(x99)-1,x99);
xlabel('n');
ylabel('Magnitude')
title('X[n]')
%% Question 3c
% Use conv to compute y[n] = (x[n] convolved with x[n]) and then use stem
% to plot the result for n from 0 to 100, inclusive

xc = conv(x,x);
xc100 = [xc, zeros(1,12)];

figure
stem(0:length(xc100)-1,xc100);
xlabel('Inputs (n)');
ylabel('Magnitude')
title('Convolved X[n]')
%% Question 3d
% Use N = 100 length DFT function to compute y[n] = (x[n] 
% convolved with x[n]) in the frequency domain, then convert back to the
% time domain. Use stem to plot the result.

xd = DFT(x99);
xd = xd.*xd;
xi = IDFT(xd);

figure
stem(0:length(xi)-1, xi);
xlabel('Inputs (n)');
ylabel('Magnitude')
title('Convolved X[n] with DFT and IDFT')
%% Question 3e
% Use N = 60, and repeat 3d

x60 = [x, zeros(1,15)];
xd = DFT(x60);
xd = xd.*xd;
xi = IDFT(xd);

figure
stem(0:length(xi)-1, xi);
xlabel('Inputs (n)');
ylabel('Magnitude')
title('Convolved X[n] with DFT and IDFT')
%% Question 3f
% Answer in your comments: What is the difference in the last two
% solutions? Why does this difference exist?

% In the second problem we are only sampling the first 60 points from the
% convolution. The original convolution is 90 points long, but when we use
% the DFT, we are only using 60. We have to wrap from point 0 to
% point 60, then it would repeat again from point 0.

%% Question 4
% This question will focus on computation time differences between the DFT
% and the FFT. Choose a song at least 3 minutes long to use in this
% problem, and include it in your submission. Load it into MATLAB using
% audioread. Note that mose audio files will be stereo, so you need to make
% Sure that you only use one column of audio data for this part of the lab
% This site has a large archive of free music that you can choose from:
% https://freemusicarchive.org/static
clear all; clc;
[x,fs] = audioread('Wisteria.mp3');
x = x(:,1);

%% Question 4a
% Use DFT to plot the magnitude of the DFT of only the first 10000 samples
% of the audio. Use tic and toc to measure the length of time it takes to 
% compute the DFT. Display the result with the disp function.

TSTART = tic;
xd = DFT(x(1:10000));
figure
stem(0:length(xd)-1, xd);
xlabel('Samples (n)');
ylabel('Magnitude')
title('DFT')
T = toc(TSTART);

fprintf('Elapsed time in seconds for DFT (10000 samples): %d seconds.\n', T);
disp(T);

%% Question 4b
% Use fft (read help fft) to plot the magnitude of the DFT of only the
% first 10000 samles of the audio. Use tic and toc to measure the length of
% time it takes to compute the DFT. Display the result with the disp 
% function.

TSTART = tic;
figure
xf = fft(x,10000);
stem(0:length(xf)-1, xf);
xlabel('Samples (n)');
ylabel('Magnitude')
title('FFT')
T = toc(TSTART);

fprintf('Elapsed time in seconds for FFT (10000 samples): %d seconds.\n', T);
disp(T);

%% Question 4c
% Answer in your comments: Are there any differences in the results? If so,
% why?

% There did not seem to be any difference in magnitude.
%% Question 4d
% Answer in your comments: How much faster is the FFT algorithm compared
% with the DFT in this scenario?

% The FFT is roughly 15-20 times faster than DFT.

%% Question 4e
% Now use fft to plot the magnitude of the DFT of the ENTIRE audio signal.

TSTART = tic;
figure
xf = fft(x);
stem(0:length(xf)-1, xf);
xlabel('Samples (n)');
ylabel('Magnitude')
title('FFT')
T = toc(TSTART);

fprintf('Elapsed time in seconds for FFT (entier song): %d seconds.\n', T);
disp(T);

%% Functions provided for the lab
function H = DTFT(x,w)
% DTFT(X,W)  compute the Discrete-time Fourier Transform of signal X
% acroess frequencies defined by W. 

    H = zeros(1, length(w));
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.*(nn-1));
    end
    
end

function X = DFT(x)
% DFT(x)  compute the N-point Discrete Fourier Transform of signal x  
% Where N is the length of signal x
    N = length(x);
    w = ((2*pi)/N)*(0:N-1);
    X = zeros(1, length(w));  
    
    for nn = 1:N    
        X = X + x(nn).*exp(-1j*w.*(nn-1));    
    end
end

function x = IDFT(X)
% IDFT(x)  compute the N-point Inverse Discrete Fourier Transform of signal
% X where N is the length of signal X
    N = length(X);
    w = ((2*pi)/N)*(0:N-1);
    x = zeros(1, length(w));    
    for nn = 1:length(X)    
        x = x + X(nn).*exp(1j*w.*(nn-1));    
    end
    x=x/N;
end
%% QUESTION 1 COMMENTING

% DO NOT REMOVE THE LINE BELOW 
% MAKE SURE 'eel3135_lab05_comment.m' IS IN SAME DIRECTORY AS THIS FILE
clear; close all;
type('eel3135_lab05_comment.m')

%% QUESTION 2 FREQUENCY FILTERING 

%% 2(a) FILL IN CODE
% ----------- Fill in FreqResponse function down below --------------

%% 2(b) CALCULATE FREQUENCY RESPONSE
w = -pi:(pi/100):pi;
b = [1 2 1];
H = FreqResponse(b,w);
figure;
subplot(2,1,1)
plot(w,abs(H)); 
grid on;
title('Magnitude Response')
xlabel('Normalized Radian Frequency');
ylabel('Amplitude');
subplot(2,1,2)
plot(w,angle(H));
grid on;
title('Phase Response')
xlabel('Normalized Radian Frequency');
ylabel('Phase');

%% 2(c) EVALUATE FREQUENCY RESPONSE FOR CERTAIN FREQUENCIES
disp("w = 0")
H1 = FreqResponse(b,0);
disp(H1);

disp("w = pi/3")
H2 = FreqResponse(b,pi/3);
disp(H2);

disp("w = 9pi/10")
H3 = FreqResponse(b,(9*pi)/10);
disp(H3);

%% 2(d) COMPUTE AND PLOT OUTPUT
n = 0:1:60;
xn = 1 + cos((pi/3)*n) + cos(((9*pi)/10)*n + pi/2);
yn = abs(H1) + (abs(H2)*cos(angle(H2)+ (pi/3)*n)) + (abs(H3)*cos(angle(H3) + ((9*pi)/10)*n + pi/2));

figure;
subplot(2,1,1)
stem(n,xn); 
grid on;
title('Input Signal')
xlabel('Samples');
ylabel('x[n]');
subplot(2,1,2)
stem(n,yn);
grid on;
title('Output Signal')
xlabel('Samples');
ylabel('y[n]');

%% 2(e) COMPARE WITH CONVOLUTION
b = [1 2 1];
zn = conv(xn,b);
n1 = 0:1:62;

figure;
subplot(2,1,1)
stem(n,xn); 
grid on;
title('Input Signal')
xlabel('Samples');
ylabel('x[n]');
subplot(2,1,2)
stem(n1,zn);
grid on;
title('Convolved Signal')
xlabel('Samples');
ylabel('z[n]');
xlim([0 60])

%% 2(f) ANSWER QUESTION
% Besides the first three points, z[n] and y[n] are indentical.


%% QUESTION 3

% DO NOT REMOVE THE LINE BELOW 
% MAKE SURE 'jingle.wav' IS IN SAME DIRECTORY AS THIS FILE
[x, fs] = audioread('jingle11k.wav');

%% 3(a) PLOT FREQUENCY RESPONSE
a = [1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9];
Ha = FreqResponse(a,w);

figure;
subplot(2,1,1)
plot(w,abs(Ha)); 
grid on;
title('Magnitude Response')
xlabel('Normalized Radian Frequency');
ylabel('Amplitude');
subplot(2,1,2)
plot(w,angle(Ha)); 
grid on;
title('Phase Response')
xlabel('Normalized Radian Frequency');
ylabel('Phase');

% <== ANSWER TO QUESTION ==>
% Lowpass filter
%
%% 3(b) APPLY FILTER
soundsc(x,fs);
%%
xa = conv(x,a);
soundsc(xa,fs)

% <==== ANSWER TO QUESTION ====>
% The filter removed higher freqeuncies from the sound.
%
%% 3(c) PLOT FREQUENCY RESPONSE
b = [1, 0, -4, 0, 6, 0, -4, 0, 1];
Hb = FreqResponse(b,w);

figure;
subplot(2,1,1)
plot(w,abs(Hb)); 
grid on;
title('Magnitude Response')
xlabel('Normalized Radian Frequency');
ylabel('Amplitude');
subplot(2,1,2)
plot(w,angle(Hb)); 
grid on;
title('Phase Response')
xlabel('Normalized Radian Frequency');
ylabel('Phase');

% <==== ANSWER TO QUESTION ====>
% Bandpass filter
%
%% 3(d) APPLY FILTER
xb = conv(x,b);
soundsc(xb,fs);

% <==== ANSWER TO QUESTION ====>
% The sound is flatter compared to the original while the freqeuncies that
% are being passed are louder.

%% 3(e) PLOT FREQUENCY RESPONSE
c = conv(b,b);
Hc = FreqResponse(c,w);

figure;
subplot(2,1,1)
plot(w,abs(Hc)); 
grid on;
title('Magnitude Response')
xlabel('Normalized Radian Frequency');
ylabel('Amplitude');
subplot(2,1,2)
plot(w,angle(Hc)); 
grid on;
title('Phase Response')
xlabel('Normalized Radian Frequency');
ylabel('Phase');
% <==== ANSWER TO QUESTION ====>
% Bandpass filter
%
%% 3(f) APPLY FILTER
xc = conv(x,c);
soundsc(xc,fs);

% <==== ANSWER TO QUESTION ====>
% The sound is similar to part e, but even more exaggerated with the
% frequencies that are allowed to pass.

%% ALL FUNCTIONS SUPPORTING THIS CODE %%

function H = FreqResponse(b,w)
%  ===> Describe function here <===
    H = zeros(1,length(w));
    for i = 1:length(b)
        H = H + b(i)*exp(-1j.*w*(i-1));
    end
end

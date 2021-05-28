
%% USER-DEFINED VARIABLES

w = -pi:(pi/100):pi;
% <-- Answer: Why is w from -pi to pi?
% Because it repeats every 2pi.

%% HIGHPASS FILTER

% FREQUENCY RESPONSE
H2 = (1-exp(-1j*w*1));
% <-- Answer: What is the difference equation for this frequency response?

% -------------------> y[n] = x[n] - x[n-1] <----------------------

% PLOT
figure;
subplot(2,1,1)
plot(w,abs(H2)); % ==> What does the abs() function do? <== Returns the absolute value of the passed in parameter
grid on;
title('Magnitude Response')
xlabel('Normalized Radian Frequency');
ylabel('Amplitude');
subplot(2,1,2)
plot(w,angle(H2)); % ==> What does the angle() function do? <== Returns the phase angle in the interval -pi to pi
grid on;
title('Phase Response')
xlabel('Normalized Radian Frequency');
ylabel('Phase');

% <-- Answer: If you input a DC value into a highpass filter, what will be
%             its amplitude?
% It should block the DC value, making the amplitude 0.
% 


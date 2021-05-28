w = -pi:(pi/100):pi;
H2 = (1-exp(-1j*w*1));
H2 = 2j.*sin(w/2).*exp(-1j.*(w/2));
%H2 = (cos(w).*exp(1j.*2.*w));

H = FreqResponse(0,w);


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

function H = FreqResponse(b,w)
%  ===> Describe function here <===
    H = zeros(1,length(w));
    for i = 1:length(b)
        H = H + b(i)*exp(-1j.*w*(i-1));
    end
end
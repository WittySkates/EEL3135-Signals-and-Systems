% Team Members: Sara Kinzbruner, Corinne Meyers, Tom Stowell, 
% Isabella Perlmutter, Connor Dupuis, Rachel Romaine

clear;
close all;

for phi = -pi:1e-3:pi
    amplitude = eel3135_lab01_sinusoid('28944', phi);
    if(amplitude < 5e-5)
        fprintf("Estimated phi: %f\nAmplitude: %f\nPhi0: %f\ntest: %f\n", phi, amplitude, phi-pi+2*pi, mod(phi-pi, 2*pi))
        break
    end
end

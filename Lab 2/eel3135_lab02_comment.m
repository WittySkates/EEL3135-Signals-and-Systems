%% ACKNOWLEDGEMENTS / REFERENCES: 
% This code uses functions written by Ken Schutte in 2019, which is used to
% read and decode midi files. The code is under a GNU General
% Public License, enabling us to run, study, share, and modify the
% software. 
%
% More info can be found at: http://www.kenschutte.com/midi


%% INITIAL SETUP
clear
close all
clc

%% DEFINE MUSIC

% INTITIAL VARIABLES
Fs = 44100;             % ==> The sampling freqeuncy <==

%  ===> Executes the fucntion midiInfo and puts the results into variables Notes and endtime <===
[Notes, endtime] = midiInfo(readmidi('gym.mid'), 0, 2);
L = size(Notes,1);      % ==> Assigns L to the length of the first dimension of Notes <==

%  ===> Passes the function build_song a column vector of ones with length L, column 3 of Notes, column 6 - 5 of notes, and Fs <===
x = build_song(ones(L,1), Notes(:,3), Notes(:,6)-Notes(:,5), Fs);

%  ===> Assigns tot_samples the ceiling of the ((sum of all the elements from the result of columns 6 - 5 of Notes) multiplied by Fs). <===
tot_samples = ceil(sum(Notes(:,6)-Notes(:,5))*Fs);

%  ===> Creates a figure with the plots x vs t one with defualt axis values and one with configured axis value<===
t = 0:1/Fs:(tot_samples-1)/Fs;  % ==> Sets t (the time) equal to the vector of 0 to Fs with increments of tot_samples-1 <==
figure(1); 
subplot(211)
plot(t, x);
xlabel('Time [s]')
ylabel('Amplitude')
subplot(212)
plot(t, x);
xlabel('Time [s]')
ylabel('Amplitude')
axis([0 0.1 -1 1])              % ==> Specifies the limits for the current axes x axis is from 0 to 0.1 and the y axis is form -1 to 1<==

%  ===> Takes any keyboard input to advance to the nextline to exeute soundsc <===
input('Click any button to play sound')
soundsc(x, Fs);



% ========
% YOU DO *NOT* NEED TO DESCRIBE THESE LINES (your free to figure it out though)
W = 0.1;    % Window size
tic;
for mm = 1:ceil(tot_samples/Fs/W)
    % PAUSE UNTIL NEXT FRAME
    xlim([(mm-1)*W+[0 W]]); % Set limits of plot
    tm = toc;                        % Check current time
    if mm*W < tm, disp(['Warning: Visualization is ' num2str(mm*W-tm)  's behind']); end
    drawnow; pause(mm*0.1-tm);       % Synchronize with clock
end
% =======



%%
% =========================================
% SUPPORTING FUNCTIONS FOUND BELOW
% Add comments appropriately below
% =========================================


function x = key_to_note(A, key, dur, fs)
% key_to_note: ========> Takes in a complex amplitude, key, duration, and sampling rate and outputs the sinusoidal waveform of the note <=========
%
% Input Args:
%     A: complex amplitude
%   key: number of the note on piano keyboard
%   dur: duration of each note (in seconds)
%    fs: A scalar sampling rate value
%
% Output:
%     x: sinusoidal waveform of the note
    
    %  ===> Sets N equal to the floor of the duraion multiplied by the sampling rate. 
    % Sets t eqaul to the vector 0 to N-1, then divided by fs. Sets freq value usign the key  <===
    N    = floor(dur*fs);
    t    = (0:(N-1)).'/fs;
    freq = (440/32)*2^((key-9)/12); 
    
    %  ===> Takes the real value of he complex equation <===
    x    = real(A*exp(1j*2*pi*freq*t));   


end


function x = build_song(As, keys, durs, fs)
% build_song:  ========> Takes in arguments As, keys, durs, fs to produce a raw audio signalof length N*fs  <=========
%
% Input Args:
%	  As: A length-N array of complex amplitudes for building notes
%	keys: A length-N array of key numbers (which key on a keyboard) for building notes
%   durs: A length-N array of durations (in seconds) for building notes
%     fs: A scalar sampling rate value
%
% Output Args: 
%      x: A length-(N*fs) length raw audio signal
%
    %  ===> Sets x equal to a column vector of zeros using the duration and sampling rate <===
    x = zeros(ceil(sum(durs)*fs), 1);      
    for k = 1:length(keys) 
        
        %  ===> Sets note equal to the converted key using the helper function key_to_note <===
        note       = key_to_note(As(k), keys(k), durs(k), fs);  
        start_time = sum(durs(1:k-1));
        
        %  ===> Sets n1, n2, and x(n1:n2) to corresponding values using the start time and sampling rate  <===
        n1         = floor(start_time*fs) + 1;               
        n2         = floor(start_time*fs) + floor(durs(k)*fs);
        x(n1:n2)   = x(n1:n2) + note;               
        
    end

end


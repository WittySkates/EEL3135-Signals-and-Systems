%% DEFINE MUSIC
% INTITIAL VARIABLES
Fs = 44100;             % ==> The sampling freqeuncy <==

%  ===> Executes the fucntion midiInfo and puts the results into variables Notes and endtime <===
[Notes, endtime] = midiInfo(readmidi('queen.mid'), 2);
L = size(Notes,1);      % ==> Assigns L to the length of the first dimension of Notes <==

% ===> Passes the function build_song a column vector of ones with length L, column 3 of Notes, column 6 - 5 of notes, and Fs <===
x = build_song_time(ones(L,1), Notes(:,3), Notes(:,5), Notes(:,6), Fs);

soundsc(x,fs)

%% HELPER FUNCTIONS
function x = key_to_note(A, key, dur, fs)
% key_to_note: Produces a sinusoidal waveform corresponding to a 
% 	given piano key number this time with rounding specifically for 
%   build_song_time
%
% Input Args:
%     A: complex amplitude
%   key: number of the note on piano keyboard
%   dur: duration of each note (in seconds)
%    fs: A scalar sampling rate value
%
% Output:
%     x: sinusoidal waveform of the note

    N    = floor(dur*fs);
    t    = (0:(N-1)).'/fs;
    freq = (1500/32)*2^((key-9)/12);
    
    Ak = [0.1155, 0.3417, 0.1789, 0.1232, 0.0678, 0.0473, 0.0260, 0.0045, 0.0020]; % Harmonic amplitudes
    phi = [-2.1299, 1.6727, -2.5454, 0.6607, -2.0390, 2.1597, -1.0467, 1.8581, -2.3925]; % Harmonic phase shifts
    
    % For loop iterating through and summing harmonics
    x = 0;
    for k = 1:length(Ak)
        x = x + Ak(k)*cos(2*pi*k*freq*t + phi(k));
    end
end

function x = build_song_time(As, keys, start_time, end_time, fs)
% build_song: Uses key_to_note and the inputted start and end time to create an output
%   of notes for a specified amount of time.
%
% Input Args:
%          As: A length-N array of complex amplitudes for building notes
%        keys: A length-N array of key numbers (which key on a keyboard) for building notes
%  start_time: A length-N array of start times (in seconds) for notes
%    end_time: A length-N array of end times (in seconds) for notes
%          fs: A scalar sampling rate value
%
% Output Args: 
%      x: A length-(N*fs) length raw audio signal
%
    x = zeros(ceil(end_time(length(end_time)))*fs, 1);
%     x = zeros(5000*fs, 1);
    for k = 1:length(keys)
        durs = end_time(k) - start_time(k);
        note       = key_to_note(As(k), keys(k), durs, fs);
        n1         = floor(start_time(k)*fs) + 1;
        n2         = floor(start_time(k)*fs) + floor(durs*fs);
        x(n1:n2)   = x(n1:n2) + note;
    end
end


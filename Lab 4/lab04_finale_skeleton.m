%% Group Member Names


%% Lab Section: 


%% SETUP
clear; clc;
load('hall.mat');

%Edit the start and end times (in seconds) of your .mp3 file
%(this is to prevent a .wav file that is 1-3 minutes long)

%Ex. I want my file to start playing at the 80 second mark and end
%at the 90 second mark; the resulting .wav file would be 9 seconds long.
start_time = 0;
end_time   = 234;

%name of mp3 file goes here
filename = 'canyoufeel.mp3';

% Make sure your chosen .mp3 file is in your project directory
[x, fs] = convertMp3(filename, start_time, end_time);

%% Create the left and right channels

y_1eft = conv(hh(:,1), x); 
y_right = conv(hh(:,2) ,x); 


%% Combine the two channels into one signal matrix (y) 
%   The output matrix should have two columns
y = [y_1eft, y_right];


%% Listen to and save sounds
soundsc(y, fs);

scaled = y/max(abs(y));

audiowrite('initial_sound.wav', x, fs);
audiowrite('hall_sound.wav', scaled, fs);


%% Question:
% Why would using convolution on each channel of hh individually have the
% same effect as a 2D convolution?

% It has the same effect because the system is time invariant.



%% All functions
function [y, Fs] = convertMp3(mp3filename, start_time, end_time)

% Inputs the name of the mp3 you want to use for your lab.
% The mp3 file must be in your project directory.
% Outputs the data from the .wav file and the sampling frequency
% Note that the .wav file is trimmed to be around 15 seconds long

%Remark: You can download an .mp3 file using YouTube. Just Google a YouTube
%to .mp3 converter and a number of websites should appear.

    [y, Fs] = audioread(mp3filename); %read the mp3 file to obtain fs
    y = y(:,1);
    %Trim the .mp3 so it is an acceptable length - to change where you want
    %to trim, edit the values multiplied by Fs
    
    %y(length(y)-118*Fs:end) = [];
    y(end_time*Fs:end) = []; %trims from 90 seconds in to the end
    y(1:(start_time)*Fs) = []; %trims the beginning 80 seconds

    %Change extension of song file
    songName = mp3filename;
    songName(length(songName)-3:length(songName)) = '.wav';
    
    %convert mp3 into .wav
    audiowrite(songName, y/max(abs(y)), Fs);
    
    % Return song information and fs from .wav file
    [y, Fs] = audioread(songName);
end


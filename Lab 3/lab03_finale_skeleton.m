
% Team Members: Sara Kinzbruner, Corinne Meyers, Tom Stowell,
% Isabella Perlmutter, Connor Dupuis, Rachel Romaine
% Section: 28944
%% DO NOT CHANGE
clear; close all;
z = VideoReader('wheel_video.mp4');
x = read(z);

% The variable x is now a 4 dimentional array, with dimensions 1 and 2 the
% m by n pixels in each color frame. The third dimention is the red,
% green, and blue colors in the image. The fourth dimension represents time


%% You can change this


% (b and c)

Dx = 270/54;
Dy = 270/90;
Dt = 60/3.75;
zs = video_sample(x, Dx, Dy, Dt);


%% Run this section to play and save video (DO NOT CHANGE)

figure(1);
for i = 1:size(zs, 4)
    tic;
    imagesc(uint8(zs(:,:,:,i)));
    axis square;
    tm = toc;
    pause(1/60-tm);
end

v = VideoWriter('output_video','MPEG-4');
open(v)
writeVideo(v,uint8(zs));
close(v)
%% ALL FUNCTIONS SUPPORTING THIS CODE %%

function zs = video_sample(z,Dx, Dy, Dt)
% Change your code from the previous sample function here so that it works
% on videos.
    
    %zs = zeros(ceil(size(z,1)/Dy),ceil(size(z,2)/Dx),size(z,3),ceil(size(z,4)/Dt));
    zs = z(1:Dy:end,1:Dx:end,:,1:Dt:end);
    
end

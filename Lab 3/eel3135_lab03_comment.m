% USER DEFINED VARIABLES
w = 20;           % Width
x = 0:1:79;       % Horiztonal Axis
y = 0:1:79;       % Vertical Axis

% ==> Equation for a 80x80 circle<==
z = round(exp(-1/w.^2*(((y.'-50)/1.5).^2+((x-20)).^2)));

% ==> Applying the two functions to the inputs z and zs respectively <==
[xs,ys,zs] = image_system1(z,3,6);
za         = image_system2(zs,90,-4);

% PLOT RESULT WITH SUBPLOT
figure(1);			
subplot(1,3,1);         % ==> Creating subplot on firgure 1 <==
imagesc(x, y, z);       % ==> Plotting z <==
axis square; axis xy;	% ==> Setting axis values <==
title('Original')
subplot(1,3,2);         % ==> Creating subplot on firgure 1 <==
imagesc(xs, ys, zs);	% ==> Plotting zs <==
axis square; axis xy;	% ==> Setting axis values <==
title('After System 1')
subplot(1,3,3);         % ==> Creating subplot on firgure 1 <==
imagesc(xs, ys, za);	% ==> Plotting za <==
axis square; axis xy;	% ==> Setting axis values <==
title('After System 2')


function [xs, ys, zs] = image_system1(z,Ux,Dy)
%IMAGE_SYSTEM1   ===> Adds Ux zero-valued pixels inserted between each pixel
%                     on the x axis (vertical lines). Samples every Dy pixels
%                     along the y axis.<===

% ==> Creates a new image of zeros of correct size based of off Ux and Dy <==
zs = zeros(ceil(size(z,2)/Dy),ceil(Ux*size(z,1)));

% ==> Creates new borders ys and xs of off correct size <==
ys = 1:ceil(size(z,1)/Dy);
xs = 1:ceil(Ux*size(z,2));

% ==> Asssigns the correct 1 values from Ux and Dy <==
zs(1:end,1:Ux:end) = z(1:Dy:end,1:end);

end

function [za] = image_system2(z,Sx,Sy)
%IMAGE_SYSTEM2   ===> Shifts the image over by Sx pixels and up/down Sy by Sy pixels <===

% ====> Creates a new image of zeros of correct size <====
za = zeros(size(z,1), size(z,2)); 

for nn = 1:size(z,1)
	for mm = 1:size(z,2)
		% ====> Checks boundry conditions so the values are applied in the correct location <====
		if nn > Sy && nn-Sy < size(z,1) && mm > Sx && mm-Sx < size(z,2)
			% ====> Assigns 1 values to correct location based off boundary condititions <====
			za(nn,mm) = 1/2*z(nn-Sy,mm-Sx);
		end
	end
end


end

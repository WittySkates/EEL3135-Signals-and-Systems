% Convolution can be seen as the mathematical expression of how one signal
% can modify the other, you will need to comment the answers to the 
% questions on the following lines of code to demonstrate the basics of 
% convolution. Also answer the question at the end in order to understand 
% the effects that convolution may have on the original signal. 


%% ONE-DIMENSIONAL CONVOLUTION

% INPUT SIGNAL
xx = [1 1 -1 -1 1 1 -1 -1 1 1 1 1 -1 -1 -1 -1]; 
% <-- Answer: What is the length of xx? 
% xx has a length of 16

% FILTER COEFFICIENTS / IMPULSE RESPONSE
bk = [1/4 1/4 1/4 1/4]; 
% <-- Answer: What is the length of bk?
% bk has a length of 4

% FILTER OUTPUT
yy = conv(bk, xx); 
% <-- Answer: What is the length of yy?
% yy has a length of 19

% PERFORM FILTERING IN ALTERNATIVE WAY
zz = 1/4*shift(xx, 0) + 1/4*shift(xx, 1) + 1/4*shift(xx, 2) + 1/4*shift(xx, 3);
% <-- Answer: What is the length of zz?
% zz has a length of 16, but as a column vector

% --> ANSWER BELOW: Explain why yy and zz are different lengths. <--
% yy uses convolution with xx and bk, while zz uses shifts which doesn't
% casue additional inputs
%


% PLOT
figure(1)
subplot(411)
stem(xx); axis([0 20 -1 1])
ylabel('Amplitude')
subplot(412)
stem(bk); axis([0 20 -1 1])
ylabel('Amplitude')
subplot(413)
stem(yy); axis([0 20 -1 1])
ylabel('Amplitude')
subplot(414)
stem(zz); axis([0 20 -1 1])
ylabel('Amplitude')
xlabel('Samples')


%% TWO-DIMENSIONAL CONVOLUTION

% INPUT SIGNAL
x2 = [ones(15) -1*ones(15)]*255; 

% FILTER COEFFICIENTS
b2 = (1/8)*[0 1 1  1  1  0; 0 1 1  1  1  0]; 
b3 =       [1 1 1 -1 -1 -1; 1 1 1 -1 -1 -1]; 

% OUTPUT
y2 = conv2(x2,b2); 
y3 = conv2(x2,b3); 

% --> ANSWER BELOW: Why does filter b2 affect the image as it does? 
%     What are possible applications of filter b2? <--
% Filter b2 blurs the imnage slightly causing the borders to blend. This
% happens due to the fractional positive value of the kernel. This filter
% could be used as a blur in practice.
%

% --> ANSWER BELOW: Why does filter b3 affect the image as it does? 
%     What are possible applications of filter b3? <--
% Filter b3 affects the image due to the negative values causing inversion
% in some of the areas. This filter could be used as an edge detection in
% practice.
%


% PLOT THE FIRST FILTER RESULTS
figure(2) 
subplot(231)
image(x2)
xlabel('x'); ylabel('y'); zlabel('z');
title('x2')
axis equal; axis tight; colormap('gray');
subplot(232)
image(b2)
xlabel('x'); ylabel('y'); zlabel('z');
title('b2')
axis equal; axis tight; colormap('gray');
subplot(233)
image(y2)
xlabel('x'); ylabel('y'); zlabel('z');
title('y2')
axis equal; axis tight; colormap('gray');

% PLOT THE SECOND FILTER RESULTS
subplot(234)
image(x2)
xlabel('x'); ylabel('y'); zlabel('z');
title('x2')
axis equal; axis tight; colormap('gray');
subplot(235)
image(b3)
xlabel('x'); ylabel('y'); zlabel('z');
title('b3')
axis equal; axis tight; colormap('gray');
subplot(236)
image(y3)
xlabel('x'); ylabel('y'); zlabel('z');
title('y3')
axis equal; axis tight; colormap('gray');


function xs = shift(x, s)
%SHIFT   ===> Shifts each elements in the input by s <===

    % ====> Initializes xs <====
    xs = zeros(length(x), 1);
    
    for n = 1:length(x)
        % ====> Sets boundry conditions <====
        if n-s > 0 && n-s < length(x)
            % ====> Assigns values to xs(n) <====
            xs(n) = x(n-s);
        end
    end

end

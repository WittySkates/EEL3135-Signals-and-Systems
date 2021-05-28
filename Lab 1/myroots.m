function r = myroots(n, a)
% myroots: Find all the nth roots of the complex number a
%
% Input Args:
%   n: a positive integer specifying the nth roots
%   a: a complex number whose nth roots are to be returned
%
% Output:
%   r: 1xn vector containing all the nth roots of a 

% n = input("Enter a positive integer specifying the nth root: ");
% a = input("Enter a complex number whose nth roots are to be returned");

poly = -1; % Defines variable -1
poly(n+1) = a; % Puts the complex number into the n+1 place of the poly matrix
r = roots(poly); % Find the roots of the poly matrix
end
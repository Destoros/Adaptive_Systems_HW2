function [y,e,c] = gd_algorithm(x,d,N,mu,Rxx,p,c0) 
% INPUTS: % x ....... input signal vector (column vector) 
% d ....... desired output signal (of same dimensions as x)
% N ....... number of filter coefficients
% mu ...... step-size parameter 
% Rxx ..... autocorrelation matrix 
% p ....... cross-correlation vector (column vector) 
% c0 ...... initial coefficient vector (optional column vector; default all zeros) 
% OUTPUTS: 
% y ....... output signal vector (same length as x) 
% e ....... error signal vector (same length as x) 
% c ....... coefficient matrix (N rows, number of columns = length of x)


if nargin < 7 %check if c0 is given, if not initialize with 0
   c0 = zeros(1,N); 
end

x = x(:); %make sure, it is a column vector
d = d(:);
c0 = c0(:);
p = p(:);

x_pad = [zeros(N-1,1); x];%pad for time instances n < 0
d_pad = [zeros(N-1,1); d];%same for d to keep things in order

y = zeros(size(x_pad)); %create placeholders; after calucation elide the appened zeroes in the beggining
e = zeros(size(x_pad));
c = zeros(N,length(x_pad));


%dont know what this sentence means: Be careful, that in this form the signal statistics are estimated beforehand and not adapted/changed during execution.

c(:,N-1) = c0; %first iteration uses c0, hence we need to write it into c
for n = N:length(x_pad)
    
    x_tap = flip(x_pad(n-N+1:n)); %flip, so the is value at time n is at the top of the vector
    y(n) = c(:,n-1)'*x_tap;
    e(n) = d_pad(n) - y(n); %' means hermitian transposed
        
    c(:,n) = c(:,n-1) + mu*(p - Rxx * c(:,n-1)); %changed to update rule for Gradient Search

end


%now delete the first entries of y,e and c which are zero, to keep the time
%indices in order
y(1:N-1) = [];
e(1:N-1) = [];
c(:,1:N-1) = [];

end


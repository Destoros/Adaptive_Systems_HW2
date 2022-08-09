function [y,e,c] = lms_algorithm(x,d,N,mu,alpha,OPTS,c0)
% INPUTS: % x ....... input signal vector (column vector) 
% d ....... desired output signal (of same dimensions as x) 
% N ....... number of filter coefficients 
% mu ...... step-size parameter 
% alpha ... algorithm dependent parameter 
% OPTS .... 0 for standard LMS, 1 for normalized LMS 
% c0 ...... initial coefficient vector (optional column vector; default all zeros) 
% OUTPUTS: 
% y ....... output signal vector (same length as x) 
% e ....... error signal vector (same length as x) 
% c ....... coefficient matrix (N rows, number of columns = length of x)


%formulas from problem class sheets, page 9

if nargin < 7 %check if c0 is given, if not initialize with 0
   c0 = zeros(1,N); 
end

if OPTS == 1
    norm_x = 1; %if the NMLS was chosen, the norm of the signal energy does not affect the update coefficient value
else
    norm_x = x(:)'*x(:);
end

if ~(0 < mu/norm_x && mu/norm_x < 2) 
    error('Step size mu causes the system to be unstable')
end

%make sure, everything is a column vector
x = x(:); 
d = d(:);
c0 = c0(:);

%pad for time instances n < 0
x_pad = [zeros(N-1,1); x];
d_pad = [zeros(N-1,1); d];%same for d to keep things in order

%create placeholders; after calucation elide the appened zeroes in the beggining
y = zeros(size(x_pad)); 
e = zeros(size(x_pad));
c = zeros(N,length(x_pad));

%intialization for loop
c(:,N-1) = c0; %first iteration uses c0, hence we need to write it into c
mu_calc = mu; %mu for standard LMS; if OPT == 1, it gets overwritten within for loop
for n = N:length(x_pad)
    
    x_tap = flip(x_pad(n-N+1:n));
    y(n) = c(:,n-1)'*x_tap;
    e(n) = d_pad(n) - y(n); %' means hermitian transposed
    
    %change mu depending on chosen OPTS (standard or normalized LMS)
    if OPTS == 1 %normalized LMS
        mu_calc = mu/(alpha + x_tap'*x_tap); %only the energy of the observed current signal   
    end
        
    c(:,n) = c(:,n-1) + mu_calc*conj(e(n))*x_tap; 

end

%now delete the first entries of y,e and c which are zero, to keep the time
%indices in order
y(1:N-1) = [];
e(1:N-1) = [];
c(:,1:N-1) = [];

end


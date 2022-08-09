%Harald Stiegler, 9330054
function c = r_dx(d,x,P)
if (length(d)~=length(x))
    warning('r_dx: length(x)=%d,length(d)=%d! Should x,d be of equal length?',length(x),length(d));
end
d=d(:).';%enforce to be a row vector
x=x(:).';%enforce to be a row vector

d=d(1,1:P);%make vector of length P
x=x(1,1:P);%make vector of length P

c = xcorr(d,x);%correlate and normalize by P
c = 1.0/P*c;
end

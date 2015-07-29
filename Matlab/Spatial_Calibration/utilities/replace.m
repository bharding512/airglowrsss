function new_y = replace(y,x,idx)
% This is simply a convenience function.
% Create a new vector new_y whose i^th element is y(i) if idx(i)==0, or is
% replaced by the next element of x if idx(i)==1.  (i.e. x(sum(idx(1:i))))
% y (Mx1) a vector that will have some values replaced by values in x
% x (Nx1) values to put into y (M >= N)
% idx (Mx1) logical vector representing which values of y should be
%           replaced
new_y = y;
new_y(idx) = x;

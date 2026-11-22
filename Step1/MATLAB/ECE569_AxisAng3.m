function [w,theta] = ECE569_AxisAng3(w_theta)
% Given the product omega times theta,
% this function will return omega and theta separately
%
% Input: w_theta = omega*theta is a 3-vector
% returns: omega is a unit-length 3-vector
%          theta is a real number
%
% You should check that the input is not a zero vector
% If so, set theta=0 and set w=zeros(size(w_theta))

if ~any(w_theta)
    theta = 0;
    w = zeros(size(w_theta));
else
    theta = sqrt(w_theta(1)^2+w_theta(2)^2+w_theta(3)^2);
    w = w_theta/theta;
end
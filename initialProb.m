function [ differential ] = initialProb( t, p2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p1 = 1 - p2;

if t < 100
    v = -100;
else
    v = 5;
end

alpha = .835 * exp (.027 * (v - 35) );	%equivalent to lambda for state 1
beta = .033 * exp (-.093 * (v - 35) );	%equivalent to lambda for state 2


differential = alpha*p1 - beta*p2;

end


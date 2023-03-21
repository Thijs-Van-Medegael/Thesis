function [congestionPrice] = calcPrice(cap,baseload,EVload,previousPrice)
%CALCPRICE Summary of this function goes here
%   Detailed explanation goes here

    w = 4; %Weighting factor to be tuned, decides sensitivity

    congestionPrice = previousPrice + w*max(0,baseload + EVload - cap)/cap; 
end


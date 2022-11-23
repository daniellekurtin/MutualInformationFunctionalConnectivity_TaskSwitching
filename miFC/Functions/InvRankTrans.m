function [newdata] = untitled2(array)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rank = tiedrank(array);
p = rank / ( length(rank) + 1 ); 
newdata = norminv( p, 0, 1 );
end


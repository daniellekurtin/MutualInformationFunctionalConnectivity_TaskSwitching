function d = pointsdist(points, x)
% calculates the distance from each point to x
difs = bsxfun(@minus, points, x);
d = sqrt(sum(difs.^2,1));
end
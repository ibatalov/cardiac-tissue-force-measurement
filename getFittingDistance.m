%% for optimal performane: x1, y1 - shorter vector, x2, y2 - longer vector
function result = getFittingDistance(x1, y1, x2, y2)
result = 0;
for n = 1 : length(x1)
    result = result + min( (x2 - x1(n)).^2 + (y2 - y1(n)).^2);
end
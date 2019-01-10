function cost = getCostFunction(c00, c01, c02, c03, c10, c11, c12, c20, c21, c30, strip_rows, strip_cols)

coefficients = [c00, c01, c02, c03; c10, c11, c12, 0; c20, c21, 0, 0; c30, 0, 0, 0];

myPolyCurve = MyPolyCurve(coefficients);
cost = 0;
for n = 1 : length(strip_rows)
    cost = cost + abs(myPolyCurve.getValue(strip_rows, strip_cols));
end

end
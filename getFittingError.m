function distance = getFittingError(x1, y1, N, L, l, alpha, r0, c0, phi)
%t0 = clock;

[x2, y2] = getFinalStripPoints(N, L, l, alpha, r0, c0, phi);
%round(etime(clock,t0) * 1000)
distance = getFittingDistance(x1, y1, x2, y2);
%round(etime(clock,t0) * 1000)

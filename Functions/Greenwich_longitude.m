function theta_G = Greenwich_longitude(date)

Y   = date(1);
M   = date(2);
D   = date(3);

J0 = 367*Y - floor(7 * (Y + floor((M + 9) / 12)) / 4) + floor(275 * M / 9) + D + 1721013.5;

T0 = (J0 - 2451545) / 36525;
theta_G = 100.4606184 + 36000.77004 * T0 + 0.000387933 * T0^2 - 2.583 * T0^3 * 1e-8;

theta_G = mod(theta_G, 360);

end
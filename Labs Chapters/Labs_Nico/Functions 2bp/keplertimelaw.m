function F = keplertimelaw(mu,e,a,E,t)

n = sqrt(mu/a^3);

F(1) = E - e*sin(E) - n*t;
end

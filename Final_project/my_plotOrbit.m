function my_plotOrbit(kep, mu, npoints, theta_i, theta_f, string)

theta = linspace(theta_i,theta_f,npoints);

R = zeros(3,npoints);

for i = 1:npoints

    kep(6) = theta(i);

    r = kep2car(kep,mu);

    R(:,i) = r;

end

plot3(R(1,:),R(2,:),R(3,:),string,'LineWidth',1.5)

end
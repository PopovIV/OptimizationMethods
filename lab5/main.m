%on border
[x1 ,x2] = meshgrid(-1.1:0.01:1.1);
z = x1.^2+2*x2.^2+exp(x1.^2+x2.^2);
contour(x1, x2, z)
hold on
%1
z = x1.^2 + 2 *x2.^2 - 1;
contour(x1,x2,z,[0,0])
hold on
%2
z = 2 * x1.^2 + x2.^2 - 0.5;
contour(x1,x2,z,[0,0])
hold on
%3
z = x1.^2 - x1;
contour(x1,x2,z,[0,0])
hold on
%answer
plot(0,0,'r*')
%colorbar
shading interp % эта команда сглаживает поверхность

%inside
[x1 ,x2] = meshgrid(-1.1:0.01:1.1);
z = x1.^2+2*x2.^2+exp(x1.^2+x2.^2);
contour(x1, x2, z)
hold on
%1
z = x1.^2 + 2 *x2.^2 - 1;
contour(x1,x2,z,[0,0])
hold on
%2
z = 2 * x1.^2 + x2.^2 - 0.5;
contour(x1,x2,z,[0,0])
hold on
%3
z = x1.^2 - 1;
contour(x1,x2,z,[0,0])
hold on
%answer
plot(0,0,'r*')
%colorbar
shading interp % эта команда сглаживает поверхность

x = -1:0.01:0;
y = x.^6 + 3*x.^2 + 6*x - 1;
[val, idx] = min(y);

plot(x,y)
hold on
plot(x(idx), y(idx), '*r')
xlabel('x')
ylabel('f(x)')
title('x^2+3x^2+6x-1')
legend('x^2+3x^2+6x-1','min', 'Location', 'NorthWest')
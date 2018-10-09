delta = 3;
Y = linspace(0,4,200);
x = zeros(size(Y));

X = Y.^3 - 2*delta*Y.^2 + (delta^2 + 1)*Y;
plot(X(Y<1),Y(Y<1), 'b'), xlim([0,6]), xlabel('X'), ylabel('Y')
hold on 
plot(X(Y >= 1),Y(Y >= 1), 'b--'), xlim([0,6]), xlabel('X'), ylabel('Y')
plot(X(Y >= 1 & Y < 2.794),Y(Y >= 1 & Y < 2.794).^1.3, 'r--'), xlim([0,6]), xlabel('X'), ylabel('Y')
plot(X(Y >= 2.794),Y(Y >= 2.794).^1.3, 'r'), xlim([0,6]), xlabel('X'), ylabel('Y')
hold off
legend('CW (stable)','CW (unstable)','MI (unstable)', 'MI (stable)')
title('Bifurcation diagram for \Delta = 2')
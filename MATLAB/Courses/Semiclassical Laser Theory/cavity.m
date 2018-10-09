delta = [41/30];
Y = linspace(0,4,200);
x = zeros(length(delta), length(Y));

for ind = 1:length(delta)
    X(ind,:) = Y.^3 - 2*delta(ind)*Y.^2 + (delta(ind)^2 + 1)*Y;
end
plot(X(X<1),Y(Y<1)), xlim([0,6]), xlabel('X'), ylabel('Y')
C=num2str(delta');
H=legend(C);
set(H,'Interpreter','none');
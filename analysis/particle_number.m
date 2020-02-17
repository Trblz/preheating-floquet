
t0 = 0;
tf = 4*pi;
n = 1000;

K = 0;
G = 28.28;

[t,chi] = ode45(@(t,x) solvechi(t,x,K,G),linspace(t0,tf,n),[5 0]);

f1 = figure(1);
f1.Position = [100 100 600 600];
plot(t,chi(:,1));
xlabel('$t$','Interpreter','latex','FontSize',23);
ylabel('$\delta\chi_k$','Interpreter','latex','FontSize',23);

w = sqrt(K.^2 + (((G.^2)./2).*(1-cos(2.*t))));
n = (w./2).*(((abs(chi(:,2)).^2)./(w.^2)) + (abs(chi(:,1)).^2)) - 0.5;

f2 = figure(2);
f2.Position = [100 100 600 600];
plot(t,log(n));
xlabel('$t$','Interpreter','latex','FontSize',23);
ylabel('$\log(n_k)$','Interpreter','latex','FontSize',23);

function chi = solvechi(t,x,K,G)

chi = [x(2) ; -(K.^2 + (((G.^2)./2).*(1-cos(2.*t)))).*x(1)];

end

t0 = 1;
tf = 60*pi;
n = 2000;

K = 0;
g = 5000;
m = 2;

[z,chi] = ode45(@(t,x) solvechi(t,x,K,g,m),linspace(t0,tf,n),[0.2 0]);

f1 = figure(1);
f1.Position = [100 100 600 600];
plot(z,chi(:,1));
xlabel('$t$','Interpreter','latex','FontSize',23);
ylabel('$\delta\chi_k$','Interpreter','latex','FontSize',23);

a = z.^(2/3);
phi = 1./(sqrt(3*pi).*m.*z);
w = sqrt((K./a).^2 + (g.*phi.*sin(m.*z)).^2);
np = (w./2).*(((abs(chi(:,2)).^2)./(w.^2)) + (abs(chi(:,1)).^2)) - 0.5;

f2 = figure(2);
f2.Position = [100 100 600 600];
plot(z,log(np),'LineWidth',2);
xlabel('$t$','Interpreter','latex','FontSize',23);
ylabel('$\log(n_k)$','Interpreter','latex','FontSize',23);

%legend('Broad-Band','Narrow-Band');

saveas(f1,'stochastic_long_b.png')
saveas(f2,'stochastic_long.png')

function chi = solvechi(t,x,K,g,m)

a = t.^(2/3);
phi = 1./(sqrt(3*pi).*m.*t);
w2 = (K./a).^2 + (g.*phi.*sin(m.*t)).^2;

chi = [x(2) ; -w2.*x(1)];

end

t0 = 1;
tf = 8*pi;
n = 1000;

K = 0;
G = 155;
% NARROW BAND: K=0.9 G=0.6325
% BROAD BAND: K=0 G=28.28

[z,chi] = ode45(@(t,x) solvechi(t,x,K,G),linspace(t0,tf,n),[2 0]);

f1 = figure(1);
f1.Position = [100 100 600 600];
plot(z,chi(:,1));
xlabel('$t$','Interpreter','latex','FontSize',23);
ylabel('$\delta\chi_k$','Interpreter','latex','FontSize',23);

a = z.^(2/3);
w = sqrt((K./a).^2 + (((G.^2)./2).*(1-cos(2.*z))));
np = (w./2).*(((abs(chi(:,2)).^2)./(w.^2)) + (abs(chi(:,1)).^2)) - 0.5;

f2 = figure(2);
f2.Position = [100 100 600 600];
plot(z,log(np),'LineWidth',2);
xlabel('$t$','Interpreter','latex','FontSize',23);
ylabel('$\log(n_k)$','Interpreter','latex','FontSize',23);


%legend('Broad-Band','Narrow-Band');

%saveas(f2,'comparison.png')

function chi = solvechi(t,x,K,G)

a = t.^(2/3);
chi = [x(2) ; -((K./a).^2 + (((G.^2)./2).*(1-cos(2.*t)))).*x(1)];

end
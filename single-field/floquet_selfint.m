% https://doi.org/10.1016/j.physletb.2019.01.020

m = 1;
t0 = 0;
tf = 20;
n = 1000;
enable_smoothing = 1;

n_res = 300;
T = zeros(1,n_res);
ybounds = [0.1,9];
kbounds = [0,1.5];
y = linspace(ybounds(1),ybounds(2),n_res);
k = linspace(kbounds(1),kbounds(2),n_res);
mu = zeros(n_res,n_res);

for i = 1:n_res
    
    [T(i),~,~] = solvebackg(y(i),m,t0,tf,n);
    
end

if enable_smoothing == 1
    T = smoothdata(T,'sgolay');
end

tic;
for i = 1:n_res
    
    [~,phib,phit] = solvebackg(y(i),m,t0,tf,n);
    
    for j = 1:n_res
        
        [z,O] = ode45(@(t,x) fmatsolve(t,x,k(j),phib,phit,m),linspace(0,T(i),n),[1 0 0 1]);
        mu(n_res-i+1,j) = max((1/T(i)).*log(abs(eig([O(end,1) O(end,2) ; O(end,3) O(end,4)]))));
        
    end
    disp(i);
end
toc;

%%

f1 = figure(1);
f1.Position = [100,100,900,800];

imagesc(mu);
xticks(linspace(1,n_res,6))
yticks(linspace(1,n_res,7))
xticklabels([0 0.3 0.6 0.9 1.2 1.5])
yticklabels([9.0 7.5 6.0 4.5 3.0 1.5 0])

colorbar;
xlabel('$k$','Interpreter','latex','FontSize',22);
ylabel('$\Phi / M$','Interpreter','latex','FontSize',22);

saveas(f1,'floquet_self.png')

%%

function [T,phi,phit] = solvebackg(phi0,m,t0,tf,n)

[t,phib] = ode45(@(t,x) bfield(t,x,m),linspace(t0,tf,n),[phi0,0]);
phi = phib(:,1);
i = find(phi < 0, 1, 'first');

T = 4.*t(i);

phit = t;

end

function phib = bfield(t,x,m)

phib = [x(2) ; -(m.^2).*x(1).*((1+(x(1).^2)).^-0.5)];

end

function O = fmatsolve(t,x,k,y,yt,m)

phi = interp1(yt,y,t);
f = -(k.^2 + (m.^2).*( ((1+(phi.^2))^-0.5) - (phi.^2).*((1+(phi.^2))^-1.5)));

O = [x(3) ; x(4) ; f.*x(1) ; f.*x(2)];

end

% V = (1/4)*(phi^4)

t0 = 0;
tf = 20;
n = 1000;

n_res = 100;
ybounds = [pi+0.1,2*pi-0.1];
kbounds = [0,3.5];
y = linspace(ybounds(1),ybounds(2),n_res);
k = linspace(kbounds(1),kbounds(2),n_res);

T = zeros(1,n_res);
mu = zeros(1,n_res);

for i = 1:n_res
    
    [~,~,T(i)] = bgsolve(t0,tf,n,y(i));
    
end

T = smoothdata(T,'sgolay');
plot(y,T);

for i = 1:n_res
    
    [phit,phib,~] = bgsolve(t0,tf,n,y(i));
    
    for j = 1:n_res
        
        [z,O] = ode45(@(t,x) fmatsolve(t,x,k(j),phib,phit),linspace(t0,T(i),n),[1 0 0 1]);
        mu(n_res-i+1,j) = max((1/T(i)).*log(abs(eig([O(end,1) O(end,2) ; O(end,3) O(end,4)]))));
        
    end
    disp(i)
end

%%

f1 = figure(1);
f1.Position = [100,100,900,800];

imagesc(mu);
xticks(linspace(1,n_res,8))
yticks(linspace(1,n_res,3))
xticklabels([0 0.5 1.0 1.5 2.0 2.5 3.0 3.5])
yticklabels({'2\pi','4/3 \pi','\pi'})

colorbar;
xlabel('$k$','Interpreter','latex','FontSize',22);
ylabel('$\Phi / f$','Interpreter','latex','FontSize',22);

saveas(f1,'floquet_cos.png')

%%

function [phit,phib,T] = bgsolve(t0,tf,n,phi0)

    [phit,phib2] = ode45(@(t,x) bgfield(t,x),linspace(t0,tf,n),[phi0,0]);
    phib = phib2(:,1);
    plot(phit,phib);
    i = find((phib-pi) < 0, 1, 'first');
    T = 4.*phit(i);

end

function eqmot = bgfield(t,x)
    
    eqmot = [x(2) ; sin(x(1))];
    
end

function O = fmatsolve(t,x,k,y,yt)

phi = interp1(yt,y,t);
f = -(k.^2 - cos(phi));

O = [x(3) ; x(4) ; f.*x(1) ; f.*x(2)];

end
% V_int = -1/2*phi*phi*chi*chi

t0 = 0;
tf = pi;
n = 100;

n_res = 1000;
K = linspace(0,3,n_res);
G = linspace(4,0,n_res);
mu = zeros(n_res,n_res);

D = parallel.pool.DataQueue;
afterEach(D,@updateprog);

parfor i = 1:n_res
    for j = 1:n_res
        
        [z,O] = ode45(@(t,x) fmatsolve(t,x,K(i),G(j)),linspace(t0,tf,n),[1 0 0 1]);
        mu(j,i) = max((1/tf).*log(abs(eig([O(end,1) O(end,2) ; O(end,3) O(end,4)]))));

    end
    send(D,i);
    %disp([num2str(i),'/',num2str(n_res)]);
end
%%

f1 = figure(1);
f1.Position = [200 200 900 800];
imagesc(mu);
colorbar;

xticks(linspace(1,n_res,7));
xticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
yticks(linspace(1,n_res,9));
yticklabels({'4.0','3.5','3.0','2.5','2.0','1.5','1.0','0.5','0'})
xlabel('$\frac{\sqrt{(k^2 + m_\chi^2)}}{m_\phi}$','Interpreter','latex','FontSize',23);
ylabel('$\frac{g\Phi}{m_\phi}$','Interpreter','latex','FontSize',23);

saveas(f1,'floquet_phi2chi2.png')

%%

function O = fmatsolve(t,x,K,G)

f = -(K.^2 + ((G.^2/2).*(1+cos(2.*t))));

O = [x(3) ; x(4) ; f.*x(1) ; f.*x(2)];

end

function updateprog(~)

persistent done
if isempty(done)
    done = 1;
end

disp(done);
done = done + 1;

end
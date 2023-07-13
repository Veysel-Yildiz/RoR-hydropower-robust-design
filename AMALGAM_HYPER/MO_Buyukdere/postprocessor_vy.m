%% post processor

d_pars = X(:,1:3);
type = round(X(:,4));
conf = round(X(:,5));

PP = [d_pars type conf ];

for i = 1:numel(type)
    if PP(i,5) ==1 
        PP(i,3) =0;
        d_pars(i,3) =0;
    end
end
clear i

%%
for i = 1:numel(type)
    if PP(i,4) ==1 
        PP2{i} = 'Kaplan';
    elseif PP(i,4) ==2 
        PP2{i} = 'Francis';
    elseif PP(i,4) ==3 
        PP2{i} = 'Pelton';
    end
    
        if PP(i,5) ==1 
        PP3{i} = 'Single';
    elseif PP(i,5) ==2 
        PP3{i} = 'Dual';
    elseif PP(i,5) ==3 
        PP3{i} = 'Triple';
    end
    
end


%%
TurbineType = PP2'; 
TurbineConfig  = PP3'; 

D = d_pars(:,1); OdL = d_pars(:,2); OdS = d_pars(:,3); 
Annual_Cost = F(:,1); Annual_Revenue = -F(:,2); WEP1 = -F(:,3);

tbl = table(TurbineType, TurbineConfig, D,OdL,OdS);

for k = 1: numel(type)
    Bc(k) = Annual_Revenue(k)/ Annual_Cost(k);
end 
Bc = Bc';

tbl = table(TurbineType, TurbineConfig, D,OdL,OdS, Annual_Cost, Annual_Revenue, WEP1, Bc);

%%
figure(10)
plot(Annual_Cost, Annual_Revenue,'o','Color','b','MarkerSize',5,'MarkerFaceColor','#D9FFFF')
xlabel('\textbf{AC [M\$]}','Interpreter','latex','fontsize',14)
%   %zlabel('\textbf{NPV [M\$]}','Interpreter','latex','fontsize',14)
ylabel('\textbf{AR [M\$]}','Interpreter','latex','fontsize',14)
grid on ;

figure(11)
%scatter3(F(:,1),-F(:,2), -F(:,3))
plot3(Annual_Cost,Annual_Revenue, WEP1,'o','Color','b','MarkerSize',5,'MarkerFaceColor','#D9FFFF')
zlabel('\textbf{WEP1 [GWh] }','Interpreter','latex','fontsize',14)
xlabel('\textbf{AC [M\$]}','Interpreter','latex','fontsize',14)
%   %zlabel('\textbf{NPV [M\$]}','Interpreter','latex','fontsize',14)
ylabel('\textbf{AR [M\$]}','Interpreter','latex','fontsize',14)
grid on ;


f = [Annual_Cost Annual_Revenue  WEP1];
figure
F = scatteredInterpolant(f(:,1),f(:,2), f(:,3),'linear','none');

sgr = linspace(min(f(:,1)),max(f(:,1)));
ygr = linspace(min(f(:,2)),max(f(:,2)));
[XX,YY] = meshgrid(sgr,ygr);
ZZ = F(XX,YY);

figure
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(f(:,1),f(:,2),f(:,3),'k.');

zlabel('$f_{dry}$ [GWh] ','Interpreter','latex','fontsize',16)
xlabel('$f_{cost}$ [M\$]','Interpreter','latex','fontsize',16)
%   %zlabel('\textbf{NPV [M\$]}','Interpreter','latex','fontsize',14)
ylabel('$f_{revenue}$  [M\$]','Interpreter','latex','fontsize',16)



figure
subplot(2,2,1)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(f(:,1),f(:,2),f(:,3),'k.');
hold off
subplot(2,2,2)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(f(:,1),f(:,2),f(:,3),'k.');
hold off
view(-148,8)
subplot(2,2,3)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(f(:,1),f(:,2),f(:,3),'k.');
hold off
view(-180,8)
subplot(2,2,4)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(f(:,1),f(:,2),f(:,3),'k.');
hold off
view(-300,8)
%% add best NPV

npv = [1.5829 4.8129 3.0616 2 2];
X(100,:) = npv;

% plot pdf of ALD with different parameters
a = 0;k = 0.5;l = 1;
% lambda
l1 = 1;l2 = 1.5;l3 = 0.5;
% kappa
k1 = 0.5;k2 = 0.2;k3 = 0.8;

x = -5:0.1:5;
y1 = pdfald(x,a,l1,k);
y2 = pdfald(x,a,l2,k);
y3 = pdfald(x,a,l3,k);
y4 = pdfald(x,a,l,k1);
y5 = pdfald(x,a,l,k2);
y6 = pdfald(x,a,l,k3);


figure('position',[1000,100,800,300])
subplot(1,2,1)
hold on
p1 = plot(x,y1,'LineWidth',3,'DisplayName','\lambda=1');
p2 = plot(x,y2,'LineWidth',3,'LineStyle','--','DisplayName','\lambda=1.5');
p3 = plot(x,y3,'LineWidth',3,'LineStyle',':','DisplayName','\lambda=0.5');
xlabel('$x$','Interpreter','latex','FontSize',15)
ylabel('$p(x;\alpha,\lambda ,\kappa )$','Interpreter','latex','FontSize',15)
text(-0.2,y1(59),'\lambda=1')
text(0.7,y2(56),'\lambda=1.5')
text(0.0,y3(66),'\lambda=0.5')
% legend('show','FontSize',15)
subplot(1,2,2)
hold on
p4 = plot(x,y4,'LineWidth',3,'DisplayName','\kappa=0.5');
p5 = plot(x,y5,'LineWidth',3,'LineStyle','--','DisplayName','\lambda=0.2');
p6 = plot(x,y6,'LineWidth',3,'LineStyle',':','DisplayName','\lambda=0.8');
xlabel('$x$','Interpreter','latex','FontSize',15)
ylabel('$p(x;\alpha,\lambda ,\kappa )$','Interpreter','latex','FontSize',15)
text(0.8,y4(56),'\kappa=0.5')
text(2,y5(64),'\kappa=0.2')
text(0.0,y6(66),'\kappa=0.8')

saveas(gcf,'ALD_3_parameters','epsc')
saveas(gcf,'ALD_3_parameters','pdf')






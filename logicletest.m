% Test logicle mapping
% Use example from Park paper
logicle=Logicle(-31.6228,4.5,10000)
x=0:.01:4.5;
y=logicle.unmap(x);
figure(gcf);
clf;
subplot(211);
semilogy(x,y,'r');
hold on;
[zcross,ind]=max(x(y<0))
slope=(y(ind)-y(ind-1))/(x(ind)-x(ind-1))
plot(x,slope*(x-zcross),'g:');
plot(x,logicle.T*10.^(x-logicle.M),'k');
legend('Logicle','Log','Linear');
xlabel('Display Position');
subplot(212);
plot(x(x<3),y(x<3),'r');
c=axis;
hold on;
plot(x,slope*(x-zcross),'g:');
plot(x,logicle.T*10.^(x-logicle.M),'k');
axis(c);
xlabel('Display Position');



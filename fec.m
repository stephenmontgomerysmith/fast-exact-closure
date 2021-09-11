fid = fopen('s.out');
stuff = fscanf(fid,'%g',[7,Inf]);
plot(stuff(1,:),stuff(2,:),'-','Color',[1,0,0]);
xlabel('t')
ylabel('a')
hold on
plot(stuff(1,:),stuff(3,:),'-','Color',[1,1,0]);
plot(stuff(1,:),stuff(4,:),'-','Color',[0,1,0]);
plot(stuff(1,:),stuff(5,:),'-','Color',[0,0,1]);
plot(stuff(1,:),stuff(6,:),'-','Color',[1,0,1]);
plot(stuff(1,:),stuff(7,:),'-','Color',[0,1,1]);
legend('a11','a12','a13','a22','a23','a33')
hold off
print('s-out.eps','-deps');
print('s-out.png');

%%%% exact mean + opt over ET domain
%%%%%%%%  plot output m fn; CBF/CPF;  and uniform over random sets
clear;close all;

xi1Bdry=[3.5056; 5;10;15;20;25;26.6857]; %[0.1;1;10];
m1=[0.0312; 0.1250; 0.4688; 0.1562; 0.1562; 0.0312]; %[0.4;0.4; 0.2];


nD=5;
lb=[0; 3.5056; 1; 0.001; 0.1];
ub=[180; 26.6857; 50; 0.02; 10]; 
 data1=csvread('flowdata_2022-11-11.csv');
 %data1=csvread('flowdata_2023-06-28.csv');
 bbb=find(data1(:,2)>=3.5056 & data1(:,2)<=26.6857);
 x_ob_ori=data1(bbb,1:5);  %3:7
for i=1:5
     x_ob1(:,i)=-1+(x_ob_ori(:,i)-lb(i))*2/(ub(i)-lb(i));
end
for umeankkk=15:15
 y_ob1=data1(bbb,umeankkk);
  data(:,1:nD)=x_ob1; data(:,nD+1)=y_ob1;
  
   
filename=sprintf('mynet_3hair_Nov05_%d.mat', umeankkk);
S=load(filename);
net=S.net;




 x1=12.5;

zz=Appximate_output(net, x1);


figure,
[fff xixi]=ksdensity(zz,'Bandwidth',0.0005, 'Function','pdf');
plot(xixi, fff, 'LineWidth',2,'Color','k');
set(gcf,'Color','w');
set(gca,'fontsize',24);
%xlabel('Mean velocity','fontsize',18);
%ylabel('Empirical PDF','fontsize',18);

% figname=sprintf('figure25hairUmeanPDF_%d.fig',umeankkk);
% saveas(gcf,figname);
figure,
[fff xixi]=ksdensity(zz, 'Bandwidth',0.0005,'Function','cdf');
plot(xixi, fff, 'LineWidth',2,'Color','k');
set(gcf,'Color','w');
set(gca,'fontsize',24);
xlim([0,0.04])
xlabel('Mean velocity','fontsize',18);
ylabel('Empirical CDF','fontsize',18);
end






function zz=Appximate_output(net, xi)

lb1=[0; 3.5056; 1; 0.001; 0.1];
ub1=[180; 26.6857; 50; 0.02; 10]; 
nS=100000;
freq=rand(nS,4);

xx=zeros(nS*length(xi(:,1)),5);
xxtemp= zeros(nS*length(xi(:,1)),4);
fvec=zeros(nS*length(xi(:,1)),1);
zz=zeros(length(xi(:,1)), nS);

for i=1:4
    if i==1
    freq(:,i)=lb1(i)+(freq(:,i))*(ub1(i)-lb1(i));
    else
        freq(:,i)=lb1(i+1)+(freq(:,i))*(ub1(i+1)-lb1(i+1));
    end
end

xx(:,2)=repelem(xi(:,1), nS);
xxtemp(:,1:4)=repmat(freq, length(xi(:,1)),1);
xx(:,1)=xxtemp(:,1);xx(:,3:5)=xxtemp(:,2:4);

for i=1:5
     xx(:,i)=-1+(xx(:,i)-lb1(i))*2/(ub1(i)-lb1(i));
 end

fvec1=forward(net, dlarray(xx, 'BC'));
fvec=extractdata(fvec1);
zz=(reshape(fvec,  nS, length(xi(:,1))))';

Apmean=(mean(zz'))';
 end





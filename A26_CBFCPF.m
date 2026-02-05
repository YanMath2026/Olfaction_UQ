%%%%%%%%  plot output m fn; CBF/CPF  %%% example for leakiness
clear;close all;

xi1Bdry=[3.5056; 5;10;15;20;25;26.6857]; %[0.1;1;10];
m1=[0.0312; 0.1250; 0.4688; 0.1562; 0.1562; 0.0312]; %[0.4;0.4; 0.2];


nD=5;
lb=[0; 3.5056; 1; 0.001; 0.1];
ub=[180; 26.6857; 50; 0.02; 10]; 
 data1=csvread('flowdata_2022-11-11.csv');
 %data1=csvread('flowdata_2023-06-28.csv');
 bbb=find(data1(:,4)>=3.5056 & data1(:,4)<=26.6857);
 x_ob_ori=data1(bbb,3:7);  %3:7


  for i=1:5
     x_ob1(:,i)=-1+(x_ob_ori(:,i)-lb(i))*2/(ub(i)-lb(i));
  end

 y_ob1=data1(bbb,18);
  data(:,1:nD)=x_ob1; data(:,nD+1)=y_ob1;


S=load(['mynet_3hair_Nov05_18.mat']);
net=S.net;



for i=1:6
 
 
        m(i)=m1(i);
        
            x1=(linspace(xi1Bdry(i),xi1Bdry(i+1), 60))'; %40);
       
 
z=Appximate_mean(net, x1);
[minZ_esp(i,1) bind(i)]=min(z); 
[maxZ_esp(i,1) bind1(i)]=max(z);

end



% % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %%% plot m
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
x(1:6,1)=minZ_esp(:,1);x(1:6,2)=maxZ_esp(:,1);x(1:6,3)=m';
x(7,1)=min(x(1:6,1)); x(7,2)=max(x(1:6,2)); x(7,3)=1-sum(m);
for i=1:6
    zz1(1,1)=0.9999*x(i,1);zz1(2:3,1)=x(i,1:2); zz1(4,1)=x(i,2)*1.0001;
    zz1(1,2)=0;zz1(2:3,2)=[x(i,3) 0]; zz1(4,2)=0;
    stairs(zz1(:,1),zz1(:,2),'k','linewidth',2);
    hold on;
end
set(gca,'fontsize',18);
xtemp=sort(x(1:6,1));
set(gca,'xtick',xtemp);
set(gca,'xticklabel',xtemp);
xlabel('Conditional expectation E(Leakiness|Gap)','fontsize',18);
ylabel('Belief mass','fontsize',18);
set(gcf,'Color','w');




%%%%%%%%%%%%%%%%%%%%% plot CBF and CPF
[a1 b1]=sort(x(:,2));
mb=x(b1,3);
[a2 b2]=sort(x(:,1));
mb2=x(b2,3);


xCBF(1)=0.95*a1(1); xCBF(2:length(b1)+1)=a1; xCBF(length(b1)+2)=1.05*a1(end); CBF(1)=0;CBF(length(b1)+2)=1; CPF(length(b1)+2)=1;
xCPF(1)=0.95*a2(1); xCPF(2:length(b1)+1)=a2; xCPF(length(b1)+2)=1.05*a1(end); CPF(1)=0;
for i=1:length(b1)
    CBF(i+1)=sum(mb(1:i));
    CPF(i+1)=sum(mb2(1:i));
end
figure,
stairs(xCBF,CBF,'k-.', 'linewidth',2);
hold on;
stairs(xCPF,CPF,'k','linewidth',2);
set(gca,'fontsize',18);
xlabel('Conditional expectation E(leakiness|Gap)','fontsize',18);
ylabel('CBF and CPF','fontsize',18);
set(gcf,'Color','w');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Apmean=Appximate_mean(net, xi)

lb1=[0; 3.5056; 1; 0.001; 0.1];
ub1=[180; 26.6857; 50; 0.02; 10]; 
nS=100000;
%rand('seed',10); 
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


%%%%%%%%%%%% This code provides surrogate w.r.t 2D inputs; and Sobol'
%%%%%%%%%%%% indices for sensitivity analysis results.
%%%%%%%% with umeankkk=18 for leakiness; umeankkk=15:17 for mean velocity 
%%%%%%%%%%% calculation SA 100 times
clear;close all;

xi1Bdry=[3.5056; 5;10;15;20;25;26.6857]; %[0.1;1;10];
xi1BdryMat = [3.5056 5; 5 10; 10 15; 15 20; 20 25; 25 26.6857];
m1=[0.0312; 0.1250; 0.4688; 0.1562; 0.1562; 0.0312];

nD=5;
lb=[0; 3.5056; 1; 0.001; 0.1];
ub=[180; 26.6857; 50; 0.02; 10]; 
 data1=csvread('flowdata_2022-11-11.csv');
 bbb=find(data1(:,2)>=3.5056 & data1(:,2)<=26.6857);


 x_ob_ori=data1(bbb,1:5);  %3:7
 for i=1:5
     x_ob1(:,i)=-1+(x_ob_ori(:,i)-lb(i))*2/(ub(i)-lb(i));
  end

 for umeankkk=18:18     %%% leakiness

 y_ob1=data1(bbb,umeankkk);
  data2(1:length(y_ob1),1:nD)=x_ob1; data2(1:length(y_ob1),nD+1)=y_ob1;
 

  dataX=dlarray(x_ob1, 'BC');
  dataY=dlarray(y_ob1, 'BC');

 filename=sprintf('mynet_3hair_Nov05_%d.mat', umeankkk);

   S=load(filename);
  net=S.net;
% 
% 
clear TestData;
kki1=1;kki2=2;
x1(1,1:100)=linspace(lb(kki1),ub(kki1),100);
x2(1,1:100)=linspace(lb(kki2),ub(kki2),100);
[p q]=meshgrid(x1(1,:),x2(1,:));
TestData111(:,1)=p(:);
TestData111(:,2)=q(:);
for i=1:5
    TestData(1:10000,i)=0; %(lb(i)+ub(i))/2;
end

    TestData(1:10000,kki1)=-1+(TestData111(:,1)-lb(kki1))*2/(ub(kki1)-lb(kki1));
TestData(1:10000,kki2)=-1+(TestData111(:,2)-lb(kki2))*2/(ub(kki2)-lb(kki2));

TestData=dlarray(TestData,'BC');
yvec=forward(net, TestData);
yvec1=extractdata(yvec);

  for i=1:length(x1)
      for j=1:length(x2)

         y_NN(i,j)=yvec1((i-1)*length(x2)+j);
      end
  end
  figure,
  surf(x2,x1,y_NN, 'EdgeColor', 'none');
 set(gca,'fontsize',18);
xlabel('Gap','fontsize',18)
ylabel('Angle','fontsize',18)
zlabel('Leakiness','fontsize',18);
set(gcf,'Color','w');



for kkkk=1:100
nS=1000000;
xx1(1:nS,1:5)=zeros(nS,5);xx2(1:nS,1:5)=zeros(nS,5); 
xx1(:,1:5)=rand(nS,5);  xx2(:,1:5)=rand(nS,5);




ind1 = find(xx1(:,2)>=0 & xx1(:,2)<=m1(1)); 
temp(ind1,1)= xi1BdryMat(1,1)+rand(length(ind1),1).*(xi1BdryMat(1,2)-xi1BdryMat(1,1));
clear ind1;
for j=2:6
ind1 = find(xx1(:,2)>sum(m1(1:j-1)) & xx1(:,2)<=sum(m1(1:j))); 
temp(ind1,1)= xi1BdryMat(j,1)+rand(length(ind1),1).*(xi1BdryMat(j,2)-xi1BdryMat(j,1)); 
clear ind1;
end
ind1 = find(xx1(:,2)>sum(m1) & xx1(:,2)<=1); 
temp(ind1,1)= xi1BdryMat(1,1)+rand(length(ind1),1).*(xi1BdryMat(6,2)-xi1BdryMat(1,1));
clear ind1;
xx1(:,2)=temp;

ind1 = find(xx2(:,2)>=0 & xx2(:,2)<=m1(1)); 
temp(ind1,1)= xi1BdryMat(1,1)+rand(length(ind1),1).*(xi1BdryMat(1,2)-xi1BdryMat(1,1));
clear ind1;
for j=2:6
ind1 = find(xx2(:,2)>sum(m1(1:j-1)) & xx2(:,2)<=sum(m1(1:j))); 
temp(ind1,1)= xi1BdryMat(j,1)+rand(length(ind1),1).*(xi1BdryMat(j,2)-xi1BdryMat(j,1)); 
clear ind1;
end
ind1 = find(xx2(:,2)>sum(m1) & xx2(:,2)<=1); 
temp(ind1,1)= xi1BdryMat(1,1)+rand(length(ind1),1).*(xi1BdryMat(6,2)-xi1BdryMat(1,1));
clear ind1;
xx2(:,2)=temp;




for i=1:4
   if i==1 
        xx1(:,i)=lb(i)+xx1(:,i)*(ub(i)-lb(i));
        xx2(:,i)=lb(i)+xx2(:,i)*(ub(i)-lb(i));
    else
        xx1(:,i+1)=lb(i+1)+xx1(:,i+1)*(ub(i+1)-lb(i+1));
        xx2(:,i+1)=lb(i+1)+xx2(:,i+1)*(ub(i+1)-lb(i+1));
   end
end


xx1=-1+(xx1-repmat(lb',nS,1))./(ub-lb)'*2;
xx2=-1+(xx2-repmat(lb',nS,1))./(ub-lb)'*2;

TestData1=xx1; TestData2=xx2; 
TestData3=xx2; TestData3(:,1)=xx1(:,1); 
TestData4=xx2; TestData4(:,2)=xx1(:,2);
TestData5=xx2; TestData5(:,3)=xx1(:,3);
TestData6=xx2; TestData6(:,4)=xx1(:,4);
TestData7=xx2; TestData7(:,5)=xx1(:,5); 
TestData8=xx1; TestData8(:,1)=xx2(:,1); 
TestData9=xx1; TestData9(:,2)=xx2(:,2);
TestData10=xx1; TestData10(:,3)=xx2(:,3);
TestData11=xx1; TestData11(:,4)=xx2(:,4);
TestData12=xx1; TestData12(:,5)=xx2(:,5); 



ddata1=forward(net, dlarray(TestData1, 'BC'));
f0=extractdata(ddata1);%(:,4);
f0(find(f0==0))=0.0;
D=sum(f0.^2)/nS-(mean(f0))^2;
var(1:nS, 1:10)=0;

ddata1=forward(net, dlarray(TestData3, 'BC'));%xlsread('Test_Output3.xlsx');
var(:,1)=extractdata(ddata1);%(:,4);
ddata1=forward(net, dlarray(TestData4, 'BC'));%xlsread('Test_Output4.xlsx');
var(:,2)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData5, 'BC'));%xlsread('Test_Output5.xlsx');
var(:,3)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData6, 'BC'));%xlsread('Test_Output6.xlsx');
var(:,4)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData7, 'BC'));%xlsread('Test_Output7.xlsx');
var(:,5)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData8, 'BC'));%xlsread('Test_Output3.xlsx');
var(:,6)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData9, 'BC'));%xlsread('Test_Output4.xlsx');
var(:,7)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData10, 'BC'));%xlsread('Test_Output5.xlsx');
var(:,8)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData11, 'BC'));%xlsread('Test_Output6.xlsx');
var(:,9)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData12, 'BC'));%xlsread('Test_Output7.xlsx');
var(:,10)=extractdata(ddata1);;%(:,4);


for i=1:10
 var(find(var(:,i)==0),i)==0;
 Dx(i)=sum(f0'.*var(:,i))/nS-(mean(f0))^2;   
end

Sol_100(kkkk,1,:)=Dx/D;

% % % % % %%%  random variable with left quadratic within each interval
nS=1000000;
xx1(1:nS,1:5)=zeros(nS,5);xx2(1:nS,1:5)=zeros(nS,5);
xx1(:,1:5)=rand(nS,5); 
xx2(:,1:5)=rand(nS,5);

ind1 = find(xx1(:,2)>=0 & xx1(:,2)<=m1(1)); 
temp(ind1,1)= xi1BdryMat(1,1)+(1-rand(length(ind1),1)).^(1/3)*(xi1BdryMat(1,2)-xi1BdryMat(1,1));
clear ind1;
for j=2:6
ind1 = find(xx1(:,2)>sum(m1(1:j-1)) & xx1(:,2)<=sum(m1(1:j))); 
temp(ind1,1)= xi1BdryMat(j,1)+(1-rand(length(ind1),1)).^(1/3)*(xi1BdryMat(j,2)-xi1BdryMat(j,1)); 
clear ind1;
end
ind1 = find(xx1(:,2)>sum(m1) & xx1(:,2)<=1); 
temp(ind1,1)= xi1BdryMat(1,1)+(1-rand(length(ind1),1)).^(1/3)*(xi1BdryMat(6,2)-xi1BdryMat(1,1));
clear ind1;
xx1(:,2)=temp;

ind1 = find(xx2(:,2)>=0 & xx2(:,2)<=m1(1)); 
temp(ind1,1)= xi1BdryMat(1,1)+(1-rand(length(ind1),1)).^(1/3)*(xi1BdryMat(1,2)-xi1BdryMat(1,1));
clear ind1;
for j=2:6
ind1 = find(xx2(:,2)>sum(m1(1:j-1)) & xx2(:,2)<=sum(m1(1:j))); 
temp(ind1,1)= xi1BdryMat(j,1)+(1-rand(length(ind1),1)).^(1/3)*(xi1BdryMat(j,2)-xi1BdryMat(j,1)); 
clear ind1;
end
ind1 = find(xx2(:,2)>sum(m1) & xx2(:,2)<=1); 
temp(ind1,1)= xi1BdryMat(1,1)+(1-rand(length(ind1),1)).^(1/3)*(xi1BdryMat(6,2)-xi1BdryMat(1,1));
clear ind1;
xx2(:,2)=temp;


for i=1:4
    if i==1
        xx1(:,i)=lb(i)+xx1(:,i)*(ub(i)-lb(i));
        xx2(:,i)=lb(i)+xx2(:,i)*(ub(i)-lb(i));
    else
        xx1(:,i+1)=lb(i+1)+xx1(:,i+1)*(ub(i+1)-lb(i+1));
        xx2(:,i+1)=lb(i+1)+xx2(:,i+1)*(ub(i+1)-lb(i+1));
    end
end




xx1=-1+(xx1-repmat(lb',nS,1))./(ub-lb)'*2;
xx2=-1+(xx2-repmat(lb',nS,1))./(ub-lb)'*2;



TestData1=xx1; TestData2=xx2; 
TestData3=xx2; TestData3(:,1)=xx1(:,1); 
TestData4=xx2; TestData4(:,2)=xx1(:,2);
TestData5=xx2; TestData5(:,3)=xx1(:,3);
TestData6=xx2; TestData6(:,4)=xx1(:,4);
TestData7=xx2; TestData7(:,5)=xx1(:,5);
TestData8=xx1; TestData8(:,1)=xx2(:,1); 
TestData9=xx1; TestData9(:,2)=xx2(:,2);
TestData10=xx1; TestData10(:,3)=xx2(:,3);
TestData11=xx1; TestData11(:,4)=xx2(:,4);
TestData12=xx1; TestData12(:,5)=xx2(:,5); 


ddata1=forward(net, dlarray(TestData1, 'BC'));
f0=extractdata(ddata1);%(:,4);
f0(find(f0==0))=0.0;
D=sum(f0.^2)/nS-(mean(f0))^2;
var(1:nS, 1:10)=0;
ddata1=forward(net, dlarray(TestData3, 'BC'));%xlsread('Test_Output3.xlsx');
var(:,1)=extractdata(ddata1);%(:,4);
ddata1=forward(net, dlarray(TestData4, 'BC'));%xlsread('Test_Output4.xlsx');
var(:,2)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData5, 'BC'));%xlsread('Test_Output5.xlsx');
var(:,3)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData6, 'BC'));%xlsread('Test_Output6.xlsx');
var(:,4)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData7, 'BC'));%xlsread('Test_Output7.xlsx');
var(:,5)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData8, 'BC'));%xlsread('Test_Output3.xlsx');
var(:,6)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData9, 'BC'));%xlsread('Test_Output4.xlsx');
var(:,7)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData10, 'BC'));%xlsread('Test_Output5.xlsx');
var(:,8)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData11, 'BC'));%xlsread('Test_Output6.xlsx');
var(:,9)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData12, 'BC'));%xlsread('Test_Output7.xlsx');
var(:,10)=extractdata(ddata1);;%(:,4);


for i=1:10
 var(find(var(:,i)==0),i)==0;
 Dx(i)=sum(f0'.*var(:,i))/nS-(mean(f0))^2;   
end

Sol_100(kkkk,2,:)=Dx/D;

% % % % % % % % % %%%  random variable with right quadratic within each interval
nS=1000000;
xx1(1:nS,1:5)=zeros(nS,5);xx2(1:nS,1:5)=zeros(nS,5);
xx1(:,1:5)=rand(nS,5); 
xx2(:,1:5)=rand(nS,5);


ind1 = find(xx1(:,2)>=0 & xx1(:,2)<=m1(1)); 
temp(ind1,1)= xi1BdryMat(1,1)+(rand(length(ind1),1)).^(1/3)*(xi1BdryMat(1,2)-xi1BdryMat(1,1));
clear ind1;
for j=2:6
ind1 = find(xx1(:,2)>sum(m1(1:j-1)) & xx1(:,2)<=sum(m1(1:j))); 
temp(ind1,1)= xi1BdryMat(j,1)+(rand(length(ind1),1)).^(1/3)*(xi1BdryMat(j,2)-xi1BdryMat(j,1)); 
clear ind1;
end
ind1 = find(xx1(:,2)>sum(m1) & xx1(:,2)<=1); 
temp(ind1,1)= xi1BdryMat(1,1)+(rand(length(ind1),1)).^(1/3)*(xi1BdryMat(6,2)-xi1BdryMat(1,1));
clear ind1;
xx1(:,2)=temp;

ind1 = find(xx2(:,2)>=0 & xx2(:,2)<=m1(1)); 
temp(ind1,1)= xi1BdryMat(1,1)+(rand(length(ind1),1)).^(1/3)*(xi1BdryMat(1,2)-xi1BdryMat(1,1));
clear ind1;
for j=2:6
ind1 = find(xx2(:,2)>sum(m1(1:j-1)) & xx2(:,2)<=sum(m1(1:j))); 
temp(ind1,1)= xi1BdryMat(j,1)+(rand(length(ind1),1)).^(1/3)*(xi1BdryMat(j,2)-xi1BdryMat(j,1)); 
clear ind1;
end
ind1 = find(xx2(:,2)>sum(m1) & xx2(:,2)<=1); 
temp(ind1,1)= xi1BdryMat(1,1)+(rand(length(ind1),1)).^(1/3)*(xi1BdryMat(6,2)-xi1BdryMat(1,1));
clear ind1;
xx2(:,2)=temp;




for i=1:4
    if i==1
        xx1(:,i)=lb(i)+xx1(:,i)*(ub(i)-lb(i));
        xx2(:,i)=lb(i)+xx2(:,i)*(ub(i)-lb(i));
    else
        xx1(:,i+1)=lb(i+1)+xx1(:,i+1)*(ub(i+1)-lb(i+1));
        xx2(:,i+1)=lb(i+1)+xx2(:,i+1)*(ub(i+1)-lb(i+1));
    end
end


xx1=-1+(xx1-repmat(lb',nS,1))./(ub-lb)'*2;
xx2=-1+(xx2-repmat(lb',nS,1))./(ub-lb)'*2;












TestData1=xx1; TestData2=xx2; 
TestData3=xx2; TestData3(:,1)=xx1(:,1); 
TestData4=xx2; TestData4(:,2)=xx1(:,2);
TestData5=xx2; TestData5(:,3)=xx1(:,3);
TestData6=xx2; TestData6(:,4)=xx1(:,4);
TestData7=xx2; TestData7(:,5)=xx1(:,5); 
TestData8=xx1; TestData8(:,1)=xx2(:,1); 
TestData9=xx1; TestData9(:,2)=xx2(:,2);
TestData10=xx1; TestData10(:,3)=xx2(:,3);
TestData11=xx1; TestData11(:,4)=xx2(:,4);
TestData12=xx1; TestData12(:,5)=xx2(:,5); 


ddata1=forward(net, dlarray(TestData1, 'BC'));
f0=extractdata(ddata1);%(:,4);
f0(find(f0==0))=0.0;
D=sum(f0.^2)/nS-(mean(f0))^2;
var(1:nS, 1:10)=0;
ddata1=forward(net, dlarray(TestData3, 'BC'));%xlsread('Test_Output3.xlsx');
var(:,1)=extractdata(ddata1);%(:,4);
ddata1=forward(net, dlarray(TestData4, 'BC'));%xlsread('Test_Output4.xlsx');
var(:,2)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData5, 'BC'));%xlsread('Test_Output5.xlsx');
var(:,3)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData6, 'BC'));%xlsread('Test_Output6.xlsx');
var(:,4)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData7, 'BC'));%xlsread('Test_Output7.xlsx');
var(:,5)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData8, 'BC'));%xlsread('Test_Output3.xlsx');
var(:,6)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData9, 'BC'));%xlsread('Test_Output4.xlsx');
var(:,7)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData10, 'BC'));%xlsread('Test_Output5.xlsx');
var(:,8)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData11, 'BC'));%xlsread('Test_Output6.xlsx');
var(:,9)=extractdata(ddata1);;%(:,4);
ddata1=forward(net, dlarray(TestData12, 'BC'));%xlsread('Test_Output7.xlsx');
var(:,10)=extractdata(ddata1);;%(:,4);


for i=1:10
 var(find(var(:,i)==0),i)==0;
 Dx(i)=sum(f0'.*var(:,i))/nS-(mean(f0))^2;   
end

Sol(1,:) =Dx/D;
Sol_100(kkkk,3,:)=Dx/D;
end


for i=1:100
for j=1:5
zzz1(i,j)=Sol_100(i,1,j);
zzz2(i,j)=Sol_100(i,2,j);
 zzz3(i,j)=Sol_100(i,3,j);
end
end





for j=1:5
sol(j,1)=median(zzz1(:,j));
sol(j,2)=median(zzz2(:,j));
sol(j,3)=median(zzz3(:,j));
end



 end



figure,
bar(sol)

set(gca,'fontsize',18);
ylabel('Sobol index','fontsize',18);
axis([0 6 -0.0 1])
 xticks([1 2 3 4 5 ]);
xticklabels({'Angle',  'Gap', 'Ant', 'Dist', 'Re'})
xtickangle(45)
set(gcf,'Color','w');
ylim([0,1]);



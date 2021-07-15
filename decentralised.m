 clc;clear;

load data33.mat
load data69.mat
load 'transdata.mat'


  installnode33=[18;22;25;33]

RR33=inv(transpose(M33))*Dr33*inv(M33);
XX33=inv(transpose(M33))*Dx33*inv(M33);
   

     RRnew33(:,1:4)=RR33(:,installnode33-1);
   Aeq33=zeros(33,42);
    Aeq33(1,1:5)=1;
     Aeq33(1,6:7)=1;
      Aeq33(1,8:9)=1;
    
    Aeq33(2:33,2:9)=[RRnew33 RRnew33];
    for k33=2:33
            Aeq33(k33,k33+9)=1;
    end

  beq33(1,1:24)=sum(P33);
 beq33(2:33,1:24)=RR33*P33(2:33,1:24)+XX33*(Q33(2:33,1:24)-QDG33(2:33,1:24));
 
  ub133=[1.1;zeros(4,1);Gendata33(3,5)*ones(2,1)*.8;Gendata33(3,5)*ones(2,1);0;(.05*ones(32,1))];
 lb133=[-1.1;zeros(4,1);Gendata33(3,4)*ones(4,1)*.8;0;(-.05*ones(32,1))];
  ub233=[1.1;Gendata33(1,5)*ones(4,1);Gendata33(3,5)*ones(2,1)*.8;Gendata33(3,5)*ones(2,1);0;(.05*ones(32,1))];
 lb233=[-1.1;Gendata33(1,4)*ones(4,1);Gendata33(3,4)*ones(4,1)*.8;0;(-.05*ones(32,1))];
  ub333=[1.1;zeros(4,1);Gendata33(3,5)*ones(2,1)*.8;Gendata33(3,5)*ones(2,1);0;(.05*ones(32,1))];
 lb333=[-1.1;zeros(4,1);Gendata33(3,4)*ones(4,1)*.8;0;(-.05*ones(32,1))];
 ub33=[ub133,ub233,ub333];
lb33=[lb133,lb233,lb333];

 A33=zeros(2,42);
   A33(1,6:7)=1;
   A33(1,8:9)=1;
   A33(2,6:7)=-1;
    A33(2,8:9)=-1;
b33=[.36*ones(1,24) ;zeros(1,24)];
 
 

H33=zeros(42,42);
for i=1:33
    H33(i+9,i+9)=2
end



  installnode=[2;3;27;64];


RR69=inv(transpose(M69))*Dr69*inv(M69);
XX69=inv(transpose(M69))*Dx69*inv(M69);


 RR69new(:,1:4)=RR69(:,installnode-1);

  Aeq69=zeros(69,78);
  
    Aeq69(1,1:5)=1;
    Aeq69(1,6:7)=1;
    Aeq69(1,8:9)=1;
    Aeq69(2:69,2:9)=[RR69new RR69new];
    for kkk=2:69
            Aeq69(kkk,kkk+9)=1;
    end
  

 beq69(1,1:24)=sum(P69); 
 
  beq69(2:69,1:24)=RR69*P69(2:69,1:24)+XX69*(Q69(2:69,1:24)-QDG69(2:69,1:24));
% (diag(-ones(23,1),-1)+eye(24))




    ub169=[6;zeros(4,1);Gendata69(3,5)*ones(2,1)*.8;Gendata69(3,5)*ones(2,1)*.8;0;(.05*ones(68,1))];
 lb169=[-6;zeros(4,1);Gendata69(3,4)*ones(4,1)*.8;0;(-.05*ones(68,1))];
  ub269=[6;Gendata69(1,5)*ones(4,1);Gendata69(3,5)*ones(2,1)*.8;Gendata69(3,5)*ones(2,1)*.8;0;(.05*ones(68,1))];
 lb269=[-6;Gendata69(1,4)*ones(4,1);Gendata69(3,4)*ones(4,1)*.8;0;(-.05*ones(68,1))];
  ub369=[6;zeros(4,1);Gendata69(3,5)*ones(2,1)*.8;Gendata69(3,5)*ones(2,1)*.8;0;(.05*ones(68,1))];
 lb369=[-6;zeros(4,1);Gendata69(3,4)*ones(4,1)*.8;0;(-.05*ones(68,1))];
 A69=zeros(2,78);
   A69(1,6:7)=1;
    A69(1,8:9)=1;
   A69(2,6:7)=-1;
    A69(2,8:9)=-1;
 b69=[1.2*ones(1,24);-.6*ones(1,24) ];

H69=zeros(78,78);
  
for i=1:69
    H69(i+9,i+9)=4;
end

Xf1=zeros(42,24);
Xf2=zeros(78,24);
Xf3=zeros(42,24);
Xf4=zeros(78,24);

beq69(1,:)=.95.*beq69(1,:);
beq33(1,:)=.95.*beq33(1,:);
 beq=.9.*beq;

    Cost=zeros(10,24);
  Cost1=zeros(10,24);
for iter=1:10
    
  
  for transtime=1:length(LSEtrans(:,3:end))
  
  
    [x(:,transtime),ff(:,transtime),EXITFLAG,OUTPUT(:,transtime),LAMBDA(:,transtime)]=quadprog(G,a,-transpose(Ciq),-biq,transpose(Ceq),beq(:,transtime));

   ff(:,transtime)=transpose(x(:,transtime))*.5*G*x(:,transtime)+transpose(a)*x(:,transtime);
       end
  
  CCC = struct2cell(LAMBDA);
      AAA=cell2mat(CCC);
  LMP=-transpose(AAA(23:27,1:24));
  
  for t=1:24
    
    %node 3
    f3(1:24,1)=[LMP(:,3)];
   f3(t,2:42)=[Gendata33(1,7)*ones(1,4),Gendata33(3,7)*ones(1,4),zeros(33,1)'];
   %node 1
       f1(1:24,1)=[LMP(:,5)];
   f1(t,2:42)=[Gendata33(1,7)*ones(1,4),Gendata33(3,7)*ones(1,4),zeros(33,1)'];
   
end

Intervals33=42;
  for zz1=1:7
 
  options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','iter','StepTolerance',1.0000e-8, 'MaxIterations',200)
   [Xf3(1:Intervals33,zz1),fval3(1,zz1),exitflag1,OUTPUTf3(:,zz1),LAMBDAf3(:,zz1)]=quadprog(H33,f3(zz1,:),A33,b33(1:2,zz1),Aeq33,beq33(1:33,zz1),lb133,ub133);
   [Xf1(1:Intervals33,zz1),fvall(1,zz1),exitflag1,OUTPUTf1(:,zz1),LAMBDAf1(:,zz1)]=quadprog(H33,f1(zz1,:),A33,b33(1:2,zz1),Aeq33,beq33(1:33,zz1),lb133,ub133);
  end
   for tt1=8:20
% 
      [Xf3(1:Intervals33,tt1),fval3(1,tt1),exitflag1,OUTPUTf3(:,tt1),LAMBDAf3(:,tt1)]=quadprog(H33,f3(tt1,:),A33,b33(1:2,zz1),Aeq33,beq33(1:33,tt1),lb233,ub233);
     [Xf1(1:Intervals33,tt1),fvall(1,tt1),exitflag1,OUTPUTf1(:,tt1),LAMBDAf1(:,tt1)]=quadprog(H33,f1(tt1,:),A33,b33(1:2,zz1),Aeq33,beq33(1:33,tt1),lb233,ub233);
%     
  end
%   
  for ff1=21:24
      [Xf3(1:Intervals33,ff1),fval3(1,ff1),exitflag1,OUTPUTf3(:,ff1),LAMBDAf3(:,ff1)]=quadprog(H33,f3(ff1,:),A33,b33(1:2,zz1),Aeq33,beq33(1:33,ff1),lb333,ub333);
       [Xf1(1:Intervals33,ff1),fvall(1,ff1),exitflag1,OUTPUTf1(:,ff1),LAMBDAf1(:,ff1)]=quadprog(H33,f1(ff1,:),A33,b33(1:2,zz1),Aeq33,beq33(1:33,ff1),lb333,ub333);
      
  end       
         
         
      for t=1:24
     %Node 2
    f2(1:24,1)=[LMP(:,2)./10];
   f2(t,2:78)=[Gendata69(1,7)*ones(1,4),Gendata69(3,7)*ones(1,4),zeros(69,1)'];
   %node 4
       f4(1:24,1)=[LMP(:,4)./10];
   f4(t,2:78)=[Gendata69(1,7)*ones(1,4),Gendata69(3,7)*ones(1,4),zeros(69,1)'];
   
      end
      
Intervals69=78;
 for zz2=1:7   
%  options = optimoptions('quadprog','Algorithm','interior-point-convex','StepTolerance',1.0000e-8, 'MaxIterations',200)
   [Xf2(1:Intervals69,zz2),fval2(1,zz2),exitflag,OUTPUTf2(:,zz2),LAMBDAf2(:,zz2)]=quadprog(H69,f2(zz2,:),A69,b69(1:2,zz2),Aeq69,beq69(1:69,zz2),lb169,ub169);
    [Xf4(1:Intervals69,zz2),fval4(1,zz2),exitflag,OUTPUTf2(:,zz2),LAMBDAf2(:,zz2)]=quadprog(H69,f4(zz2,:),A69,b69(1:2,zz2),Aeq69,beq69(1:69,zz2),lb169,ub169);   
  end

   for tt2=8:20
    [Xf2(1:Intervals69,tt2),fval2(1,tt2),exitflag,OUTPUTf2(:,tt2),LAMBDAf2(:,tt2)]=quadprog(H69,f2(tt2,:),A69,b69(1:2,tt2),Aeq69,beq69(1:69,tt2),lb269,ub269);   
     [Xf4(1:Intervals69,tt2),fval4(1,tt2),exitflag,OUTPUTf2(:,tt2),LAMBDAf2(:,tt2)]=quadprog(H69,f4(tt2,:),A69,b69(1:2,tt2),Aeq69,beq69(1:69,tt2),lb269,ub269); 
  end
  
  for ftr=21:24
      [Xf2(1:Intervals69,ftr),fval2(1,ftr),exitflag,OUTPUTf2(:,ftr),LAMBDAf2(:,ftr)]=quadprog(H69,f2(ftr,:),A69,b69(1:2,ftr),Aeq69,beq69(1:69,ftr),lb369,ub369);
       [Xf4(1:Intervals69,ftr),fval4(1,ftr),exitflag,OUTPUTf2(:,ftr),LAMBDAf2(:,ftr)]=quadprog(H69,f4(ftr,:),A69,b69(1:2,ftr),Aeq69,beq69(1:69,ftr),lb369,ub369);

  end 
     beq=[zeros(1,24);Xf2(1,:)*10;Xf3(1,:)*100;Xf4(1,:)*10;Xf1(1,:)*100];
for t=1:24
  Cost(iter,:)=ff./1000;
   Cost1(iter,:)=fval2
end
end

for t=1:24

   f1dist(:,t)=.5*transpose(Xf1(:,t))*H33*Xf1(:,t)+f1(t,:)*Xf1(:,t)
    f2dist(:,t)=.5*transpose(Xf2(:,t))*H69*Xf2(:,t)+f2(t,:)*Xf2(:,t)
     f3dist(:,t)=.5*transpose(Xf3(:,t))*H33*Xf3(:,t)+f3(t,:)*Xf3(:,t)
      f4dist(:,t)=.5*transpose(Xf4(:,t))*H69*Xf4(:,t)+f4(t,:)*Xf4(:,t)
      
      
end

ft=[ff./1000;f1dist;f2dist;f3dist;f4dist]
% 
% save('X','Xf11','Xf22','Xf33','Xf44')
 Pf2opt69(1:24,1:9)=[transpose(Xf2(1:9,:)*10) ]% PU to MW
    Pf4opt69(1:24,1:9)=[transpose(Xf4(1:9,:)*10) ]% PU to MW
    Vbest269=1-Xf2(10:78,:);
   Vbest469=1-Xf4(10:78,:);
 Pf1opt33(1:24,1:9)=[transpose(Xf1(1:9,:)*100) ]% PU to MW
    Pf3opt33(1:24,1:9)=[transpose(Xf3(1:9,:)*100) ]% PU to MW
    Vbest133=1-Xf1(10:42,:);
   Vbest333=1-Xf3(10:42,:);
  newload=beq-[zeros(1,24);sum(Pf1opt33,2)';sum(Pf2opt69,2)';sum(Pf3opt33,2)';sum(Pf4opt69,2)']
   
  Ppv=[zeros(1,24);sum(Pf1opt33(:,1:9),2)';sum(Pf2opt69(:,1:9),2)';sum(Pf3opt33(:,1:9),2)';sum(Pf4opt69(:,1:9),2)']
    %DC OPF results
  %PG1 PG2 PG3 PG4 PG5 sigma2 sigma3 sigma5 sigma5
  X=[transpose(1:24),transpose(x)];

%Branch km multiplier (Thermal limit inequlaity constraint multiplier)
% Hour 12 14 15 23 34 45 21 41 51 32 43 54 
Brkmmulti=[transpose(1:24) transpose(AAA(1:12,:))];

% Lower and Upper Production Inequality Constraint Multipliers for Each Generator

% Hour PG1l PG12  PG13  PG14  PG15 PGUl PGU2  PGU3  PGU4  PGU5
LowUP=[transpose(1:24) transpose(AAA(13:22,:))];

dailycost=sum(ft,2)

P12=b12*(0-X(:,7));
P14=b14*(0-X(:,9));
P15=b15*(0-X(:,10));
P23=b23*(X(:,7)-X(:,8));
P34=b34*(X(:,8)-X(:,9));
P45=b45*(X(:,9)-X(:,10));

Pkm=[P12 P14 P15 P23 P34 P45]


pgen=transpose(x(1:5,:))
 clc;clear;

 
 %%%%%%%%%This code has run on 24th of April 2021  
 
 
load data33.mat
load data69.mat
load 'transdata.mat'





% load 'LMPbefore.mat'

% %first we build the matrixes off all systems and then we combine them to a
% %big one



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
    H69(i+9,i+9)=2;
end

Xf1=zeros(42,24);
Xf2=zeros(78,24);
Xf3=zeros(42,24);
Xf4=zeros(78,24);

%%%BIG optimzation problem parameters 
Aeqt=zeros(209,249);
Aeqt(1:33,1:42)=Aeq33;
Aeqt(34:34+68,43:120)=Aeq69;
Aeqt(103:135,121:162)=Aeq33;
Aeqt(136:204,163:240)=Aeq69;
Aeqt(205:209,241:249)=transpose(Ceq);

At=zeros(30,249);
At(1:2,1:42)=A33;
At(3:4,43:120)=A69;
At(5:6,121:162)=A33;
At(7:8,163:240)=A69;
At(9:30,241:249)=-transpose(Ciq);
bt=[b33;b69;b33;b69;-biq*ones(1,24)];



 lbt1=[lb133;lb169;lb133;lb169;-biq(13:17,1);-2.5*ones(4,1)];
%  lbt1(241:249,1)=0;
 lbt2=[lb233;lb269;lb233;lb269;-biq(13:17,1);-2.5*ones(4,1)];
%  lbt2(241:249,1)=0;
 lbt3=[lb333;lb369;lb333;lb369;-biq(13:17,1);-2.5*ones(4,1)];
%  lbt3(241:249,1)=0;
 ubt1=[ub133;ub169;ub133;ub169;-biq(18:22,1);ones(4,1)];
%   ubt1(241:249,1)=0;
 ubt2=[ub233;ub269;ub233;ub269;-biq(18:22,1);ones(4,1)];
%   ubt2(241:249,1)=0;
 ubt3=[ub333;ub369;ub333;ub369;-biq(18:22,1);ones(4,1)];
%   ubt2(241:249,1)=0;

 Ht=zeros(249,249);
 Ht(1:42,1:42)=H33;
  Ht(43:120,43:120)=H69;
   Ht(121:162,121:162)=H33;
    Ht(163:240,163:240)=H69;
     Ht(241:249,241:249)=G;
     
      LMPt=LMPbefore;
      ft=zeros(24,249);
       for t=1:24
    ft(1:24,1)=[LMPt(:,6)]; % Fcost of generation at at node 5
    ft(t,2:42)=[Gendata33(1,7)*ones(1,4),Gendata33(3,7)*ones(1,4),zeros(33,1)'];
    ft(1:24,43)=[LMPt(:,3)./10];% Fcost of generation at at node 2
   ft(t,44:120)=[Gendata69(1,7)*ones(1,4),Gendata69(3,7)*ones(1,4),zeros(69,1)'];
    ft(1:24,121)=[LMPt(:,4)];% Fcost of generation at at node 3
    ft(t,122:162)=[Gendata33(1,7)*ones(1,4),Gendata33(3,7)*ones(1,4),zeros(33,1)'];
      ft(1:24,163)=[LMPt(:,5)./10];  % Fcost of generation at at node 3
   ft(t,164:240)=[Gendata69(1,7)*ones(1,4),Gendata69(3,7)*ones(1,4),zeros(69,1)'];
   ft(t,241:249)=a;
       end
     

       
      Intervalst=249;
      
  cost=zeros(45,24)
  
   
        for w=101
    w1=(w-1)/100;
    w2=(100-(w-1))/100;
%  w1=.5
%  w2=1
   ftt(1:24,1:240)=w1.* ft(1:24,1:240)
   ftt(1:24,241:249)=w2.*ft(1:24,241:249);
    
   Htt(1:240,1:240)=w1.* Ht(1:240,1:240);
   Htt(241:249,241:249)= w2.* Ht(241:249,241:249); 

     
for iter =1:3
        
       beqt=[beq33;beq69;beq33;beq69;beq];
       
       for zz2=1:7   
%  options = optimoptions('quadprog','Algorithm','interior-point-convex','StepTolerance',1.0000e-8, 'MaxIterations',200)
   [X(1:Intervalst,zz2),fvalt(1,zz2),exitflag1,OUTPUTft(:,zz2),LAMBDAzz2(:,zz2)]=quadprog(Htt,ftt(zz2,:),At,bt(1:30,zz2),Aeqt,beqt(1:209,zz2),lbt1(1:245,:),ubt1(1:245,:));
   
  end
      
   for tt2=8:20
    [X(1:Intervalst,tt2),fvalt(1,tt2),exitflag2,OUTPUTf2(:,tt2),LAMBDAtt2(:,tt2)]=quadprog(Htt,ftt(tt2,:),At,bt(1:30,tt2),Aeqt,beqt(1:209,tt2),lbt2(1:245,:),ubt2(1:245,:));   
 
  end
% %   
  for ftr=21:24
      [X(1:Intervalst,ftr),fvalt(1,ftr),exitflag3,OUTPUTf2(:,ftr),LAMBDAftr(:,ftr)]=quadprog(Htt,ftt(ftr,:),At,bt(1:30,ftr),Aeqt,beqt(1:209,ftr),lbt3(1:245,:),ubt3(1:245,:));
  
  end 

      end
  

for t=1:24
   ftransm(:,t)=.5*transpose(X(241:249,t))*Htt(241:249,241:249)*X(241:249,t)+(ftt(1,241:249))*X(241:249,t)
   f1dist(:,t)=.5*transpose(X(1:42,t))*Htt(1:42,1:42)*X(1:42,t)+ftt(t,1:42)*X(1:42,t)
   f2dist(:,t)=.5*transpose(X(43:120,t))*Htt(43:120,43:120)*X(43:120,t)+ftt(t,43:120)*X(43:120,t)
   f3dist(:,t)=.5*transpose(X(121:162,t))*Htt(121:162,121:162)*X(121:162,t)+ftt(t,121:162)*X(121:162,t)
   f4dist(:,t)=.5*transpose(X(163:240,t))* Htt(163:240,163:240)*X(163:240,t)+ftt(t,163:240)*X(163:240,t)
   
   
    end
                    
  ftottr(1,1:24,w)=ftransm./1000
 ftotf1(1,1:24,w)=f1dist;
 ftotf2(1,1:24,w)=f2dist;
 ftotf3(1,1:24,w)=f3dist;
 ftotf4(1,1:24,w)=f4dist;

         end
        
  
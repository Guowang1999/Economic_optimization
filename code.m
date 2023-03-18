clear;
clc;
fivebus
bus1 = busdata(:,1);
bus2 = busdata(:,2);
bus3 = busdata(:,3);
Theta = busdata(:,4)*pi/180;
PL = busdata(:,5);
QL = busdata(:,6);
PG = busdata(:,7);
QG = busdata(:,8);
Nbus = length(bus1);
v=zeros(Nbus,2);
v(:,1)=bus3;
v(:,2)=[0 0 0 0 0]';
Br1 = linedata(:,1);
Br2 = linedata(:,2);
LR  = linedata(:,3);
LX  = linedata(:,4);
LB  = linedata(:,5);
Tap = linedata(:,6);
Nline= length(Br1);
%------------形成导纳矩阵Y----------------
Y = zeros(Nbus,Nbus);
   %---------导纳矩阵的对角线元素---------
for ii = 1:Nbus
    for jj = 1:Nline
        if Br1(jj) == ii
            Y(ii,ii) = Y(ii,ii)+1/(LR(jj)+LX(jj)*1i)*(1/Tap(jj)^2)+1i*LB(jj);
        elseif Br2(jj) == ii
            Y(ii,ii) = Y(ii,ii)+1/(LR(jj)+LX(jj)*1i)+1i*LB(jj);
        end
    end
end
   %---------导纳矩阵的非对角线元素--------
for ii = 1:Nline
    Y(Br1(ii),Br2(ii)) = Y(Br1(ii),Br2(ii))-1/(LR(ii)+LX(ii)*1i)*(1/Tap(ii));
    Y(Br2(ii),Br1(ii)) = Y(Br1(ii),Br2(ii));
end
%------------导纳矩阵Y-end----------------
G = real(Y);
B = imag(Y);
e = v(:,1);
f = v(:,2);
%-----------利用直角坐标求解潮流------------
P=zeros(Nbus,1);
Q=zeros(Nbus,1);
U2=zeros(Nbus,1);%PV节点电压的平方
for ii=1:Nbus
    if bus2(ii)==2
        U2(ii)=bus3(ii)^2;
    end
end
delta_e=ones(Nbus,1);
delta_f=ones(Nbus,1);
delta_e(5)=0;
delta_f(5)=0;
Delta_PQ=zeros(2*Nbus,1);
delta_P=PG-PL;
delta_Q=QG-QL;
Miter=100;    %迭代次数限制
k=0;
tic;
while  k==0.0001 && max(abs(delta_f))>=0.0001
    k=k+1;
    U3=zeros(Nbus,1);%PV节点电压的平方
    T1_P=zeros(Nbus,1);
    T2_P=zeros(Nbus,1);
    T1=zeros(Nbus,1);
    T2=zeros(Nbus,1);
 for ii=1:Nbus
      if  bus2(ii)==0
        for jj=1:Nbus
           T1_P(ii)=T1_P(ii)+G(ii,jj)*e(jj)-B(ii,jj)*f(jj);
           T2_P(ii)=T2_P(ii)+G(ii,jj)*f(jj)+B(ii,jj)*e(jj);
        end
        P(ii)=T1_P(ii)*e(ii)+T2_P(ii)*f(ii);
        Q(ii)=T1_P(ii)*f(ii)-T2_P(ii)*e(ii);
    elseif bus2(ii)==2
         for jj=1:Nbus
           T1_P(ii)=T1_P(ii)+G(ii,jj)*e(jj)-B(ii,jj)*f(jj);
           T2_P(ii)=T2_P(ii)+G(ii,jj)*f(jj)+B(ii,jj)*e(jj);
         end
        P(ii)=T1_P(ii)*e(ii)+T2_P(ii)*f(ii);
        U3(ii)=e(ii)^2+f(ii)^2;
      end
 end
 for ii=1:Nbus
       if bus2(ii)==0
         Delta_PQ(2*ii-1)=delta_P(ii)-P(ii);
         Delta_PQ(2*ii)=delta_Q(ii)-Q(ii);
     elseif bus2(ii)==2
         Delta_PQ(2*ii-1)=delta_P(ii)-P(ii);
         Delta_PQ(2*ii)=U2(ii)-U3(ii);
     end
 end
        %-----------形成雅克比矩阵JJ-----------
  J=zeros(2*Nbus,2*Nbus);
 for ii=1:Nbus
      if bus2(ii)==0
        for jj=1:Nbus
            if ii~=jj
                J(2*ii-1,2*jj-1)=-B(ii,jj)*e(ii)+G(ii,jj)*f(ii);
                J(2*ii-1,2*jj)=G(ii,jj)*e(ii)+B(ii,jj)*f(ii);
                J(2*ii,2*jj-1)=-J(2*ii-1,2*jj);
                J(2*ii,2*jj)=J(2*ii-1,2*jj-1);
            elseif ii==jj
                for kk=1:Nbus
                    T1(jj)=T1(jj)+G(jj,kk)*f(kk)+B(jj,kk)*e(kk);%T1=b
                    T2(jj)=T2(jj)+G(jj,kk)*e(kk)-B(jj,kk)*f(kk);%T2=a
                end
                J(2*ii-1,2*ii-1)=-B(ii,ii)*e(ii)+G(ii,ii)*f(ii)+T1(ii);
                J(2*ii-1,2*ii)=G(ii,ii)*e(ii)+B(ii,ii)*f(ii)+T2(ii);
                J(2*ii,2*ii-1)=-G(ii,ii)*e(ii)-B(ii,ii)*f(ii)+T2(ii);
                J(2*ii,2*ii)=-B(ii,ii)*e(ii)+G(ii,ii)*f(ii)-T1(ii);
            end
        end
    elseif bus2(ii)==2
        for jj=1:Nbus
            if ii~=jj
                J(2*ii-1,2*jj-1)=-B(ii,jj)*e(ii)+G(ii,jj)*f(ii);
                J(2*ii-1,2*jj)=G(ii,jj)*e(ii)+B(ii,jj)*f(ii);
                J(2*ii,2*jj-1)=0;
                J(2*ii,2*jj)=0;
            elseif ii==jj
                for kk=1:Nbus
                    T1(jj)=T1(jj)+G(jj,kk)*f(kk)+B(jj,kk)*e(kk);%T1=b
                    T2(jj)=T2(jj)+G(jj,kk)*e(kk)-B(jj,kk)*f(kk);%T2=a
                end
                J(2*ii-1,2*ii-1)=-B(ii,ii)*e(ii)+G(ii,ii)*f(ii)+T1(ii);
                J(2*ii-1,2*ii)=G(ii,ii)*e(ii)+B(ii,ii)*f(ii)+T2(ii);
                J(2*ii,2*ii-1)=2*f(ii);
                J(2*ii,2*ii)=2*e(ii);
            end
        end
    end
 end
 Delta_S=Delta_PQ(1:(2*Nbus-2));
 JJ=J(1:(2*Nbus-2),1:(2*Nbus-2));  
      %-----------形成雅克比矩阵JJ--end---------
 Delta_u=JJ\Delta_S;
for ii=1:(Nbus-1)
     delta_e(ii)=Delta_u(2*ii);
     delta_f(ii)=Delta_u(2*ii-1);
 end
    e=delta_e+e;
    f=delta_f+f;
end
u=e+1i*f;U=abs(u);
V=abs(u)*230;
theta=180*imag(log(u))/pi;%节点母线电压相角
04

% PQ潮流法程序编写：

clear;
clc;
fivebus
bus1 = busdata(:,1);
bus2 = busdata(:,2);
bus3 = busdata(:,3);
% Theta = busdata(:,4)*pi/180;
Theta = busdata(:,4);
PL = busdata(:,5);
QL = busdata(:,6);
PG = busdata(:,7);
QG = busdata(:,8);
Nbus = length(bus1);
v=zeros(Nbus,2);
v(:,1)=bus3;
v(:,2)=Theta;
Br1 = linedata(:,1);
Br2 = linedata(:,2);
LR  = linedata(:,3);
LX  = linedata(:,4);
LB  = linedata(:,5);
Tap = linedata(:,6);
Nline= length(Br1);
%------------形成导纳矩阵Y----------------
Y = zeros(Nbus,Nbus);
   %---------导纳矩阵的对角线元素---------
for ii = 1:Nbus
    for jj = 1:Nline
        if Br1(jj) == ii
            Y(ii,ii) = Y(ii,ii)+1/(LR(jj)+LX(jj)*1i)*(1/Tap(jj)^2)+1i*LB(jj);
        elseif Br2(jj) == ii
            Y(ii,ii) = Y(ii,ii)+1/(LR(jj)+LX(jj)*1i)+1i*LB(jj);
        end
    end
end
   %---------导纳矩阵的非对角线元素--------
for ii = 1:Nline
    Y(Br1(ii),Br2(ii)) = Y(Br1(ii),Br2(ii))-1/(LR(ii)+LX(ii)*1i)*(1/Tap(ii));
    Y(Br2(ii),Br1(ii)) = Y(Br1(ii),Br2(ii));
end
%------------导纳矩阵Y-end----------------
G = real(Y);
B = imag(Y);
B1=B(1:(Nbus-1),1:(Nbus-1));
B2=B(1:(Nbus-2),1:(Nbus-2));
U = v(:,1);            %电压幅值
D = v(:,2);            %电压相角
%-----------利用极坐标求解潮流------------
P=PG-PL;
Q=QG-QL;
k=0;   %k为迭代次数
kp=0;  %计算P不平衡量deltaPi的收敛标志(0表示不收敛,1表示收敛)
kq=0;  %计算U不平衡量deltaQi的收敛标志(0表示不收敛,1表示收敛)
deltaPi=zeros(Nbus-1,1);%deltaPi为4阶矩阵
deltaQi=zeros(Nbus-2,1);%deltaQi为pqnum*1阶矩阵
e=0.00001;
tic;
while((kq==0)||(kp==0))&&(k<1000) %效果一样
% while((~kq|~kp)&(k<100))
    k=k+1;
    %-------------------------P--------------------------
    for ii=1:(Nbus-1)
        T1_P=0;
       for jj=1:Nbus
             DD=D(ii)-D(jj);
             T1_P=T1_P+U(ii)*U(jj)*(G(ii,jj)*cos(DD)+B(ii,jj)*sin(DD));
       end
       deltaPi(ii)=P(ii)-T1_P;
    end
    Up=U(1:(Nbus-1));
    B11=inv(B1);
    deltaD=(-B11*(deltaPi./Up))./Up;%求相角a的不平衡量(参考书上式4-59a)
    for ii=1:(Nbus-1)%求相角D的新迭代值矩阵
        D(ii)=D(ii)+deltaD(ii);
    end
    max1=max(abs(deltaPi./Up));%求deltaP/U绝对值的最大值
    if max1<=e   %如果最大值满足要求,则kp置为"1",表示收敛
         kp=1;
    end
   %---------------------P   END---------------------
   %-----------------------Q----------------------
    for ii=1:(Nbus-2)  %书上式4-45(b)求deltaQi
         T1_Q=0;
         for jj=1:Nbus
             DD=D(ii)-D(jj);
             T1_Q=T1_Q+U(ii)*U(jj)*(G(ii,jj)*sin(DD)-B(ii,jj)*cos(DD));
         end
         deltaQi(ii)=Q(ii)-T1_Q;
    end
    Uq=U(1:(Nbus-2));
    B22=inv(B2);
    deltaU=-B22*(deltaQi./Uq);   %求U的不平衡量deltaU
    for ii=1:(Nbus-2)%求幅值U的新迭代值矩阵
        U(ii)=U(ii)+deltaU(ii);
    end
    max2=max(abs(deltaQi./Uq)); %求deltaQ/U绝对值的最大值
    if max2<=e     %如果最大值满足要求,则kp置为"1",表示收敛
         kq=1;
    end
%------------------------------Q   END----------------------------
end    
V=abs(U)*230; %幅值
theta=D*180/pi;

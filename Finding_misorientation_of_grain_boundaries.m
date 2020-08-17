%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Matlab code to find misorientation of grain boundaries from Euler Angles

%%%%%%%%%%%%%%%%%%%%%%%%%%

%input
p1=input('\nEnter phi1 value  (in degrees) of 1st grain:  ');
p1=p1*pi/180;
p=input('\nEnter phi value (in degrees) of 1st grain:  ');
p=p*pi/180;
p2=input('\nEnter phi2 (in degrees) of 1st grain:  ');
p2=p2*pi/180;
fprintf('\n');
q1=input('\nEnter phi1 (in degrees) of 2nd grain:  ');
q1=q1*pi/180;
q=input('\nEnter phi (in degrees) of 2nd grain:  ');
q=q*pi/180;
q2=input('\nEnter phi2 (in degrees) of 2nd grain:  ');
q2=q2*pi/180;
t1=zeros(24);
t2=zeros(24);
t3=zeros(24);
theta=zeros(1,24);
g1=zeros(3,3,1);
g2=zeros(3,3,1);
gp=zeros(3,3,1);
gp1=zeros(3,3,1);
gp2=zeros(3,3,1);
gq=zeros(3,3,1);
gq1=zeros(3,3,1);
gq2=zeros(3,3,1);
m=zeros(3,3,24);
%converting in the form of matrices for both grains
gp1(1,1,1)=cos(p1);
gp1(2,1,1)=-sin(p1);
gp1(1,2,1)=sin(p1);
gp1(2,2,1)=cos(p1);
gp1(3,3,1)=1;
gp2(1,1,1)=cos(p2);
gp2(2,1,1)=-sin(p2);
gp2(1,2,1)=sin(p2);
gp2(2,2,1)=cos(p2);
gp2(3,3,1)=1;
gp(1,1,1)=1;
gp(2,2,1)=cos(p);
gp(2,3,1)=sin(p);
gp(3,2,1)=-sin(p);
gp(3,3,1)=cos(p);
gq1(1,1,1)=cos(q1);
gq1(2,1,1)=-sin(q1);
gq1(1,2,1)=sin(q1);
gq1(2,2,1)=cos(q1);
gq1(3,3,1)=1;
gq2(1,1,1)=cos(q2);
gq2(2,1,1)=-sin(q2);
gq2(1,2,1)=sin(q2);
gq2(2,2,1)=cos(q2);
gq2(3,3,1)=1;
gq(1,1,1)=1;
gq(2,2,1)=cos(q);
gq(2,3,1)=sin(q);
gq(3,2,1)=-sin(q);
gq(3,3,1)=cos(q);
g1=gp2*gp*gp1;
g2=gq2*gq*gq1;

%symmetry matrices considering the 24 symmteries for cubic system

T=zeros(3,3,24);
T(:,:,1)=[1 0 0; 0 1 0; 0 0 1];
T(:,:,2)=[0 0 -1;  0 -1 0; -1 0 0];
T(:,:,3)=[0 0 -1; 0 1 0; 1 0 0];
T(:,:,4)=[-1 0 0; 0 1 0; 0 0 -1];
T(:,:,5)=[0 0 1; 0 1 0; -1 0 0];
T(:,:,6)=[1 0 0; 0 0 -1; 0 1 0];
T(:,:,7)=[1 0 0; 0 -1 0; 0 0 -1];
T(:,:,8)=[1 0 0; 0 0 1; 0 -1 0];
T(:,:,9)=[0 -1 0; 1 0 0; 0 0 1];
T(:,:,10)=[-1 0 0; 0 -1 0; 0 0 1];
T(:,:,11)=[0 1 0; -1 0 0; 0 0 1];
T(:,:,12)=[0 0 1; 1 0 0; 0 1 0];
T(:,:,13)=[0 1 0; 0 0 1; 1 0 0];
T(:,:,14)=[0 0 -1; -1 0 0; 0 1 0];
T(:,:,15)=[0 -1 0; 0 0 1; -1 0 0];
T(:,:,16)=[0 1 0; 0 0 -1; -1 0 0];
T(:,:,17)=[0 0 -1; 1 0 0; 0 -1 0];
T(:,:,18)=[0 0 1; -1 0 0; 0 -1 0];
T(:,:,19)=[0 -1 0; 0 0 -1; 1 0 0];
T(:,:,20)=[0 1 0; 1 0 0; 0 0 -1];
T(:,:,21)=[-1 0 0; 0 0 1; 0 1 0];
T(:,:,22)=[0 0 1; 0 -1 0; 1 0 0];
T(:,:,23)=[0 -1 0; -1 0 0; 0 0 -1];
T(:,:,24)=[-1 0 0; 0 0 -1; 0 -1 0];


%finding the 24 misorientation matrices(also can be calculated for 576 matrices)

for i=1:24
m(:,:,i)=inv(T(:,:,i)*g1(:,:,1))*g2(:,:,1);
t1(i)=m(1,1,i);
t2(i)=m(2,2,i);
t3(i)=m(3,3,i);
theta(i)=acos(0.5*(t1(i)+t2(i)+t3(i)-1));
end

%minimum of 24 angles is taken as miorientation angle

ansRad=min(theta);
ansTheta=ansRad*180/pi;
disp('The misorientation angle is ');
disp(ansTheta);

%finding axis of rotation

for n=1:24
    if theta(n)==ansRad
        if theta(n)==0 || theta(n)==pi
            r1=sqrt(0.5*(m(1,1,n)+1));
            r2=sqrt(0.5*(m(2,2,n)+1));
            r3=sqrt(0.5*(m(3,3,n)+1));
        else
            r1=(m(2,3,n)-m(3,2,n))/(2*sin(ansRad));
            r2=(m(3,1,n)-m(1,3,n))/(2*sin(ansRad));
            r3=(m(1,2,n)-m(2,1,n))/(2*sin(ansRad));
        end
    end
end
disp('The axis of rotation is ');
disp(r1);
disp(r2);
disp(r3);

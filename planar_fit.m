
function[Adv_x,Adv_y,Adv_z,Tilt_x_z_plane,Tilt_y_z_plane,Z_shift] = Eddy_Planar_Fit(Adv_x1,Adv_y1,Adv_z1) 

V           = ones(length(Adv_x1),3);
V(:,1)      = Adv_x1;
V(:,2)      = Adv_y1;
V(:,3)      = Adv_z1;
V1          = V;

%see: Wilczak 2000 SONIC ANEMOMETER TILT CORRECTION ALGORITHMS------
subsamples  = 5000;
bins        = floor(length(Adv_x1)/subsamples);

u=ones(1,bins);
v=ones(1,bins);
w=ones(1,bins);
    
for k=1:1:bins;                               
 
Row                 = k*subsamples;                                
w(k)       = mean(Adv_z1((Row-subsamples+1):Row));         
v(k)       = mean(Adv_y1((Row-subsamples+1):Row));         
u(k)       = mean(Adv_x1((Row-subsamples+1):Row));         
end

% calculate b coefficients
flen=length(u);
su=sum(u); %sums of velocities
sv=sum(v);
sw=sum(w);
suv=sum(u*v'); %sums of velocity products
suw=sum(u*w');
svw=sum(v*w');
su2=sum(u*u');
sv2=sum(v*v');
H=[flen su sv; su su2 suv; sv suv sv2]; %create 3 * 3 matrix
g=[sw;suw;svw]; % transpose of g
x=H\g; %matrix left division

p31=-x(2)/sqrt((x(2))^2+(x(3))^2+1);
p32=-x(3)/sqrt((x(2))^2+(x(3))^2+1);
p33=1/sqrt((x(2))^2+(x(3))^2+1);

C=zeros(3,3);
D=zeros(3,3);
C(1,1)=1;
C(2,2)=(p33/sqrt(p32^2+p33^2));
C(3,3)=(p33/sqrt(p32^2+p33^2));
C(3,2)=(-p32/sqrt(p32^2+p33^2));
C(2,3)=(p32/sqrt(p32^2+p33^2));
D(1,1)=sqrt(p32^2+p33^2);
D(2,2)=1;
D(3,3)=sqrt(p32^2+p33^2);
D(1,3)=p31;
D(3,1)=-p31;
P=D*C;
V1(:,3)=V1(:,3)-x(1);
V1=V1*P;

Adv_x=V1(:,1);
Adv_y=V1(:,2);
Adv_z=V1(:,3);

% Adv_x=Adv_x1;
% Adv_y=Adv_y1;
% Adv_z=Adv_z1;


Z_shift = x(1);
Tilt_x_z_plane = acos(D(1,1))/pi/2*360;
Tilt_y_z_plane  = acos(C(2,2))/pi/2*360;
end

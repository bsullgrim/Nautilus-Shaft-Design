%Grimsley_Project2 11 May 15 - Modiied Goodman Analysis 
%==========================================================================
%This script takes given operating conditions of a propeller shaft for a 
%Nautilus class nuclear submarine and calculates the Shear, Bending Moment,
%Torque, and Angle of Twist along the shaft's axis.

%This code also creates design charts for factor of safety, shaft twist
%ratio, and the mass based on material properties for various diameters

clear 
clc
close all
format compact
format shortg
%% Nautilus Data From Project Description:
Pprop = 5.0e6;  % propeller power requirement (W)       
Pgen  = 100e3;  % Generator Power Requirement (W)
D1    = 0.33;    % Shaft Diameter D1 (m)
D2    = 0.34;    % Shaft Diameter D2 (m)
E     = 200e9;  % Young's Modulus (N/m2)
pois  = 0.3;    % Poisson's Ratio (-)
dens = 7900 ; % material density (kg/m3)
SigmaU = 570e6 ; % ultimate stress (N/m2) 
Nf = 5e6 ;      % Desired Shaft Life (Cycles)

%% Other Given Data:
L1 = 0.5;      % Axial Location 1 (m)
L2 = 0.5;      % Axial Location 2 (m)
L3 = 1.0;      % Axial Location 3 (m)
L4 = 3.0;      % Axial Location 4 (m)
L5 = 2.0;      % Axial Location 5 (m)
L6 = 1.0;      % Axial Location 6 (m)
Rg = 0.4;      % Gear Radius (m)
v  = 43*1000/3600;      % Submarine Linear Speed (km/hr)
w=290*2*pi/60;      %Propeller Shaft Rotational Speed (Rad/s)

%% S-N Data
Nf1 = 10e3;    % 0.9SigmaU
Nf2 = 10e6;    % 0.5SigmaU

%% Initial Calculations:
Tout     = Pprop/w;              % Torque at the Propeller (N-m) 
Tgear    = Pgen/w;               % Torque at the Gear (N-m)
Fz       = Tgear/Rg;             % Force acting on the gear (N)
Fthrust  = Pprop/v;              % Thrust Force (N)
R2z = (-Fz*(L2+L3+L4))/(L2+L3);    % z Reaction Force at Bearing 2 (N)
R1z = -Fz - R2z;                 % z Reaction Force at Bearing 1 (N)
Tin = -(Tout+Tgear);               % Torque at the Engine (N-m)
G = E/(2*(1+pois));              % Modulus of Rigidity (N/m2)
% J1 = (pi/2)*((D1/2)^4);          % Polar Moment of Inertia for D1
% J2 = (pi/2)*((D2/2)^4);          % Polar Moment of Inertia for D2

%External Loads, Shear Forces, & Bending Moments Represented as a Singularity Functions:
%q(z) = R1z<x-L1>^-1 + R2z<x-(L1+L2+L3)>^-1 + Fz<x-(L1+L2+L3+L4)>^-1
%V(z) = R1z<x-L1>^0 + R2z<x-(L1+L2+L3)>^0 + Fz<x-(L1+L2+L3+L4)>^0
%M(z) = R1z<x-L1>^1 + R2z<x-(L1+L2+L3)>^1 + Fz<x-(L1+L2+L3+L4)>^1
D1T = linspace(0.1,0.6,50);
D2T = linspace(0.35,0.55,5);

%% Calculate Shear along the Shaft
x=linspace(0,8,80);   
V=zeros(1,length(x));
V_n=0;
for n=1:80
    if x(n)<L1
        V(n)=0;
   elseif x(n)<(L1+L2+L3)
        V(n)=R1z;
   elseif x(n)<(L1+L2+L3+L4)
        V(n)=R1z+R2z;
         
    else
        V(n)=R1z+R2z+Fz;
       
    end
end
% figure1= figure('position',[50 50 900 900]);
% %% Create text
% 
% annotation(figure1,'textbox',[0.00208768267223382 0.897442455242967 1 0.1],...
%     'String',{'Grimsley Project1 12 Mar 15','This script takes given operating conditions of a propeller shaft for a Nautilus class nuclear','submarine and calculates the Shear, Bending Moment, Torque, and Angle of Twist along the shaft''''s axis'},...
%     'HorizontalAlignment','center',...
%     'FitBoxToText','on',...
%     'EdgeColor','none');
% 
% subplot(2,2,1)
% stairs(x,V)
% grid on
% 
% xlabel('Axial Location (m)')
% ylabel('Shear Force (N)')

%% Calculate Bending Moment along the Shaft
M=zeros(1,length(x));
M_n=0;
for n=1:80
    if x(n)<L1
        M_n=R1z*0;
         M(n)=M_n;   
    elseif x(n)<(L1+L2+L3)
        M_n=R1z*((x(n)-L1));
       M(n)=M_n;   
    elseif x(n) <(L1+L2+L3+L4)
        M_n=R1z*(x(n)-L1)+(R2z*(x(n)-(L1+L2+L3)));
        M(n)=M_n; 
    else
        M_n=(R1z*(x(n)-L1))+(R2z*(x(n)-(L1+L2+L3)))+(Fz*(x(n)-(L1+L2+L3+L4)));
% M_n=0;
         M(n)=M_n;
    end
end

% 
% subplot(2,2,2)
% plot(x,M)
% grid on
% xlabel('Axial Location (m)')
% ylabel('Bending Moment (N-m)')

%% Calculate Torque along the Shaft
T=zeros(1,length(x));
T_n=0;
for n=0:79
    if n<(5+5+10+30)
        T_n=-Tout;
         T(n+1)=T_n;
    elseif n<(5+5+10+30+20+10)
        T_n=-Tout-Tgear;
        T(n+1)=T_n;   
    else
        T_n=-Tout-Tgear-Tin;
        T(n+1)=T_n;  
    end
end
% for n=1:79
%     if x(n)<(L1+L2+L3+L4)
%         T(n)=-Tout;
%     elseif x(n)<(L1+L2+L3+L4+L5+L6)
%         T(n)=-Tout-Tgear;
%     else
%         T(n)=Tout-Tgear-Tin;
%     end
% end
% 
% subplot(2,2,3)
%     stairs(x,T)
%     grid on
% xlabel('Axial Location (m)')
% ylabel('Torque(N-m)')


% % Calculate Twist Angle along the Shaft
% Phi=zeros(1,length(x));
% Phi_n=0;
% for n=1:80
%     if x(n)<L1+L2
%         Phi_n=Phi_n+((-Tout*0.1))/(J1*G);
%          Phi(n)=Phi_n;
%     elseif x(n)<L1+L2+L3+L4+L5
%         Phi_n=Phi_n+((-Tout*0.1))/(J2*G);
%         Phi(n)=Phi_n;   
%     else
%         Phi_n=Phi_n+((-Tgear-Tout)*0.1)/(J1*G);
%         Phi(n)=Phi_n;  
%     end
% end

% 
% 
% subplot(2,2,4)
%     plot(x,Phi)
%     grid on
% xlabel('Axial Location (m)')
% ylabel('Twist Angle (Rad)')
% Xsm = zeros(5,10);
m = zeros(5,10);
for z =1:5
% D2 = D2T(z);


for j = 1:50
%     D1 = D1T(j);
 %% Calculate Twist Angle along the Shaft
 J1 = (pi/2)*((D1T(j)/2)^4);         
J2 = (pi/2)*((D2T(z)/2)^4); 
Phi=zeros(1,length(x));
Phi_n=0;
for n=1:80
    if x(n)<L1+L2
        Phi_n=Phi_n+((-Tout*0.1))/(J1*G);
         Phi(n)=Phi_n;
    elseif x(n)<L1+L2+L3+L4+L5
        Phi_n=Phi_n+((-Tout*0.1))/(J2*G);
        Phi(n)=Phi_n;   
    else
        Phi_n=Phi_n+((-Tgear-Tout)*0.1)/(J1*G);
        Phi(n)=Phi_n;  
    end
end
%% Calculate SigmaX

for n =1:80   
if x(n)<L1+L2
    SigmaX(n) = 32*(M(n))/(pi*D1T(j)^3);
elseif x(n)<L1+L2+L3+L4+L5
    SigmaX(n) = 32*(M(n))/(pi*D2T(z)^3);
else
    SigmaX(n) = 32*(M(n))/(pi*D1T(j)^3);

end
end

%% Calculate TauXZ
% if x<L1+L2
%     TauXZ = (T*16)/(pi*D1^3);
% elseif x<L1+L2+L3+L4+L5
%     TauXZ = (T*16)/(pi*D2^3);
% else
%     TauXZ = (T*16)/(pi*D1^3);
% end

%% Calculate SigmaM
SigmaM = zeros(1,80);

for n = 1:80
if x(n)<L1+L2
    SigmaM(n) = Fthrust/(pi*(D1T(j)/2)^2); 
elseif x(n)<L1+L2+L3
    SigmaM(n) = Fthrust/(pi*(D2T(z)/2)^2); 
else
    SigmaM(n) = 0;
end
end

%% Calculate SigmaA
% SigmaA = (1/sqrt(2))*(((SigmaZ.^2)+(SigmaZ.^2)+ 6*(TauXZ)).^0.5);
%% Calculate SigmaAR
SigmaAR = SigmaX./(1-(SigmaM./SigmaU));
% figure(2)
% plot(x,SigmaAR)
%% B & A Coefficients 
B = (log((0.9*SigmaU)/(0.5*SigmaU)))/(log(Nf1/Nf2));
A = (0.9*SigmaU)/(Nf1^B);

%% Stress due to expected life
Sigma1 = A*(Nf^B);

%% Factor of Safety
Xs = Sigma1./SigmaAR;

% XSs = Xs(1:50);
% S = x(1:50);
Xsm(z,j)= min(Xs(1:50));
%% Calculate Shaft Mass
m(z,j) = dens*pi*((((D1T(j)/2)^2)*(L1+L2+L6))+((((D2T(z)/2)^2)*(L3+L4+L5))));


%% R Values
PC1 = 0.0175*L1/(20*D1T(j));
PC2 = 0.0175*(L2)/(20*D1T(j));
PC3 = 0.0175*(L3)/(20*D2T(z));
PC4 = 0.0175*(L4)/(20*D2T(z));
PC5 = 0.0175*(L5)/(20*D2T(z));
PC6 = 0.0175*L6/(20*D1T(j));

P1 = max(abs(Phi(1:5)-Phi(1)));
P2 = max(abs(Phi(5:10)-Phi(5)));
P3 = max(abs(Phi(10:20)-Phi(10)));
P4 = max(abs(Phi(20:50)-Phi(20)));
P5 = max(abs(Phi(50:70)-Phi(50)));
P6 = max(abs(Phi(70:80)-Phi(70)));

R1 = PC1/(P1);
R2 = PC2/(P2);
R3 = PC3/(P3);
R4 = PC4/(P4);
R5 = PC5/(P5);
R6 = PC6/(P6);

Rarray = [R1 R2 R3 R4 R5 R6];
R(z,j) = min(Rarray);

end


end
s1 = num2str(D2T(1));
s2 = num2str(D2T(2));
s3 = num2str(D2T(3));
s4 = num2str(D2T(4));
s5 = num2str(D2T(5));
% 
% figure(3)
% plot(D1T,Xsm,'linewidth',2)
% hold on
% title ('Nautilus Design Data','fontweight','bold')
% xlabel('D1(m)','fontweight','bold')
% ylabel('Minimum Factor of Safety, Xs','fontweight','bold')
% g=line('XData',[0.1 0.6],'YData',[3 3],'LineStyle','-.','Linewidth',3,'Color','k');
% h=legend('show','location','best',s1,s2,s3,s4,s5,'Desired Xs');
% htitle = get(h,'Title');
% set(htitle,'string','D2(m)','fontweight','bold');
% grid on
% 
% 
% figure(4)
% plot(D1T,R,'linewidth',2)
% title('Nautilus Design Data','fontweight','bold')
% xlabel('D1(m)','fontweight','bold')
% ylabel('Twist Ratio, R','fontweight','bold')
% g=line('XData',[0.1 0.6],'YData',[1.3 1.3],'LineStyle','-.','Linewidth',3,'Color','k');
% h=legend('show','location','best',s1,s2,s3,s4,s5,'Desired R');
% htitle = get(h,'Title');
% set(htitle,'string','D2(m)','fontweight','bold');
% grid on
% 
% figure(5)
% plot(D1T,m,'linewidth',2)
% title ('Nautilus Design Data','fontweight','bold')
% xlabel('D1(m)','fontweight','bold')
% ylabel('Mass of the Shaft (kg)','fontweight','bold')
% h=legend('show','location','best',s1,s2,s3,s4,s5);
% htitle = get(h,'Title');
% set(htitle,'string','D2(m)','fontweight','bold');
% grid on

 %% Calculate Twist Angle along the Shaft
 J1 = (pi/2)*((D1/2)^4);         
J2 = (pi/2)*((D2/2)^4); 
Phi=zeros(1,length(x));
Phi_n=0;
for n=1:80
    if x(n)<L1+L2
        Phi_n=Phi_n+((-Tout*0.1))/(J1*G);
         Phi(n)=Phi_n;
    elseif x(n)<L1+L2+L3+L4+L5
        Phi_n=Phi_n+((-Tout*0.1))/(J2*G);
        Phi(n)=Phi_n;   
    else
        Phi_n=Phi_n+((-Tgear-Tout)*0.1)/(J1*G);
        Phi(n)=Phi_n;  
    end
end
%% Calculate SigmaX

if x<L1+L2
    SigmaX = 32.*(M)/(pi*D1^3);
elseif x<L1+L2+L3+L4+L5
    SigmaX = 32.*(M)/(pi*D2^3);
else
    SigmaX = 32.*(M)/(pi*D1^3);

end

%% Calculate TauXZ
% if x<L1+L2
%     TauXZ = (T*16)/(pi*D1^3);
% elseif x<L1+L2+L3+L4+L5
%     TauXZ = (T*16)/(pi*D2^3);
% else
%     TauXZ = (T*16)/(pi*D1^3);
% end

%% Calculate SigmaM
SigmaM = zeros(1,80);


if x<L1+L2
    SigmaM = Fthrust/(pi*(D1/2)^2); 
elseif x<L1+L2+L3
    SigmaM = Fthrust/(pi*(D2/2)^2); 
else
    SigmaM = 0;
end

%% Calculate SigmaA
% SigmaA = (1/sqrt(2))*(((SigmaZ.^2)+(SigmaZ.^2)+ 6*(TauXZ)).^0.5);
%% Calculate SigmaAR
SigmaAR = SigmaX./(1-(SigmaM./SigmaU));
% figure(2)
% plot(x,SigmaAR)
%% B & A Coefficients 
B = (log((0.9*SigmaU)/(0.5*SigmaU)))/(log(Nf1/Nf2));
A = (0.9*SigmaU)/(Nf1^B);

%% Stress due to expected life
Sigma1 = A*(Nf^B);

%% Factor of Safety
Xs = Sigma1./SigmaAR;
% XSs = Xs(1:50);
% S = x(1:50);
Xsm = min(Xs(1:50));
%% Calculate Shaft Mass
m = dens*pi*((((D1/2)^2)*(L1+L2+L6))+((((D2/2)^2)*(L3+L4+L5))));


%% R Values
PC1 = 0.0175*L1/(20*D1);
PC2 = 0.0175*(L2)/(20*D1);
PC3 = 0.0175*(L3)/(20*D2);
PC4 = 0.0175*(L4)/(20*D2);
PC5 = 0.0175*(L5)/(20*D2);
PC6 = 0.0175*L6/(20*D1);

P1 = max(abs(Phi(1:5)-Phi(1)));
P2 = max(abs(Phi(5:10)-Phi(5)));
P3 = max(abs(Phi(10:20)-Phi(10)));
P4 = max(abs(Phi(20:50)-Phi(20)));
P5 = max(abs(Phi(50:70)-Phi(50)));
P6 = max(abs(Phi(70:80)-Phi(70)));

R1 = PC1/(P1);
R2 = PC2/(P2);
R3 = PC3/(P3);
R4 = PC4/(P4);
R5 = PC5/(P5);
R6 = PC6/(P6);

Rarray = [R1 R2 R3 R4 R5 R6];
R = min(Rarray);

fprintf('Nautilus Design Data \n')
fprintf('Minimum Factor of Safety:  \n'), disp(Xsm)
fprintf('Shaft Twist Ratio: \n'), disp(R)
fprintf('Mass of the Shaft (kg):  \n'), disp(m)

%% Create textbox strings
Pp = ['Pprop =' num2str(Pprop) '[W]'];
Pg = ['Pgen =' num2str(Pgen) '[W]'];
Es = ['E =' num2str(E) '[Pa]'];
Lam = ['\lambda =' num2str(pois) '[-]'];
r = ['\rho =' num2str(dens) '[kg/m^{3}]'];
su = ['\sigma_{u} =' num2str(SigmaU) '[Pa]'];
nf = ['N_{f} =' num2str(Nf) '[Cycles]'];
d1 = ['D_{1} =' num2str(D1) '[m]'];
d2 = ['D_{2} =' num2str(D2) '[m]'];
sm = {Pp Pg Es Lam r su nf d1 d2}';

figure6 = figure('position',[50 50 800 600]);
plot(x,SigmaAR,'linewidth',2)
title ('Nautilus Design Data - \sigma_{ar}','fontweight','bold')
xlabel('Distance along the Shaft (m)','fontweight','bold')
ylabel('Equivalent Alternating Stress, \sigma_{ar} (N/m^{2})','fontweight','bold')
h = annotation('textbox',[0.6 0.7 0.3 0.15],'string',sm,'fitboxtotext','on','BackgroundColor',[0.9  0.9 0.9]);
grid on

figure7 = figure('position',[900 50 800 600]);
plot(x(1:50),Xs(1:50),'linewidth',2)
xlim([0 8])
ylim([0 600])
title ('Nautilus Design Data - X_{s}','fontweight','bold')
xlabel('Distance along the Shaft (m)','fontweight','bold')
ylabel('Factor of Safety, X_{s}','fontweight','bold')
h = annotation('textbox',[0.6 0.7 0.3 0.15],'string',sm,'fitboxtotext','on','BackgroundColor',[0.9  0.9 0.9]);
grid on
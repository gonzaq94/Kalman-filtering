% Ejercicio 4: se tiene un sesgo en posición, por lo que se expande el
% sistema para seguir utilizando el Kalman lineal. Se ve que este sesgo en
% posición no puede detectarse pues el sistema es no observable.
close all;

GyroStruct=load('Gyro.mat');
RadarStruct=load('Radar-sesgo-pos.mat');
AcelStruct=load('Acel.mat');
TrayectoriaStruct=load('trayectoria.mat');

Gyro=table2array(struct2table(GyroStruct));
Radar=table2array(struct2table(RadarStruct));
Acel=table2array(struct2table(AcelStruct));
Trayectoria=table2array(struct2table(TrayectoriaStruct));

X=zeros(8,length(Trayectoria));
S=zeros(2,length(Trayectoria));
Xexpandido=[X;S];

% Condiciones iniciales
deltaP=100;
deltaV=0.2;
deltaT=40;

Pinicial=Trayectoria(1,2:3)+deltaP/4*randn(1);
Vinicial=Trayectoria(1,5:6)+deltaV/4*randn(1);
Thetainicial=Trayectoria(1,8)+deltaT/4*randn(1);

C11_0=cos(Thetainicial);
C12_0=-sin(Thetainicial);
C21_0=-C12_0;
C22_0=C11_0;
C_0=[C11_0,C12_0;
    C21_0,C22_0;];

deltaC11=abs(sin(Thetainicial)*deltaT);
deltaC12=abs(cos(Thetainicial)*deltaT);
deltaC21=deltaC12;
deltaC22=deltaC11;

X0_0=[Pinicial';Vinicial';C11_0;C12_0;C21_0;C22_0];
S0_0=zeros(2,1);
Xexpandido0_0=[X0_0;S0_0];

deltas=[deltaP;deltaP;deltaV;deltaV;deltaC11;deltaC12;deltaC21;deltaC22;deltaP;deltaP];
Pexpandido0_0=deltas*deltas';
%Pexpandido0_0=[deltas;deltaP;deltaP]*[deltas;deltaP;deltaP]';

Xexpandido(:,1)=Xexpandido0_0;
Pkexpandido=Pexpandido0_0;

%Matrices del sistema
deltaPruido=100;
deltaVruido=0.1;
deltaCruido=10^(-10);

I2=eye(2,2);
z=zeros(2,2); 

Rkexpandido= diag([deltaPruido^2 deltaPruido^2 deltaVruido^2 deltaVruido^2 deltaCruido^2 deltaCruido^2]);

Ck = [1 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0;
      0 0 0 1 0 0 0 0;
      0 0 0 0 1 0 0 -1;
      0 0 0 0 0 1 1 0];
Ckexpandido=[Ck,eye(6,2)];

Bk=diag(ones(1,8)/100);
Bkexpandido=[Bk;zeros(2,6),eye(2,2)/10000];

Qkexpandido=diag(ones(1,8));
%Qkexpandido=[Qk,zeros(8,2)];

Tk=1/100;

j=1;

for k=1:length(Trayectoria)-1
    
    if (Radar(j,1)==Acel(k,1))
        yk=[RadarStruct.Pradar(j,2);RadarStruct.Pradar(j,3);RadarStruct.Vradar(j,2);RadarStruct.Vradar(j,3);0;0];
        Kkexpandido=Pkexpandido*Ckexpandido'/(Rkexpandido+Ckexpandido*Pkexpandido*Ckexpandido');
        Xexpandido(:,k)=Xexpandido(:,k)+Kkexpandido*(yk-Ckexpandido*Xexpandido(:,k));
        Pkexpandido=(eye(10,10)-Kkexpandido*Ckexpandido)*Pkexpandido;
        j=j+1;
    end
    
ax_b=Acel(k,2);
ay_b=Acel(k,3);
w_b=Gyro(k,2);

M_a=[ax_b,ay_b,0,0;
    0,0,ax_b,ay_b;];
M_w=[0,w_b,0,0;
    -w_b,0,0,0;
    0,0,0,w_b;
    0,0,-w_b,0;];

Ak=[I2,Tk*I2,(Tk^2/2)*M_a;
    zeros(2,2),I2,Tk*M_a;
    zeros(4,4),M_w];
Akexpandido=[Ak,zeros(8,2); zeros(2,8), eye(2,2)];

Xexpandido(:,k+1)=Akexpandido*Xexpandido(:,k);
Pkexpandido=Akexpandido*Pkexpandido*Akexpandido'+Bkexpandido*Qkexpandido*Bkexpandido';

end

sesgo=Xexpandido(9:10,:);

figure (1)
plot(Trayectoria(:,2),Trayectoria(:,3));
hold on;
plot(Xexpandido(1,:),Xexpandido(2,:));
xlabel('Eje x(metros)');
ylabel('Eje y(metros)');
title('Trayectorias con sesgo en posición');
grid on;
legend('Trayectoria real','Trayectoria estimada');
print('trayectorias-sesgo-pos.png','-dpng');

error=sqrt(Trayectoria(:,2).^2+Trayectoria(:,3).^2)-sqrt((Xexpandido(1,:))'.^2+(Xexpandido(2,:))'.^2);
figure (4)
plot(Trayectoria(:,1),error);
grid on;
title('Error con sesgo en posición');
xlabel('Tiempo(segundos)');
ylabel('Error (metros)');
print('error-sesgo-pos.png','-dpng');
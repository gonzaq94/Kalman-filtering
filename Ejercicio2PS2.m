% Ejercicio 2: se implementa el filtro Kalman lineal. 
close all;

GyroStruct=load('Gyro.mat');
RadarStruct=load('Radar.mat');
AcelStruct=load('Acel.mat');
TrayectoriaStruct=load('trayectoria.mat');

Gyro=table2array(struct2table(GyroStruct));
Radar=table2array(struct2table(RadarStruct));
Acel=table2array(struct2table(AcelStruct));
Trayectoria=table2array(struct2table(TrayectoriaStruct));

X=zeros(8,length(Trayectoria));

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

deltas=[deltaP;deltaP;deltaV;deltaV;deltaC11;deltaC12;deltaC21;deltaC22];
P0_0=deltas*deltas';

X(:,1)=X0_0;
Pk=P0_0;

%Matrices del sistema
deltaPruido=100;
deltaVruido=0.1;
deltaCruido=10^(-10);

I2=eye(2,2);
z=zeros(2,2); 

Rk = diag([deltaPruido^2 deltaPruido^2 deltaVruido^2 deltaVruido^2 deltaCruido^2 deltaCruido^2]);

Ck = [1 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0;
      0 0 0 1 0 0 0 0;
      0 0 0 0 1 0 0 -1;
      0 0 0 0 0 1 1 0];

Bk=diag(ones(1,8));

Qk=diag(ones(1,8));

Tk=1/100;

j=1;

for k=1:length(Trayectoria)-1
    
    if (Radar(j,1)==Acel(k,1))
        yk=[RadarStruct.Pradar(j,2);RadarStruct.Pradar(j,3);RadarStruct.Vradar(j,2);RadarStruct.Vradar(j,3);0;0];
        Kk=Pk*Ck'/(Rk+Ck*Pk*Ck');
        X(:,k)=X(:,k)+Kk*(yk-Ck*X(:,k));
        Pk=(eye(8,8)-Kk*Ck)*Pk;
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

X(:,k+1)=Ak*X(:,k);
Pk=Ak*Pk*Ak'+Bk*Qk*Bk';

end
figure (1)
plot(Trayectoria(:,2),Trayectoria(:,3));
hold on;
plot(X(1,:),X(2,:));
xlabel('Eje x(metros)');
ylabel('Eje y(metros)');
grid on;
title('Trayectorias real y estimada');
legend('Trayectoria real','Trayectoria estimada');
print('trayectorias.png','-dpng');

error=sqrt(Trayectoria(:,2).^2+Trayectoria(:,3).^2)-sqrt(X(1,:)'.^2+X(2,:)'.^2);
figure (4)
plot(Trayectoria(:,1),error);
grid on;
title('Error de estimación');
xlabel('Tiempo(segundos)');
ylabel('Error (metros)');
print('error.png','-dpng');

% Ej. 7)Si no se dispone de los datos de la trayectoria real, la manera de
% verificar que lo que se esta haciendo esta bien es calcular la función de
% autocorrelacion de las innovaciones gk y verificar que el resultado sea
% una delta discreta, lo que indicaria que se verifica que el principio de
% ortogonalidad. Otra manera sería graficar la evolución de la matriz P y
% verificar que la misma se achique con el tiempol, indicando que el error
% de estimación disminuye.
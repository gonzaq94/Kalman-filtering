% Ejercicio 6: Como se tiene un sesgo en aceleración, se debe recurrir al
% filtro de Kalman extendido

close all;
clear all;

GyroStruct=load('Gyro.mat');
RadarStruct=load('Radar.mat');
AcelStruct=load('Acel.mat');
TrayectoriaStruct=load('trayectoria.mat');

Gyro=table2array(struct2table(GyroStruct));
Radar=table2array(struct2table(RadarStruct));
Acel=table2array(struct2table(AcelStruct));
Trayectoria=table2array(struct2table(TrayectoriaStruct));

X=zeros(7,length(Trayectoria));
sesgox=0.2;
sesgoy=0.5;

% Condiciones iniciales
deltaP=100;
deltaV=0.2;
deltaT=40;
deltaB=0.2;

Pinicial=Trayectoria(1,2:3)+deltaP/4*randn(1);
Vinicial=Trayectoria(1,5:6)+deltaV/4*randn(1);
Thetainicial=Trayectoria(1,8)+deltaT/4*randn(1);
Binicial=zeros(2,1);

X0_0=[Pinicial';Vinicial';Thetainicial;Binicial];

P0_0 = [eye(2,2)*deltaP,zeros(2,5);
      zeros(2,2),eye(2,2)*deltaV,zeros(2,3);
      zeros(1,4),deltaT*pi/180,zeros(1,2);
      zeros(2,5),eye(2,2)*deltaB];

X(:,1)=X0_0;
Pk=P0_0;

%Matrices del sistema
deltaPruido=100;
deltaVruido=0.1;

I2=eye(2,2);
z=zeros(2,2); 

Rk = diag([deltaPruido^2 deltaPruido^2 deltaVruido^2 deltaVruido^2]);

Ck = [ eye(4,4),zeros(4,3)];

Bk = [eye(5),zeros(5,2);
          zeros(2,7)];

Qk = diag([0.001,0.001,0.01,0.01,0.1,0.1,0.1]);

Tk=1/100;

j=1;

for k=2:length(Trayectoria)
    
    ax=Acel(k-1,2)+sesgox; 
    ay=Acel(k-1,3)+sesgoy;
    w=Gyro(k-1,2);
    
    posx=X(1,k-1);
    posy=X(2,k-1);
    vx=X(3,k-1);
    vy=X(4,k-1);
    theta=X(5,k-1);
    bx=X(6,k-1);
    by=X(7,k-1);
    
    uk=[ posx+Tk*vx+Tk^2/2*ax ;
         posy+Tk*vy+Tk^2/2*ay ;
         Tk*cos(theta)*(ax-bx)-Tk*sin(theta)*(ay-by)+vx ;
         Tk*sin(theta)*(ax-bx)+Tk*cos(theta)*(ay-by)+vy ;
         Tk*w+theta ;
         bx ;
         by ;];
     
    Maux = [-sin(theta)*(ax-bx)-cos(theta)*(ay-by),-cos(theta),sin(theta);
            cos(theta)*(ax-bx)-sin(theta)*(ay-by),-sin(theta),-cos(theta)];
     
    Ak = [eye(2,2),Tk*eye(2,2),(Tk^2/2)*Maux;
          zeros(2,2),eye(2,2),Tk*Maux;
          zeros(3,4),eye(3,3)];                  
      
    X(:,k)=uk;  
    Pk=Ak*Pk*Ak'+Bk*Qk*Bk';
    
    if (Radar(j,1)==Acel(k-1,1))
      
        yk=[RadarStruct.Pradar(j,2);RadarStruct.Pradar(j,3);RadarStruct.Vradar(j,2);RadarStruct.Vradar(j,3);];
        Kk=Pk*Ck'/(Rk+Ck*Pk*Ck');
        X(:,k)=X(:,k)+Kk*(yk-Ck*X(:,k));
        Pk=(eye(7,7)-Kk*Ck)*Pk;
        j=j+1;
    else
    end

end
figure (1)
plot(Trayectoria(:,2),Trayectoria(:,3));
hold on;
plot(X(1,:),X(2,:));
xlabel('Eje x(metros)');
ylabel('Eje y(metros)');
grid on;
title('Trayectorias con sesgo en acelerómetro');
legend('Trayectoria real','Trayectoria estimada');
print('trayectoriasExt.png','-dpng');

error=sqrt(Trayectoria(:,2).^2+Trayectoria(:,3).^2)-sqrt(X(1,:)'.^2+X(2,:)'.^2);
figure (4)
plot(Trayectoria(:,1),error);
grid on;
title('Error con sesgo en acelerómetro');
xlabel('Tiempo(segundos)');
ylabel('Error (metros)');
print('errorExt.png','-dpng');

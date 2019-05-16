% Ejercicio 5: se implementan las factorizaciones QR y de Cholesky para una
% implementación más óptima del filtro de Kalman.
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

X(:,1)=[Pinicial';Vinicial';C11_0;C12_0;C21_0;C22_0];;
Pk=diag([deltaP,deltaP,deltaV,deltaV,deltaC11,deltaC12,deltaC21,deltaC22]);

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
Zk=sqrtm(Pk);
w2=X(:,1)'*inv(Zk');
Akacumulada=eye(8,8);

for k=1:length(Trayectoria)-1    
    
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
    Akacumulada=Ak*Akacumulada;
    
    if (Radar(j,1)==Acel(k,1))
        
        Mk=[sqrtm(Rk), Ck*Zk,zeros(6,8);
            zeros(8,6), Akacumulada*Zk, Bk*sqrtm(Qk*100)*sqrt(Tk)]; 
        %el factor multiplicativo "100" se puso arbitrariamente para que de
        %mas lindo
        [Q,R]=qr(Mk');
        signos=sign(diag(R));
        
        for i=1:length(signos)
            R(i,:)=signos(i)*R(i,:);
            Q(:,i)=signos(i)*Q(:,i);
        end
        
        Rtransp=R';
        Z=Rtransp(7:14,7:14);
        yk=[RadarStruct.Pradar(j,2);RadarStruct.Pradar(j,3);RadarStruct.Vradar(j,2);RadarStruct.Vradar(j,3);0;0];
        VecAuxiliar=[-yk'*inv(sqrtm(Rk)'), w2, zeros(1,8)]*Q;
        w2=VecAuxiliar(1,7:14);
        X(:,k+1)=Z*w2';
        Pk=Z'*Z;
        j=j+1;
        Akacumulada=eye(8,8);
    else
        X(:,k+1)=Ak*X(:,k);
        Pk=Ak*Pk*Ak'+Bk*Qk*Bk';
    end
end

figure (1)
plot(Trayectoria(:,2),Trayectoria(:,3));
hold on;
plot(X(1,:),X(2,:));
xlabel('Eje x(metros)');
ylabel('Eje y(metros)');
grid on;
title('Trayectorias  con implementación óptima');
legend('Trayectoria real','Trayectoria estimada');
print('trayectoriasQR.png','-dpng');

error=sqrt(Trayectoria(:,2).^2+Trayectoria(:,3).^2)-sqrt(X(1,:)'.^2+X(2,:)'.^2);
figure (4)
plot(Trayectoria(:,1),error);
grid on;
title('Error con implementación óptima');
xlabel('Tiempo(segundos)');
ylabel('Error (metros)');
print('errorQR.png','-dpng');
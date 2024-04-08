clear all
%close all
%% Add functions
 addpath('../');


%% SDOF system

m=1;
k=1;
ep = 0.02;
xi_na =0.3;
xi = 0.5;
xi_gd = 0*1.14;


omega0 = sqrt(k/m);
P = 0.61;

kappa= 0;


% Tuning parameters, these were changed for the other simulations

%% Nonlinear absorber
% Simulation parameters

    Tl =2000;

    t=0:0.1:Tl;
    N = length(t);
    f = 0:10/N:(N-1)*10/N;
    
    w = 0.98:0.001:1.02;
for i=1:length(w)
    x0_real= 0;

    x0 = [0 sqrt(-kappa)];
    x0_dot = [0 0];

 f1 = @(t,y)[y(3);y(4);...
    -(y(1) + ep*xi*y(3)- ep*kappa*(y(2)-y(1)) - ep*(y(2)-y(1))^3 -ep*xi_na*(y(4)-y(3))-ep*P*cos(w(i)*t))/m;...
    -(ep*kappa*(y(2)-y(1)) + ep*(y(2)-y(1))^3 + ep*xi_na*(y(4)-y(3)))/ep;...
    ];

    Prec = 1e-8;
    % Actual numerical simulation

    options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
          [T1,Y1] = ode45(f1,t,[x0 x0_dot],options);

    % Plot simulation result      

    Y1filt = bandpass(Y1,[0.9*w(i)/(2*pi) 1.1*w(i)/(2*pi)],10,'ImpulseResponse','iir','Steepness',0.999999);
%     A_main(i)=mean(envY(round(3*length(Y1)/4):end));
% 
%     B_main(i)=mean(envB(round(3*length(Y1)/4):end));
    envY  = envelope(Y1(:,1),150,'peaks'); %sqrt(Y1(:,1).^2+Y1(:,3).^2);
envB  = envelope(abs(Y1(:,1)-Y1(:,2)),150,'peaks'); %sqrt((Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2);

     A_main(i)=rms(Y1(round(3*length(Y1)/4):end,1))*sqrt(2);
    B_main(i)=rms(Y1(round(3*length(Y1)/4):end,2)-Y1(round(3*length(Y1)/4):end,1))*sqrt(2);
    
    TF = islocalmin(envY,'MinSeparation',100,'SamplePoints',T1);
    if(length(TF)>0)
            envY_min = envY(TF)
     A_main_min(i)= min(envY_min(round(3*length(envY_min)/4):end))
     if(A_main_min(i)/ A_main(i) > 0.9)
         A_main_min(i)=NaN;
     end
                   % A_main_min(i)=min(envY(round(3*length(Y1)/4):end));
                A_main_max(i)=max(envY(round(3*length(Y1)/4):end));
                    B_main_min(i)=min(envB(round(3*length(Y1)/4):end));
                B_main_max(i)=max(envB(round(3*length(Y1)/4):end));
                
                 if(A_main_max(i)/ A_main(i) < 1.1)
         A_main_max(i)=NaN;
                 end
     
    else
     A_main_min(i)= NaN;
                   % A_main_min(i)=min(envY(round(3*length(Y1)/4):end));
                A_main_max(i)= NaN;
                    B_main_min(i)= NaN;
                B_main_max(i)= NaN;
    end
    



% % 
figure(100)
plot(T1,Y1(:,1))
hold on
plot(T1,envY,T1(TF),envY(TF),'r*')
% 
xlim([500, 2000])
     pause(0.8)
 clf
figure(101)
plot(T1,Y1(:,2)-Y1(:,1))
hold on
plot(T1,envB)
xlim([500, 2000])
title(num2str(w(i)))
     pause(0.2)
clf

    i
end
%%
figure(1)
% plot(w,A_main/(ep*m*omega0^2*F)*sqrt(2),'-o','LineStyle', 'none' )
 hold on
plot(w,A_main)%,'-o','LineStyle', 'none' )
 hold on
plot(w,A_main_min,'-o','LineStyle', 'none' )
 hold on
plot(w,A_main_max,'-o','LineStyle', 'none' )

figure(2)
% hold on
% plot(w,B_main*sqrt(2),'-o','LineStyle', 'none' )
hold on
plot(w,B_main,'-o','LineStyle', 'none' )
 hold on
plot(w,B_main_min,'-o','LineStyle', 'none' )
 hold on
plot(w,B_main_max,'-o','LineStyle', 'none' )



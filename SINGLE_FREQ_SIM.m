clear all
%close all
%% Add functions
 addpath('../');


%% SDOF system

m=1;
k=1;
ep = 0.02;
xi_na =0.2;
xi = 0.5;

omega0 = sqrt(m/k);
P = 0.61;

kappa= 0;

% Tuning parameters, these were changed for the other simulations

%% Nonlinear absorber
% Simulation parameters
    Tl =2000;

    t=0:0.05:Tl;
    N = length(t);
    f = 0:10/N:(N-1)*10/N;
    
    w = 1;

    x0_real= 0;

    x0 = [0 sqrt(-kappa)];
    x0_dot = [0 0];

 f1 = @(t,y)[y(3);y(4);...
    -(y(1) + ep*xi*y(3)- ep*kappa*(y(2)-y(1)) - ep*(y(2)-y(1))^3 -ep*xi_na*(y(4)-y(3)) -ep*P*cos(w*t))/m;...
    -(ep*kappa*(y(2)-y(1)) + ep*(y(2)-y(1))^3 + ep*xi_na*(y(4)-y(3)))/ep;...
    ];

    Prec = 1e-8;
    % Actual numerical simulation

    options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
          [T1,Y1] = ode45(f1,t,[x0 x0_dot],options);
    % Plot simulation result      

    Y1filt = bandpass(Y1,[0.90/(2*pi) 1.1/(2*pi)],10,'ImpulseResponse','iir','Steepness',0.999999);
envY  = envelope(Y1(:,1),200,'peaks'); %sqrt(Y1(:,1).^2+Y1(:,3).^2);
envB  = envelope(abs(Y1(:,1)-Y1(:,2)),200,'peaks'); %sqrt((Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2);

figure
box on

subplot(2,1,1);
plot(T1,Y1(:,1),'k','LineWidth',2)
hold on
plot(T1,envY,'LineWidth',2)
%xlim([1000, 1300])
   ax = gca; 
ax.FontSize = 15; 
ylabel('$A$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
 
subplot(2,1,2);
box on

plot(T1,Y1(:,2)-Y1(:,1),'k','LineWidth',2)
hold on
plot(T1,envB,'LineWidth',2)
%xlim([1000, 1300])
   ax = gca; 
ax.FontSize = 15; 
ylabel('$B$','FontSize',16.5,'interpreter','latex') 
xlabel('$\tau$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis

     Zb = 0:0.001:2.5;  
     Za = (Zb.*(w^2*xi_na^2+(w^2-kappa-3/4*Zb).^2))/(w^4);

     phase = asin(xi_na*w*sqrt(Zb)./(w^2*sqrt(Za)));
     
     for i=2:length(Zb)
            m11= stability_m11( sqrt(Zb(i))*exp(1i*phase(i)),xi_na,kappa,w^2,xi);
            m22 = conj(m11);
            m12 = stability_m12( sqrt(Zb(i))*exp(1i*phase(i)),xi_na,kappa,w^2,xi);
            m21 = conj(m12);
            M = [m11, m12;
               m21, m22];
         if(any(real(eig(M))>0))
               B_unstable(i) =sqrt(Zb(i)); 
               A_unstable(i) =sqrt(Za(i)); 
               B_stable(i) = NaN; 
               A_stable(i) = NaN;
            else
           B_stable(i) =sqrt(Zb(i)); 
               A_stable(i) =sqrt(Za(i));  
               B_unstable(i) = NaN; 
               A_unstable(i) = NaN;
            end
     end
     

%%
figure(501)
box on
hold on
plot(B_stable,A_stable,'k','LineWidth',2 )
plot(B_unstable,A_unstable,'k--','LineWidth',2 )
% plot(w,A_main/(ep*m*omega0^2*F)*sqrt(2),'-o','LineStyle', 'none' )
plot(envB,envY,'r' ,'LineWidth',2 )

   ax = gca; 
ax.FontSize = 15; 
xlabel('$B$','FontSize',16.5,'interpreter','latex') 
ylabel('$A$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
%legend('$\Omega=1$','$\Omega=0.965$','$\Omega=0.965$ $x(0)''=0.5$','Interpreter','latex')


function m11  = stability_m11(B,xi_na,kappa,X,xi)
        A = conj(B);
           m11  = -3/2*A.*B-kappa;
          m11 = (-1i*X - xi_na*sqrt(X) -1i*m11);

       % m11 = (-1i*X - xi_na*sqrt(X) -1i*m11-1i*X*X/(sigma+1i*xi*sqrt(X)-X))/( sqrt(X)*2*(sigma+1i*xi*sqrt(X)) )*(sigma+1i*xi*sqrt(X)-X);
end

function m12  = stability_m12(B,xi_na,kappa,X,xi)
    A = conj(B);
       m12  = -3/4*B^2;
       m12 = -1i*m12
   % m12 = -1i*m12/( sqrt(X)*2*(sigma+1i*xi*sqrt(X)) )*(sigma+1i*xi*sqrt(X)-X);
end


%==========================================================
clearvars;
%close all;
clc;
%% Define system
m=1;
k=1;
ep = 0.02;
xi_na = 0.2;
xi = 0.2;

m=1;
omega0 = sqrt(k/m);
P = 0.61;
kappa=0;
Om=linspace(0.95,1.05,10000);

for i=1:length(Om)
        X = Om(i)^2;
    sigma = (1-X)/ep;
    a =  ((X-sigma)*(3/4))^2 + X*(xi*3/4)^2;
    b = 2*((X-sigma)*3/4)*((X-sigma)*(-X+kappa)+X*xi_na*xi+X^2)+2*X*(xi*3/4)*(xi*(-X+kappa)+(sigma-X)*xi_na);
    c = ((X-sigma)*(-X+kappa)+X*xi_na*xi+X^2)^2+X*(xi*(-X+kappa)+(sigma-X)*xi_na)^2;
    d = -X^2*P^2;
    
    
    %%TODO: Find minimum of the SIM and the Melkinov!!
    [r] = roots([a,b,c,d]);

    if(~isreal(r))
         r_real=r(imag(r)==0);
          zo(i,1) = ( X*xi_na^2+ ( X-kappa-3/4*r_real)^2)*r_real/(X^2);

        zo(i,2)=NaN;
        zo(i,3)=NaN;    
        zna(i,1) = r_real;
        zna(i,2)=NaN;
        zna(i,3)=NaN;
        
                            alpha1(i,1) = asin( xi*sqrt(X)*sqrt(zo(i,1))/P); %- asin( (xi*sqrt(X)*sqrt(zo(i,1))+xi_na*X*sqrt(zna(i,1))/(sqrt(X)*sqrt(zo(i,1))))/P );

                    beta1(i,1) = asin(xi_na* sqrt(zna(i,1)/ sqrt(zo(i,1))/sqrt(X)))-   alpha1(i,1);

           
           
           M2= [-ep*sigma - 1i*xi*ep*sqrt(X), 0, 1i*ep*xi_na*sqrt(X)+6/4*ep*zna(i,1)+ep*kappa, 3/4*ep*zna(i,1)*exp(1i*beta1(i,1))^2;
            0,  ep*sigma - 1i*xi*ep*sqrt(X), -conj(3/4*ep*zna(i,1)*exp(1i*beta1(i,1))^2),  1i*ep*xi_na*sqrt(X)-6/4*ep*zna(i,1)-ep*kappa;
            ep*sigma + 1i*xi*ep*sqrt(X)+ X , 0, -1i*xi_na*sqrt(X)*(1+ep)+X-(1+ep)*(6/4*zna(i,1)+kappa),  -(1+ep)*3/4*zna(i,1)*exp(1i*beta1(i,1))^2;
            0,   -ep*sigma + 1i*xi*ep*sqrt(X) - X,  -conj( -(1+ep)*3/4*zna(i,1)*exp(1i*beta1(i,1))^2), -conj(-1i*xi_na*sqrt(X)*(1+ep)+X-(1+ep)*(6/4*zna(i,1)+kappa));]/(sqrt(X)*2*1i );
        eig(M2);
            if(any(real(eig(M2))>0))
               B_unstable(i,1) =sqrt(zna(i,1)); 
               A_unstable(i,1) =sqrt(zo(i,1)); 
               B_stable(i,1) = NaN; 
               A_stable(i,1) = NaN;
            else
           B_stable(i,1) =sqrt(zna(i,1)); 
               A_stable(i,1) =sqrt(zo(i,1));  
               B_unstable(i,1) = NaN; 
               A_unstable(i,1) = NaN;
            end
                 B_stable(i,2) = NaN; 
               A_stable(i,2) = NaN;
                    B_stable(i,3) = NaN; 
               A_stable(i,3) = NaN;
            B_unstable(i,2) = NaN; 
               A_unstable(i,2) = NaN;
               B_unstable(i,3) = NaN; 
               A_unstable(i,3) = NaN;
    else
     zna(i,1) = r(1);
        zna(i,2)= r(2);
        zna(i,3)=r(3);
        zna(i,:) =sort(zna(i,:),'ascend');
                r =sort(r,'ascend');

        zo(i,1)= ( X*(xi_na)^2+ ( X-kappa-3/4*r(1))^2)*r(1)/(X^2);
        zo(i,2)= ( X*(xi_na)^2+ ( X-kappa-3/4*r(2))^2)*r(2)/(X^2);
        zo(i,3)= ( X*(xi_na)^2+ ( X-kappa-3/4*r(3))^2)*r(3)/(X^2);
        
       
       % zo(i,:) =sort(zo(i,:),'ascend') ;

        
        

        
        
        for j=1:3
                                        alpha1(i,j) =  asin( xi*sqrt(X)*sqrt(zo(i,j))/P); %-asin( (xi*sqrt(X)*sqrt(zo(i,j))+xi_na*X*sqrt(zna(i,j))/(sqrt(X)*sqrt(zo(i,j))))/P );

                    beta1(i,j) = asin(xi_na* sqrt(zna(i,j)/ sqrt(zo(i,j))/sqrt(X)))-   alpha1(i,j);

%            
           
             M2= [-ep*sigma - 1i*xi*ep*sqrt(X), 0, 1i*ep*xi_na*sqrt(X)+6/4*ep*zna(i,j)+ep*kappa, 3/4*ep*zna(i,j)*exp(1i*beta1(i,j))^2;
            0,  ep*sigma - 1i*xi*ep*sqrt(X), -conj(3/4*ep*zna(i,j)*exp(1i*beta1(i,j))^2),  1i*ep*xi_na*sqrt(X)-6/4*ep*zna(i,j)-ep*kappa;
            ep*sigma + 1i*xi*ep*sqrt(X)+ X , 0, -1i*xi_na*sqrt(X)*(1+ep)+X-(1+ep)*(6/4*zna(i,j)+kappa),  -(1+ep)*3/4*zna(i,j)*exp(1i*beta1(i,j))^2;
            0,   -ep*sigma + 1i*xi*ep*sqrt(X) - X,  -conj( -(1+ep)*3/4*zna(i,j)*exp(1i*beta1(i,j))^2), -conj(-1i*xi_na*sqrt(X)*(1+ep)+X-(1+ep)*(6/4*zna(i,j)+kappa));]/(sqrt(X)*2*1i );
        eig(M2);
        
            if(any(real(eig(M2))>0))
               B_unstable(i,j) =sqrt(zna(i,j)); 
               A_unstable(i,j) =sqrt(zo(i,j)); 
               B_stable(i,j) = NaN; 
               A_stable(i,j) = NaN;
            else
           B_stable(i,j) =sqrt(zna(i,j)); 
               A_stable(i,j) =sqrt(zo(i,j));  
               B_unstable(i,j) = NaN; 
               A_unstable(i,j) = NaN;
            end
        end
        

    end
    
    % Werk nog niet voor nonlinear geometric damping
        a2 = 27/16;
    b2= -3*(X-kappa);
    c2 = ((X-kappa)^2+X*xi_na^2);
    [r] = roots([a2,b2,c2]);
    za_mel = -2*kappa;
 if(isreal(r))
         r_min=min(r);
         zna_min(i) = r_min;
          zo_max(i) =  ( X*(xi_na)^2+ ( X-kappa-3/4*r_min)^2)*r_min/(X^2);
         r_max=max(r);
                  zna_max(i) = r_max;
          zo_min(i) = ( X*(xi_na)^2+ ( X-kappa-3/4*r_max)^2)*r_max/(X^2);
 end
    
    
end

    zeta = xi/ep*2;
    
figure(1)
hold on
plot(Om,A_stable,'k','LineWidth',2)
hold on
plot(Om,A_unstable,'k--','LineWidth',2)
plot(Om,sqrt(zo_max),'k--')
plot(Om,sqrt(zo_min),'k--')
hold on
plot(Om,P*ep./sqrt( (1-Om.^2).^2+4*zeta^2*Om.^2),'k:')

   ax = gca; 
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
xlim([min(Om) max(Om)])
figure(2)
hold on
hold on
plot(Om,B_stable,'k','LineWidth',2)
hold on
plot(Om,B_unstable,'k--','LineWidth',2)
plot(Om,sqrt(zna_max),'k--')
plot(Om,sqrt(zna_min),'k--')


   ax = gca; 
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$b$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
xlim([min(Om) max(Om)])


function m11  = stability_m11(B,xi_na,kappa,X,xi,sigma)
        A = conj(B);
           m11  = -3/2*A.*B-kappa;
          m11 = (-1i*X - xi_na*sqrt(X) -1i*m11);

       % m11 = (-1i*X - xi_na*sqrt(X) -1i*m11-1i*X*X/(sigma+1i*xi*sqrt(X)-X))/( sqrt(X)*2*(sigma+1i*xi*sqrt(X)) )*(sigma+1i*xi*sqrt(X)-X);
end

function m12  = stability_m12(B,xi_na,kappa,X,xi,sigma)
    A = conj(B);
       m12  = -3/4*B^2;
       m12 = -1i*m12
   % m12 = -1i*m12/( sqrt(X)*2*(sigma+1i*xi*sqrt(X)) )*(sigma+1i*xi*sqrt(X)-X);
end




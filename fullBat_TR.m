%Code Begins Finally

%defining Computation parameters
dt=1e-9;
dx=1E-6;    %INCREASE OR DECREASE BY A FACTOR OF 10, ALWAYS
t_end=10;
t=0-dt;

%defining geometric parameters
%the lengths of each region
l1=9E-6;    %al current collector 
l2=59E-6;   %positive active material 
l3=20E-6;   %separator
l4=92E-6;   %negative active material 
l5=16E-6;   %Cu current collector 
L=l1+l2+l3+l4+l5; 
Larr=[l1 l2 l3 l4 l5];
%defining thermophysical properties of each material  
%structure P_i=[rho  cp  k] 
P_cu=[8900; 385; 398];        %cu
P_gr=[2660; 1437.4; 1.04];    %graphite
P_ad=[1750; 1120; 0.12];      %additives
P_sep=[492; 1978; 0.334];     
P_al=[1500; 903; 238];
P_ncm=[2380; 710; 1.5];
P_elec=[1290; 133.9; 0.45];

%PROPERTIES OF REGION
%for this weights are used which depend on volumetric fraction of that
%material in the region
w1= 1;                      %only Al
w2=[0.277; 0.523; 0.2];     %[elec ncm ad]
w3=[0.4; 0.6];              %[elec sep]
w4= [0.301; 0.523; 0.176];  %[elec gr ad] 
w5=1;                       %only Cu
%computing weighed properties 
p1=P_al*w1;
p2=[P_elec P_ncm P_ad]*w2;
p3=[P_elec P_sep]*w3;
p4=[P_elec P_gr P_ad]*w4;
p5=P_cu*w5;

n=(L/dx)-4; %number of nodes
%initialiseing temperature matrix
Ti=zeros(1,n);        %spatial temperature array in previous temperature step 
Tf=zeros(1,n);        %updated spatial temperature array

T_initial=80+273;   %Initial constant battery temperature
Tamb=25+273;        %ambiant temperature

Tf=T_initial*ones(n,1);
dT=0;
dTmax=0;    
dTlim=1;    %maximum allowed delta T per timestep

%rection initialisation 
yo=[0.999;0.999;1;0.15;1;1];%
y=yo;
for i=1:1:n-1
    y=[y yo];%each column of y will store the updated reaction progression
end
ynew=y; %will be used to store the new conc in every time step

h=100; %HTC

while t<t_end
disp(['t ' num2str(t) ' dt ' num2str(dt)]) ;
t=t+dt;
Ti=Tf;
y=ynew;
for j=1:1:n
    if (j==1)
        rho=p1(1);
        cp=p1(2);
        k=p1(3);
        dT=((k*(Ti(j+1)-Ti(j))/dx)-h*(Ti(j)-Tamb))*(dt/(rho*cp*dx));
        if (dT>dTmax)
            dTmax=dT;
        end
        Tf(j)=Ti(j)+dT;
    elseif (j==n)
        rho=p5(1);
        cp=p5(2);
        k=p5(3);
        dT=((k*(Ti(j-1)-Ti(j))/dx)-h*(Ti(j)-Tamb))*(dt/(rho*cp*dx));
        if (dT>dTmax)
            dTmax=dT;
        end
        Tf(j)=Ti(j)+dT;
    else
        [rho, cp, k, reg , inter]=prop(j, dx, Larr, p1, p2, p3, p4, p5); %reg is region number (1 to 5) and 0 for interface, inter is interface number (1 to 4) and 0 for internal region
        %code for computing reaction heat and updating conc
        [ynew(:,j),Qreac]=reactions(y(:,j),Ti(j),dt, reg, inter, p2, p3, p4);
        if (inter==0)
            dT=(k*dt/((dx^2)*rho*cp))*(Ti(j+1)-2*Ti(j)+Ti(j-1))+Qreac*dt/(rho*cp);        
        elseif (inter==1)
            rhoCp=p1(1)*p1(2)+p2(1)*p2(2);
            dT=(dt/(1.5*dx*dx*rhoCp))*(-p1(3)*(Ti(j)-Ti(j-1))+p2(3)*(Ti(j+1)-Ti(j)))+Qreac*dt/rhoCp;
        elseif (inter==2)
            rhoCp=p2(1)*p2(2)+p3(1)*p3(2);
            dT=(dt/(1.5*dx*dx*rhoCp))*(-p2(3)*(Ti(j)-Ti(j-1))+p3(3)*(Ti(j+1)-Ti(j)))+Qreac*dt/rhoCp;
        elseif (inter==3)
            rhoCp=p3(1)*p3(2)+p4(1)*p4(2);
            dT=(dt/(1.5*dx*dx*rhoCp))*(-p3(3)*(Ti(j)-Ti(j-1))+p4(3)*(Ti(j+1)-Ti(j)))+Qreac*dt/rhoCp;
        elseif (inter==4)
            rhoCp=p4(1)*p4(2)+p5(1)*p5(2);
            dT=(dt/(1.5*dx*dx*rhoCp))*(-p4(3)*(Ti(j)-Ti(j-1))+p5(3)*(Ti(j+1)-Ti(j)))+Qreac*dt/rhoCp;
        end
        if (dT>dTmax)
            dTmax=dT;
        end
        Tf(j)=Ti(j)+dT;
    end
end
% if (dTmax>=dTlim)
%     if (dt>1E-7)
%         t=t-dt;
%         dt=0.1*dt;
%         Tf=Ti;
%         ynew=y;
%     end
% elseif (dTmax>0&&dTmax<0.1*dTlim)
%     dt=10*dt;
% end

    if (t>0.0001 && t<0.00011)
        plot(Tf)
        disp('ho gaya')
    end
end
function [rho, cp, k, reg, inter ]=prop(i, dx, L, p1, p2, p3, p4, p5)
%decides the location of the control volume aand spits out the relevant properties
%if the node is at interface the properties are set to zero and the
%interfdace type is decided
if (i<(L(1)/dx))
    rho=p1(1);
    cp=p1(2);
    k=p1(3);
    reg=1;
    inter=0;
elseif (i==(L(1)/dx))
    rho=0;
    cp=0;
    k=0;
    reg=0;
    inter=1;
elseif ((L(1)/dx)<i&&i<(((L(1)+L(2))/dx)-1))
    rho=p2(1);
    cp=p2(2);
    k=p2(3);
    reg=2;
    inter=0;
elseif (i==(((L(1)+L(2))/dx)-1))
    rho=0;
    cp=0;
    k=0;
    reg=0;
    inter=2;    
elseif ((((L(1)+L(2))/dx)-1)<i&&i<(((L(1)+L(2)+L(3))/dx)-2))
    rho=p3(1);
    cp=p3(2);
    k=p3(3);
    reg=3;
    inter=0;
elseif (i==(((L(1)+L(2)+L(3))/dx)-2))
    rho=0;
    cp=0;
    k=0;
    reg=0;
    inter=3;
elseif ((((L(1)+L(2)+L(3))/dx)-2)<i&&i<(((L(1)+L(2)+L(3)+L(4))/dx)-3))
    rho=p4(1);
    cp=p4(2);
    k=p4(3);
    reg=4;
    inter=0;
elseif (i==(((L(1)+L(2)+L(3)+L(4))/dx)-3))
    rho=0;
    cp=0;
    k=0;
    reg=0;
    inter=4;
elseif ((((L(1)+L(2)+L(3)+L(4))/dx)-3)<i)
    rho=p5(1);
    cp=p5(2);
    k=p5(3);
    reg=5;
    inter=0;
end
end

function [y_new,Qr]=reactions(y,T,dt, reg, inter, p2, p3, p4)
%decides the amount of heat released by the reactions 

%heat of different reactions
hcat1=77*1000;
hcat2=84*1000;
han1=1714*1000;
hSEI=257*1000;
hsep=-233.2*1000;
hele=800*1000;

%densities of materials
rho_cat=1244.74;
rho_an=1391.18;
rho_sep=295.2;
rho_ele_an=388.9;
rho_ele_cat=357.33;
rho_ele_sep=516;

%heat of reactions
ydot=rate(y,T); %ydot is generally negative
Qcat1=(-hcat1*ydot(1));
Qcat2=(-hcat2*ydot(2));
Qan1=(-han1*ydot(3));
if(T<533)
    SEI_dec=ydot(4)+5*ydot(3);
else
    SEI_dec=0;
end
QSEI=(-hSEI*SEI_dec);
Qsep=(-hsep*ydot(5));
Qele=(-hele*ydot(6));
if (T<1123)&&(T>533)
    Q_ISC=31720.7/0.72;
else
    Q_ISC=0;
end

% y_new = y*0.9;
% Qr = 1000;

if (reg==1)
    y_new=y;
    Qr=0;
    disp('i was here')
elseif ( reg==2 || ((reg==0)&&(inter==1)) )
    y_new=y+dt*ydot;
    Qr=rho_cat*Qcat1+rho_cat*Qcat2+rho_ele_cat*Qele+p2(1)*Q_ISC;
    disp('i was here 2')
elseif (((reg==0)&&(inter==2)))
    y_new=y+dt*ydot;
    Qr=rho_cat*Qcat1+rho_cat*Qcat2+rho_ele_cat*Qele+p2(1)*Q_ISC+rho_ele_sep*Qele+rho_sep*Qsep+p3(1)*Q_ISC;
    disp('i was here 3')
elseif (reg==3)
    y_new=y+dt*ydot;
    Qr=rho_ele_sep*Qele+rho_sep*Qsep+p3(1)*Q_ISC;
    disp('i was here4')
elseif (((reg==0)&&(inter==3)))
    y_new=y+dt*ydot;
    Qr=rho_ele_sep*Qele+rho_sep*Qsep+p3(1)*Q_ISC+rho_an*(Qan1+QSEI)+rho_ele_an*Qele+p4(1)*Q_ISC;
    disp('i was here5')
elseif (reg==4||((reg==0)&&(inter==4)))
    y_new=y+dt*ydot;
    Qr=rho_an*(Qan1+QSEI)+rho_ele_an*Qele+p4(1)*Q_ISC;
    disp('i was here6')
elseif (reg==5)
    y_new=y;
    Qr=0;
    disp('i was here7')
end 

end

function ydot=rate(y,T)
%takes in the concentration and temperature and spits out the rate at this
%concentration value for ALL the reactions. The computation of heat of
%reactions and the selection of those heats depending on the
%region/interface is done by another function called reactions

   c1=y(1); %cathode 1 reaction progress
   c2=y(2); %cathode 2 reaction progress
   a1=y(3); %anode decomposition progress
   aSEI=y(4); %SEI decomposition progress
   sep=y(5); %seperator meltdown progress
   el=y(6); %electrolyte decomposition progress
    
   R=8.314; %Gas COnstant
%    m_an=100.58; %anode Mass
%    m_cat=179.12;%Cathode mass
%    m_el=108;    %electrolyte mass
%    m_sep=17.6;  %separator mass
   
   ydot=zeros(6,1);
   
   %REACTION C1
   Ac1=1.75E12;
   Ec1=1.495E5;
   hc1=77;
   if (0<c1)&&(c1<1)&&(T>453)
       ydot(1)=-Ac1*(c1)*(1-c1)*exp(-Ec1/(R*T));
   else 
       ydot(1)=0;
   end
  % Qc1=-m_cat*(hc1*ydot(1));
   
   %REACTION C2
   Ac2=1.077E12;
   Ec2=1.5888E5;
   hc2=84;
   if (0<c2)&&(c2<1)&&(T>493)
       ydot(2)=-Ac2*c2*((1-c2))*exp(-Ec2/(R*T));
   else 
       ydot(2)=0;
   end
   %Qc2=-m_cat*(hc2*ydot(2));
   
   %REACTION Anode decomposition
   if (T<533)
       Aa1=0.035;
   else
       Aa1=5;
   end
   Ea1=3.3E4;
   ha1=1714;
   if (0<a1)&&(T>323)
       ydot(3)=-Aa1*a1*exp(-Ea1/(R*T))*exp(-aSEI);
   else 
       ydot(3)=0;
   end
   %Qa1=-m_an*(ha1*ydot(3));
  
   %REACTION SEI decompostion and regeneration
   ASEI=1.667E15;
   ESEI=1.3508E5;
   hSEI=257;
   if (0<aSEI)&&(T>323)&&(T<533)
       SEI_dec=-ASEI*aSEI*exp(-ESEI/(R*T));
       SEI_gen=5*ydot(3);
   else 
       SEI_dec=0;
       SEI_gen=0;
   end
   ydot(4)=SEI_dec-SEI_gen;
  % QSEI=-m_an*(-hSEI*SEI_dec);
   
   %REACTION seperator
   Asep=1.5E50;
   Esep=4.2E5;
   hsep=-233.2;
   if (0<sep)&&(T>393)
       ydot(5)=-Asep*sep*exp(-Esep/(R*T));
   else 
       ydot(5)=0;
   end
  % Qsep=-m_sep*(hsep*ydot(5));
   
   %REACTION electrolyte 
   Ael=3E15;
   Eel=1.7E5;
   hel=800;
   if (0<el)&&(T>433)
       ydot(6)=-Ael*el*exp(-Eel/(R*T));
   else 
       ydot(6)=0;
   end
  %Qel=-m_el*hel*ydot(6);
   
   if(T>533)&&(T<1123)
       Q_ISC=31720.7/0.72;  %heat released per unit mass of battery
   else
       Q_ISC=0;
   end

end




%finally beginning the code

mbat=360;   %has to be adjusted for different datasets
Cp=1.100;     %(J/g/K) Specific taken from Feng et al (2015) 

T_in=80+273;    %(K) ie 80 deg C
dt=1;           %time step (s) 
ti=dt;           %initial time
dTdt=10/60;     %(C/s) external temperature increase rate = 10C/min
tend=120000;  

T=T_in;     %the variable temperature
delT=0;     %temperature change at each step
tf=ti;      %the variable time

c1=0.999;    %conc for reaciton1 at cathode 
c2=0.999;    %conc for reaciton2 at cathode 
cSEI=0.15;  %conc for SEId reaction at anode
can=1;
csep=1;
celec=1;

Q1=0;
Q2=0;
QSEI=0;
QAn=0;
Qsep=0;
Qelec=0;
j=1;
while (tf<=tend)
    if (tf>ti) 
        T(j)=T(j-1)+delT;
    else 
        T(j)=T_in;
    end

    % ----- CATHODE Heat ------  
    [Q1(j),c1]=reaC1(T(j),dt,c1);
    [Q2(j),c2]=reaC2(T(j),dt,c2);
    
    %----- ANODE Heat --------
    [QSEI(j),QAn(j), cSEI, can]=reac_an(T(j),dt,cSEI,can);
    
    %----- Separator meltdown ----
    [Qsep(j),csep]=reacSep(T(j),dt,csep);
    
    %----- Electrolyte decomposition ------
    [Qelec(j),celec]=reacElec(T(j),dt,celec);
    
    delT=dt*(Q1(j)+Q2(j)+QSEI(j)+QAn(j)+Qsep(j)+Qelec(j))/(mbat*Cp);
%     if ((delT/dt)<dTdt) 
%         delT=dTdt*dt; 
%         tf=tf+dt;
%     else 
%         tf=tf+round((delT/dTdt),0);
%     end
    tf=tf+dt;
    j=j+1;
end
%**************************************************************************
%......................Plotting and Stuff..................................
%**************************************************************************
time=(dt/60):(dt/60):(tend/60);
figure 

plot(time,T)
figure
hold on
plot(T,Q1)
plot(T,Q2)
plot(T,QSEI)
plot(T,QAn)
plot(T,Qsep)
plot(T,Qelec)
legend('Q1','Q2','QSEI','QAn','Qsep','Qelec')

%##########################################################################
%.......................Cathode reaction Functions.........................
%##########################################################################

function [Q, conc] =reaC1(T,dt,a) 
%Takes the temerature and the time for which the reaction has been going on
%as the input 
%This is the "Cathode oxidation Step 1" acc to Bilyaz et al (2020) 
A=1.75E9;
m=1;
n=1;
E=1.1495E5;
R=8.314;
% t1=0;   %Start time
% t2=tf;  %end time
% ao=0.03;%initial concn 
mass_Cat=179.12; %(g) assuming battery to be 100 g 
h_reac=77; %(J/g)
%[time,conc]=ode45(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)), [t1 t2], ao);
if (T>453)
    if (a>1||a<0) 
        rate=0;
    else
        rate=A*(a^m)*((1-a)^n)*exp(-E/(R*T));
    end
    conc=a-dt*rate;
    if (conc<0)
        conc=0;
    end
    Q=mass_Cat*h_reac*rate;
else
    Q=0;
    conc=a;
end
if (imag(Q)~=0)
    disp('culprit is C1')
end 
end

function [Q, conc]=reaC2(T,dt,a)  
%Takes the temerature and the time for which the reaction has been going on
%as the input 
%This is the "Cathode oxidation Step 2" acc to Bilyaz et al (2020) 
A=1.077E12;
m=1;
n=1;
E=1.5888E5;
R=8.314;
mass_Cat=179.12; %(g) assuming battery to be 100 g 
h_reac=84; %(J/g)
%[time,conc]=ode45(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)), [t1 t2], ao);
if (T>493)
    if (a>1||a<0) 
        rate=0;
    else
        rate= A*(a^m)*((1-a)^n)*exp(-E/(R*T));
    end
    conc=a-dt*rate;
    if (conc<0)
        conc=0;
    end
    Q=mass_Cat*h_reac*rate;
else
    Q=0;
    conc=a;
end
end

function [Q_SEI,Q_an,conc_SEI,conc_an]=reac_an(T,dt,a,c)  
%Takes the temerature and the time for which the reaction has been going on
%as the input 
%This is the "SEI decomposition reaction" acc to Bilyaz et al (2020) and is
A1=1.667E15;
E1=1.3508E5;
R=8.314;
mass_an=100.58; %(g) assuming battery to be 100 g 
h_reac1=257; %(J/g)
if (T>323)
if (a>0) 
    rate_dec=A1*(a)*exp(-E1/(R*T));
else
    rate_dec=0;
    a=0;
end
if (T<533)
    A2=0.035;
else
    A2=5;
end
E2=3.3E4;
h_reac2=1714; %(J/g)

if (c>0) 
    rate_an=A2*(c)*exp(-E2/(R*T))*exp(-a/1);
else
    rate_an=0;
    c=0;
end

rate_gen=5*rate_an;
conc_SEI=a-dt*(rate_dec-rate_gen);
if (conc_SEI<0)
    conc_SEI=0;
    rate_dec=a/dt;
end

conc_an=c-dt*rate_an;
if (conc_an<0)
    conc_an=0;
    %rate_an=c/dt;
end
Q_SEI=mass_an*h_reac1*rate_dec;
Q_an=mass_an*h_reac2*rate_an;
else
    conc_SEI=a;
    conc_an=c;
    Q_SEI=0;
    Q_an=0;
end
end

function [Q, conc]=reacSep(T,dt,a)  
%Takes the temerature and the time for which the reaction has been going on
%as the input 
%This is the "Cathode oxidation Step 2" acc to Bilyaz et al (2020) 
A=1.5E50;
m=1;
n=0;
E=4.2E5;
R=8.314;
mass_sep=17.6; %(g) assuming battery to be 100 g 
h_reac=-233.2; %(J/g)
%[time,conc]=ode45(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)), [t1 t2], ao);
if (T>393)
if (a>1||a<0) 
    rate=0;
else
    rate= A*(a^m)*((1-a)^n)*exp(-E/(R*T));
end
conc=a-dt*rate;
if (conc<0)
    conc=0;
end
Q=mass_sep*h_reac*rate;
else
    Q=0;
    conc=a;
end
end

function [Q, conc]=reacElec(T,dt,a)  
%Takes the temerature and the time for which the reaction has been going on
%as the input 
%This is the "Cathode oxidation Step 2" acc to Bilyaz et al (2020) 
A=3E15;
m=1;
n=0;
E=1.7E5;
R=8.314;
mass_elec=108; %(g) assuming battery to be 100 g 
h_reac=140; %(J/g)
%[time,conc]=ode45(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)), [t1 t2], ao);
if (T>413)
if (a>1||a<0) 
    rate=0;
else
    rate= A*(a^m)*((1-a)^n)*exp(-E/(R*T));
end
conc=a-dt*rate;
if (conc<0)
    conc=0;
end
Q=mass_elec*h_reac*rate;
else
    Q=0;
    conc=a;
end
end
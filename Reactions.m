%finally beginning the code
% ----- ANODE REACITONS --------

%Go for it Yashraj 

% ----- CATHODE REACTIONS ------ 


%Cathode+electrolyte (here LCO) 

T=800; %Putting a fixed temperature that i know will get a good graph for reaction C1 and C2
        %C3 needs a higher temperature around 700-800
       
[t_C1, conc_C1, rate_C1] =reaC1(T);
[t_C2, conc_C2, rate_C2] =reaC2(T);
[t_C3, conc_C3, rate_C3] =reaC3(T);

figure 
hold on
legend
%plot(t_C1,rate_C1)
%plot(t_C2,rate_C2)
plot(t_C3,rate_C3)

function [time, conc, rate] =reaC1(T) 
%Didnt take the time of reaction and inital conc^n as an input variable as
%we will be probably fixing it for a particular reaction
A=1.27E10;
m=0;
n=1;
p=2/3;
E=1.138E5;
R=8.314;
t1=0;   %Start time
t2=10;  %end time
ao=0.03;%initial concn 
[time,conc]=ode45(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)), [t1 t2], ao);
rate=arrayfun(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)),time, conc);
end

function [time, conc, rate]=reaC2(T)  
A=2.33E12;
m=0;
n=1;
p=2/3;
E=1.407E5;
R=8.314;
t1=0;   %Start time
t2=10;  %end time
ao=0.03;%initial concn
[time,conc]=ode45(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)), [t1 t2], ao);
rate=arrayfun(@(t,a) A*(a^m)*((1-a)^n)*((-log(1-a))^p)*exp(-E/(R*T)),time, conc);
end

function [time, conc, rate]=reaC3(T)  
A=5.71E9;
E=1.593E5;
R=8.314;
t1=0;   %Start time
t2=10;  %end time
ao=0;%initial concn
[time,conc]=ode45(@(t,a) A*(1-a)*exp(-E/(R*T)), [t1 t2], ao);
rate=arrayfun(@(t,a) A*(1-a)*exp(-E/(R*T)),time, conc);
end
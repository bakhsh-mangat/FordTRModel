function [T,Y] = RK45(dydt, trange, yinit, tol)
% function solves ordinary differential equation just like ode45.
% dydt = system of first order differential equations
% trange = range of the independent variable
% yinit = initial condition
% tol = error absolute tolerance, default value = 1e-6;
% Example:
% [t,y] = RK45(@(t,y)[y(2); (1-y(1)^2)*y(2)-y(1)],[0 20],[2; 0]);
% subplot(1,2,2);
% plot(t,y,'-o')
% title('Solution of van der Pol Equation (\mu = 1) with RK45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2');
% 
% [t2,y2] = ode45(@(t,y)[y(2); (1-y(1)^2)*y(2)-y(1)],[0 20],[2; 0]);
% subplot(1,2,1);
% plot(t2,y2,'-o')
% title('Solution of van der Pol Equation (\mu = 1) with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2');


    global K1 K3 K4 K5
    if(nargin < 4) 
        tol = 1e-6;
    end
    if size(yinit,2) == 1
        yinit = yinit';
    end
    dtmax = diff(trange)/1e-2;
    dtmin = diff(trange)/1e20;
    dt = dtmax;
    Y = yinit; T = trange(1);Yn = yinit'; t = T(1);
    stop = 0; 
    while(~stop)
        residual = Res(dydt, dt,t,Yn);
        if(residual < tol)
            t = t + dt;

            Yn = Yn + ((25 * K1 / 216) + (1408 * K3 / 2565) + (2197 * K4 / 4104) - (K5 / 5));
            T = [T;t];
            Y = [Y;Yn'];
        end
        delta = min([4, max([0.1, 0.84 * (tol /residual)^0.25])]);
        dt = min([dtmax, delta * dt]);
        if (dt < dtmin)
            error('stiff equation')
        end
        if (t >= trange(2)) 
            stop = 1;
        elseif(t+dt>trange(2))
            dt = trange(2)-t;
        end
    end
end

function res = Res(dydx, dx, x, y)
     global K1 K3 K4 K5
     XX1 = x; YY1 = y; dxx = dx;
     K1 = dxx * dydx(XX1, YY1);
     XX2 = XX1 + dxx / 4;
     YY2 = YY1 + K1 / 4;
     K2 = dxx * dydx(XX2, YY2);
     XX3 = XX1 + 3 * dxx / 8;
     YY3 = YY1 + (3 * K1 / 32) + (9 * K2 / 32);
     K3 = dxx * dydx(XX3, YY3);
     XX4 = XX1 + 12 * dxx / 13;
     YY4 = YY1 + (1932 * K1 / 2197) - (7200 * K2 / 2197) + (7296 * K3 / 2197);
     K4 = dxx * dydx(XX4, YY4);
     XX5 = XX1 + dxx;
     YY5 = YY1 + (439 * K1 / 216) - 8 * K2 + (3680 * K3 / 513) - (845 * K4 / 4104);
     K5 = dxx * dydx(XX5, YY5);
     XX6 = XX1 + dxx / 2;
     YY6 = YY1 - (8 * K1 / 27) + 2 * K2 - (3544 * K3 / 2565) + (1859 * K4 / 4104) - (11 * K5 / 40);
     K6 = dxx * dydx(XX6, YY6);
     res = norm((K1 / 360) - (128 * K3 / 4275) - (2197 * K4 / 75240) + (K5 / 50) + (2 * K6 / 55));
end
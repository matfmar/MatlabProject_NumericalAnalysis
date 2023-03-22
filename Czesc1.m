function [result] = Czesc1(n) %n=1,2,3,4,5 dla kolejnych wymuszeń. n=5 dla e=sin(t)
    global licznikWykresow;
    licznikWykresow = 0;
    R1 = 0.1;
    R2 = 10;
    C = 0.5;
    L1 = 3;
    L2 = 5;
    M = 0.8;
    di1dt = @(i1,i2,uC,E) ((1/((L1/M)-(M/L2)))*(((-R1/M)*i1)+((R2/L2)*i2)-((1/M)*uC)+((1/M)*E)));
    di2dt = @(i1,i2,uC,E) ((1/((M/L1) - (L2/M)))*(((-R1/L1)*i1)+((R2/M)*i2)-((1/L1)*uC)+((1/L1)*E)));
    %di1dt = @(i1,i2,uC,E) ((-1/M) * (R2*i2 + L2 * di2dt(i1,i2,uC,E)));
    duCdt = @(i1) ((1/C)*i1);
    global h;
    h = 0.0001;
    global tmax;
    tmax = 30;
    t = 0:h:tmax;
    E = zeros(1, length(t));
    i1 = zeros(1, length(t));
    i2 = zeros(1, length(t));
    uC = zeros(1, length(t));
    i = 1;
    E(1) = fE(n, t(1)); %funkcja fE liczy E(t) dla roznych (n) podpunktow tego zadania
    while (t(i) < tmax) %implementacja pętli opisanej powyżej
        E(i+1) = fE(n, t(i+1));
        i1(i+1) = i1(i) + (h * di1dt(i1(i),i2(i),uC(i),E(i)));
        i2(i+1) = i2(i) + (h * di2dt(i1(i),i2(i),uC(i),E(i)));
        uC(i+1) = uC(i) + (h * duCdt(i1(i)));
        i = i + 1;
    end
    %odtąd rysowanie wykresów
    Rysuj(t, uC, E, i1, i2);
    LiczIRysujAnal(n);%liczy i rysuje weryfikację półanalityczną
    LiczIRysujMacierze(n);%liczy i rysuje weryfikację drugą metodą numeryczną
    result = 'wykonano';
end

function [E] = fE(n, t)
    global h;
    switch (n(1))
        case 1
            t2 = (1 / h) * t(1);
            T = (1 / h) * 3;
            div = rem(t2, T);
            if (div < (T/2))    %dobieram takie h żeby T/2 było całkowite
                E = 120;
            else
                E = 0;
            end
        case 2
            E = 240*sin(t(1));
        case 3
            E = 210*sin(2 * pi * 5 * t(1));
        case 4
            E = 120*sin(2 * pi * 50 * t(1));
        case 5
            E = sin(t(1));
        otherwise
            E = 1000;
    end
end

function [result] = LiczIRysujAnal(n)   %liczy i rysuje met.półanalityczną
    global h;
    global tmax;
    t = 0:h:tmax;
    i1Anal = zeros(1,length(t));
    i2Anal = zeros(1,length(t));
    uAnal = zeros(1,length(t));
    if (n ~= 5)
        return;
    end
    for i=1:length(t)
        uAnal(i) = -0.0029943*exp(-2.0779*t(i)) + 3.52233*exp(-0.023*t(i))*sin(0.8184*t(i)) - 2.87857*sin(t(i)) + 0.449056*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.446062*cos(t(i)) + 0;
        i1Anal(i) = 0.00311093*exp(-2.0779*t(i)) - 0.224261*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.223031*sin(t(i)) + 1.43617*exp(-0.023*t(i))*cos(0.8184*t(i)) - 1.43928*cos(t(i)) + 0;
        i2Anal(i) = 0.0132769*exp(-2.0779*t(i)) - 0.0870447*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.0992511*sin(t(i)) + 0.0185062*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.0317831*cos(t(i)) + 0;
    end
    RysujAnal(t, uAnal, i1Anal, i2Anal);
    result = 'obliczono analitycznie';
end

function [result] = Rysuj(t, uC, E, i1, i2) %rysuje obliczenia met. Eulera
    global licznikWykresow;
    licznikWykresow = licznikWykresow + 1;
    figure(licznikWykresow);
    hold on;
    plot(t, i1, 'DisplayName', 'i1(t)');
    plot(t, i2, 'DisplayName', 'i2(t)');
    %plot(t, uC, 'DisplayName', 'uC(t)');
    yline(0);
    xlabel('t');
    ylabel('i');
    %ylabel('i(t) lub u(t)');
    title('i(t) - Euler');
    %title('i(t) i u(t) - Euler');
    legend;
    hold off;
    licznikWykresow = licznikWykresow + 1;
    figure(licznikWykresow);
    hold on;
    plot(t, E, 'DisplayName', 'E(t)');
    plot(t, uC, 'DisplayName', 'uC(t)');
    yline(0);
    xlabel('t');
    ylabel('u');
    title('u(t) - Euler');
    legend;
    hold off;
    result = 'narysowano';
end

function [result] = RysujAnal(t, uAnal, i1Anal, i2Anal)
    global licznikWykresow;
    licznikWykresow = licznikWykresow + 1;
    figure(licznikWykresow);
    hold on;
    plot(t, i1Anal, 'DisplayName', 'i1(t)');
    plot(t, i2Anal, 'DisplayName', 'i2(t)');
    plot(t, uAnal, 'DisplayName', 'E(t)');
    yline(0);
    xlabel('t');
    %ylabel('i');
    ylabel('i(t) lub u(t)');
    %title('i(t) - półanalitycznie');
    title('i(t) lub u(t) - półanalitycznie');
    legend;
    hold off;
    licznikWykresow = licznikWykresow + 1;
    figure(licznikWykresow);
    hold on;
    plot(t, uAnal, 'DisplayName', 'E(t)');
    yline(0);
    xlabel('t');
    ylabel('u');
    title('u(t) - półanalitycznie');
    legend;
    hold off;
    result = 'narysowano';
end

function [result] = LiczIRysujMacierze(n)   %liczy i rysuje met. Rungego-Kutty
    global licznikWykresow;
    global tmax;
    global h;
    if (n(1) == 5)
        [t,y] = ode45(@odefun, [0 tmax], [0 0 0]);
        licznikWykresow = licznikWykresow + 1;
        figure(licznikWykresow);
        plot(t, y);
        yline(0);
        xlabel('t');
        ylabel('i(t) lub u(t)');
        title('i(t) i u(t) - Runge-Kutta');
        legend('i1(t)','i2(t)','u(t)');
    end
    result = 'wykonano obliczenia macierzowe';
end

function dydt = odefun(t, y)
    R1 = 0.1;
    R2 = 10;
    C = 0.5;
    L1 = 3;
    L2 = 5;
    M = 0.8;
    dydt = [(1/(L1/M - M/L2))*((-R1/M)*y(1) + (R2/L2)*y(2) - (1/M)*y(3) + (1/M)*sin(t));
        (1/(M/L1 - L2/M))*((-R1/L1)*y(1) + (R2/M)*y(2) - (1/L1)*y(3) + (1/L1)*sin(t));
        (1/C)*y(1)];
end

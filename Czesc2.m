function [result] = Czesc2(n) %n=3,4,(5)
    wezly = [20,50,100,150,200,250,280,300];
    wartosci = [0.46,0.64,0.78,0.68,0.44,0.23,0.18,0.18];
    global h;
    h = 0.0001;
    global tmax;
    tmax = 30;
    t = 0:h:tmax;
    global sp;
    global figCounter;
    figCounter = 1;

    %Trial(wezly, wartosci); %liczy M(uL) dla roznych metod
    
    %wielomianowa interp (lagrange)
    sp = 1;
    [uC_lagr, i1_lagr, i2_lagr, E, R1, R2, uL1_lagr, uR2_lagr] = ObliczNumerycznie(n, wezly, wartosci, 0);
    if (n == 5)
        [tM_lagr, yM_lagr] = LiczMacierze(wezly, wartosci, 0);%metoda weryfikacyjna dla e=sin(t)
        i1M_lagr = yM_lagr(:,1);
        i2M_lagr = yM_lagr(:,2);
        uM_lagr = yM_lagr(:,3);
        figure(figCounter);
        figCounter = figCounter + 1;
        hold on;
        plot(tM_lagr, yM_lagr);
        title('mac,lagrange');
        legend('i1','i2','uC');
        hold off;
    end
%     figure(figCounter);
%     figCounter = figCounter+1;
%     hold on;
%     plot(t, i1_lagr);
%     plot(t, i2_lagr);
%     plot(t, uC_lagr);
%     title('interpolacja wielomianowa');
%     hold off;
    
    %spline interp
    sp = 2;
    [uC_spl, i1_spl, i2_spl, E, R1, R2, uL1_spl, uR2_spl] = ObliczNumerycznie(n, wezly, wartosci, 0);
    if (n == 5)
        [tM_spl, yM_spl] = LiczMacierze(wezly, wartosci, 0);
        i1M_spl = yM_spl(:,1);
        i2M_spl = yM_spl(:,2);
        uM_spl = yM_spl(:,3);
        figure(figCounter);
        figCounter = figCounter + 1;
        hold on;
        plot(tM_spl, yM_spl);
        title('mac,spline');
        legend('i1','i2','uC');
        hold off;
    end
%     figure(figCounter);
%     figCounter = figCounter+1;
%     hold on;
%     plot(t, i1_spl);
%     plot(t, i2_spl);
%     plot(t, uC_spl);
%     title('interpolacja sklejana');
%     hold off;

    %aprox 3
    sp = 3;
    A3 = AproksymacjaStopien(wezly, wartosci, 3);%najpierw współczynniki
    [uC_apr3, i1_apr3, i2_apr3, E, R1, R2, uL1_apr3, uR2_apr3] = ObliczNumerycznie(n, wezly, wartosci, A3);
    if (n == 5)
        [tM_apr3, yM_apr3] = LiczMacierze(wezly, wartosci, A3);
        i1M_apr3 = yM_apr3(:,1);
        i2M_apr3 = yM_apr3(:,2);
        uM_apr3 = yM_apr3(:,3);
        figure(figCounter);
        figCounter = figCounter + 1;
        hold on;
        plot(tM_apr3, yM_apr3);
        title('mac,apr3');
        legend('i1','i2','uC');
        hold off;
    end
%     figure(figCounter);
%     figCounter = figCounter+1;
%     hold on;
%     plot(t, i1_apr3);
%     plot(t, i2_apr3);
%     plot(t, uC_apr3);
%     title('interpolacja aproks. 3 st.');
%     hold off;

    %aprox 5
    sp = 4;
    A5 = AproksymacjaStopien(wezly, wartosci, 5);%najpierw współczynniki
    [uC_apr5, i1_apr5, i2_apr5, E, R1, R2, uL1_apr5, uR2_apr5] = ObliczNumerycznie(n, wezly, wartosci, A5);
    if (n == 5)
        [tM_apr5, yM_apr5] = LiczMacierze(wezly, wartosci, A5);
        i1M_apr5 = yM_apr5(:,1);
        i2M_apr5 = yM_apr5(:,2);
        uM_apr5 = yM_apr5(:,3);
        figure(figCounter);
        figCounter = figCounter + 1;
        hold on;
        plot(tM_apr5, yM_apr5);
        title('mac,apr5');
        legend('i1','i2','uC');
        hold off;
    end
%     figure(figCounter);
%     figCounter = figCounter+1;
%     hold on;
%     plot(t, i1_apr5);
%     plot(t, i2_apr5);
%     plot(t, uC_apr5);
%     title('interpolacja aproks. 5 st.');
%     hold off;

    % poniższe funkcje rysują wykresy
    RysowanieRoznic(t, i1_lagr, i1_spl, i1_apr3, i1_apr5, 'i1');
    RysowanieRoznic(t, i2_lagr, i2_spl, i2_apr3, i2_apr5, 'i2');
    RysowanieRoznic(t, uC_lagr, uC_spl, uC_apr3, uC_apr5, 'uC');
    RysowanieRoznic(t, uR2_lagr, uR2_spl, uR2_apr3, uR2_apr5, 'uR2');
    RysowanieRoznic(t, uL1_lagr, uL1_spl, uL1_apr3, uL1_apr5, 'uL1');

    result = 'wykonano';
end

function [uC_2, i1_2, i2_2, E_2, r1, r2, uL1_2, uR2_2] = ObliczNumerycznie(n, wezly, wartosci, A)
    global h;
    global tmax;
    R1 = 0.1;
    R2 = 10;

    C = 0.5;
    L1 = 3;
    L2 = 5;
    M = 0.8;

    di1dt = @(i1,i2,uC,E,mm) ((1/((L1/mm)-(mm/L2)))*(((-R1/mm)*i1)+((R2/L2)*i2)-((1/mm)*uC)+((1/mm)*E)));
    di2dt = @(i1,i2,uC,E,mm) ((1/((mm/L1) - (L2/mm)))*(((-R1/L1)*i1)+((R2/mm)*i2)-((1/L1)*uC)+((1/L1)*E)));
    %di1dt = @(i1,i2,uC,E) ((-1/M) * (R2*i2 + L2 * di2dt(i1,i2,uC,E)));
    duCdt = @(i1) ((1/C)*i1);

    t = 0:h:tmax;
    E = zeros(1, length(t));
    i1 = zeros(1, length(t));
    i2 = zeros(1, length(t));
    uC = zeros(1, length(t));
    uL1 = zeros(1, length(t));
    uR2 = zeros(1, length(t));
    i = 1;
    E(1) = fE(n, 0); %funkcja fE liczy E(t) dla roznych (n) podpunktow tego zadania (na dole)
    while (t(i) < tmax)
        mm = fM(wezly, wartosci, uL1(i), A);% to daje M
        E(i+1) = fE(n, t(i+1));
        poch = di1dt(i1(i),i2(i),uC(i),E(i),mm);
        i1(i+1) = i1(i) + (h * poch);
        i2(i+1) = i2(i) + (h * di2dt(i1(i),i2(i),uC(i),E(i),mm));
        uC(i+1) = uC(i) + (h * duCdt(i1(i)));
        uR2(i+1) = i2(i+1)*R2;
        uL1(i+1) = abs(L1 * poch);
        i = i + 1;
    end
    uC_2 = uC;
    i1_2 = i1;
    i2_2 = i2;
    E_2 = E;
    r1 = R1;
    r2 = R2;
    uL1_2 = uL1;
    uR2_2 = uR2;
end

function [E] = fE(n, t)%funkcja wyliczająca wymuszenie
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
            E = 240*sin(2 * t(1));
        case 4
            E = 120*sin(2 * t(1));
        case 5
            E = sin(t(1));
        otherwise
            E = 1000;
    end
end

function [t, y] = LiczMacierze(wezly, wartosci, A)
    global tmax;
    global h;
    [t,y] = ode45(@(t,y) odefun(t,y,wezly,wartosci, A), [0 tmax], [0 0 0]);
end

function dydt = odefun(t, y, wezly, wartosci, A)
    R1 = 0.1;
    R2 = 10;
    C = 0.5;
    L1 = 3;
    L2 = 5;
    M = 0.8;

    dydt = zeros(3,1);
    dydt(1) = (1/((L1/fM(wezly, wartosci, abs(L1*dydt(1)), A)) - (fM(wezly, wartosci, abs(L1*dydt(1)), A)/L2)))*((-R1/fM(wezly, wartosci, abs(L1*dydt(1)), A))*y(1) + (R2/L2)*y(2) - (1/fM(wezly, wartosci, abs(L1*dydt(1)),A))*y(3) + (1/fM(wezly, wartosci, abs(L1*dydt(1)),A))*sin(t));
    dydt(2) =  (1/(fM(wezly, wartosci, abs(L1*dydt(1)),A)/L1 - L2/fM(wezly, wartosci, abs(L1*dydt(1)),A)))*((-R1/L1)*y(1) + (R2/fM(wezly, wartosci, abs(L1*dydt(1)),A))*y(2) - (1/L1)*y(3) + (1/L1)*sin(t));
    dydt(3) = (1/C)*y(1);

end

function [M] = fM(wezly, wartosci, uL, A)%w zależności od sposobu inter/apr odpowiednia funkcja jest wywoływana
    global sp;
    if (sp == 1)
        M = Lagrange(wezly, wartosci, uL);
    elseif (sp == 2)
        M = Spline(wezly, wartosci, uL);
    elseif (sp == 3 || sp == 4)
        M = LiczAproksymacje(uL, A);
    else
        M = 0.8;
    end
end

function [y] = Lagrange(wezly, wartosci, x)%klasyczna funkcja licząca współczynniki met.Lagrange'a
    wynik = 0;
    for i=1:length(wezly)
        iloczyn = 1;
        for j=1:length(wezly)
            if (j == i)
                continue;
            end
            iloczyn = iloczyn * ((x-wezly(j)) / (wezly(i)-wezly(j)));
        end
        wynik = wynik + (iloczyn * wartosci(i));
    end
    y = wynik;
end

function [a,b,skr] = ZnajdzPrzedzialSpline(wezly, x)%kwalifikuje uL do przedziału
    a=1;
    b=1;
    skr=0;
    if (x(1) < wezly(1))
        a = 1;
        b = 2;
        skr = 0;
        return;
    end
    if (x(1) >= wezly(length(wezly)))
        a = length(wezly);
        b = a;
        skr = 1;
        return;
    end
    for i=1:length(wezly)
        if (x(1) == wezly(i))
            a = i;
            b = i;
            skr = 1;
            return;
        end
        if (x(1) > wezly(i) && x(1) < wezly(i+1))
            a = i;
            b = i + 1;
            skr = 0;
            return;
        end
    end
end

function [y] = Spline(wezly, wartosci, x)
    [a, b, skr] = ZnajdzPrzedzialSpline(wezly, x);%aby wyliczyć M trzeba 'zakwalifikować' uL do odpowiedzniego przedziału
    if (skr == 1)
        y = wartosci(a);
        return;
    end
    rozw = zeros(2, 1);
    left = [wezly(a), 1; wezly(b), 1];
    right = [wartosci(a); wartosci(b)];
    rozw = linsolve(left, right);%rozwiązuje prosty układ 2x2
    A = rozw(1);
    B = rozw(2);
    y = A * x(1) + B;
end

function [A] = AproksymacjaStopien(wezly, wartosci, stopien)
    st = stopien(1);
    M = zeros(length(wezly), st+1);%macierz funkcji psi(wezel)
    for i=1:length(wezly)
        for j=1:(st+1)
            M(i, j) = wezly(i)^(j-1);
        end
    end
    Y = transpose(wartosci);
    A = zeros(st+1, 1);
    left = transpose(M) * M;
    right = transpose(M) * Y;
    A = linsolve(left, right);
end

function [y] = LiczAproksymacje(x, A)%mnozy wspólczynniki przez x i sumuje
    y = 0;
    for i=1:length(A)
        y = y + A(i)*x(1)^(i-1);
    end
end

function [result] = Trial(wezly, wartosci)
    global figCounter;
    x = 0:1:400;
    yLagrange = zeros(1, length(x));
    ySpline = zeros(1, length(x));
    yApr3 = zeros(1, length(x));
    yApr5 = zeros(1, length(x));
    A3 = AproksymacjaStopien(wezly, wartosci, 3);
    A5 = AproksymacjaStopien(wezly, wartosci, 5);
    for i=1:length(x)
        yLagrange(i) = Lagrange(wezly, wartosci, x(i));
        ySpline(i) = Spline(wezly, wartosci, x(i));
        yApr3(i) = LiczAproksymacje(x(i), A3);
        yApr5(i) = LiczAproksymacje(x(i), A5);
    end
    figure(figCounter);
    figCounter = figCounter + 1;
    hold on;
    plot(x, yLagrange);
    plot(x, ySpline);
    plot(x, yApr3);
    plot(x, yApr5);
    title('M(uL)');
    legend('lagr','spl','apr3','apr5');
    hold off;  
    result = 'wykonano';
end

function [result] = RysowanieRoznic(t, a1, a2, a3, a4, nazwa)
    global figCounter;
    figure(figCounter);
    figCounter = figCounter + 1;
    hold on;
    plot(t, a1);
    plot(t, a2);
    plot(t, a3);
    plot(t, a4);
    title(nazwa);
    xlabel('t');
    ylabel(nazwa);
    legend('inter. wielom.','inter. sklej.', 'aproks. 3st.', 'aproks. 5st.');
    result = 'wykonano';
end
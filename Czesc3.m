function [result] = Czesc3(n)
    global sp;%dla zmiennego sprzeżenia.
    sp = 1;
    global h;
    h = 0.00001;
    global Ctmax;
    tmax = 30;
    t = 0:h:tmax;
    wezly = [20,50,100,150,200,250,280,300];%dla zmiennego sprzężenia
    wartosci = [0.46,0.64,0.78,0.68,0.44,0.23,0.18,0.18];
    %start liczenia pradow
    [uC, i1, i2, E, R1, R2] = ObliczNumerycznie(n);
    %[uC, i1, i2, E, R1, R2] = ObliczNumerycznie(n, wezly, wartosci, 0);
    [uCAnal, i1Anal, i2Anal, energia] = ObliczAnal(n);
    [tM, yM] = LiczMacierze(n);
    i1M = yM(:,1);
    i2M = yM(:,2);
    %koniec liczenia pradow
    
    %całka - prostokąt - półanalitycznie i numerycznie
    calkaProstokat1 = 0;
    calkaProstokat2 = 0;
    calkaProstokat = 0;
    calkaParabola1 = 0;
    calkaParabola2 = 0;
    calkaParabola = 0;
    calkaProstokatAnal1 = 0;
    calkaProstokatAnal2 = 0;
    calkaProstokatAnal = 0;
    for i=1:length(t)
        calkaProstokat1 = calkaProstokat1 + h*(i1(i)^2);
        calkaProstokat2 = calkaProstokat2 + h*(i2(i)^2);
        calkaProstokatAnal1 = calkaProstokatAnal1 + h*(i1Anal(i)^2);
        calkaProstokatAnal2 = calkaProstokatAnal2 + h*(i2Anal(i)^2);
    end
    calkaProstokat1 = calkaProstokat1 * R1;
    calkaProstokat2 = calkaProstokat2 * R2;
    calkaProstokat = calkaProstokat1 + calkaProstokat2;
    calkaProstokatAnal1 = calkaProstokatAnal1 * R1;
    calkaProstokatAnal2 = calkaProstokatAnal2 * R2;
    calkaProstokatAnal = calkaProstokatAnal1 + calkaProstokatAnal2;
    
    %całka - parabola
    i = 1;
    while (i < length(t) - 1)
        calkaParabola1 = calkaParabola1 + (h/3)*((i1(i)^2) + 4*(i1(i+1)^2) + (i1(i+2)^2));
        calkaParabola2 = calkaParabola2 + (h/3)*((i2(i)^2) + 4*(i2(i+1)^2) + (i2(i+2)^2));
        i = i + 2;
    end
    calkaParabola1 = calkaParabola1 * R1;
    calkaParabola2 = calkaParabola2 * R2;
    calkaParabola = calkaParabola1 + calkaParabola2;

    %całka - metodą weryfikacyjną
    calkaProstMac1 = 0;
    calkaProstMac2 = 0;
    calkaProstMac = 0;
    for i=1:(length(tM)-1)
        calkaProstMac1 = calkaProstMac1 + ((tM(i+1)-tM(i))*i1M(i)^2);
        calkaProstMac2 = calkaProstMac2 + ((tM(i+1)-tM(i))*i2M(i)^2);
    end
    calkaProstMac1 = calkaProstMac1 * R1;
    calkaProstMac2 = calkaProstMac2 * R2;
    calkaProstMac = calkaProstMac1 + calkaProstMac2;

    mocProstokat = calkaProstokat / tmax;
    mocParabola = calkaParabola / tmax;
    mocProstokatAnal = calkaProstokatAnal / tmax;
    mocIntegralAnal = energia / tmax;
    mocMacierze = calkaProstMac / tmax;

    result = sprintf('Energia[J]:\nprostokaty: %.3f\nparabole: %.3f\nnumerycznie: %.3f\nintegral(): %.3f\nmacierze: %.3f\nMoc[W]:\nprostokaty: %.3f\nparabole: %.3f\nnumerycznie: %.3f\nintegral(): %.3f\nmacierze: %.3f\n',calkaProstokat, calkaParabola, calkaProstokatAnal, energia, calkaProstMac, mocProstokat, mocParabola, mocProstokatAnal, mocIntegralAnal, mocMacierze);
end

% function [uC_2, i1_2, i2_2, E_2, r1, r2] = ObliczNumerycznie(n, wezly,
% wartosci, A)%funkcja dla zmiennego sprzężenia
%     global h;
%     global tmax;
%     R1 = 0.1;
%     R2 = 10;
% 
%     C = 0.5;
%     L1 = 3;
%     L2 = 5;
%     M = 0.8;
% 
%     di1dt = @(i1,i2,uC,E,mm) ((1/((L1/mm)-(mm/L2)))*(((-R1/mm)*i1)+((R2/L2)*i2)-((1/mm)*uC)+((1/mm)*E)));
%     di2dt = @(i1,i2,uC,E,mm) ((1/((mm/L1) - (L2/mm)))*(((-R1/L1)*i1)+((R2/mm)*i2)-((1/L1)*uC)+((1/L1)*E)));
%     %di1dt = @(i1,i2,uC,E) ((-1/M) * (R2*i2 + L2 * di2dt(i1,i2,uC,E)));
%     duCdt = @(i1) ((1/C)*i1);
% 
%     t = 0:h:tmax;
%     E = zeros(1, length(t));
%     i1 = zeros(1, length(t));
%     i2 = zeros(1, length(t));
%     uC = zeros(1, length(t));
%     uL1 = zeros(1, length(t));
%     uR2 = zeros(1, length(t));
%     i = 1;
%     E(1) = fE(n, 0); %funkcja fE liczy E(t) dla roznych (n) podpunktow tego zadania (na dole)
%     while (t(i) < tmax)
%         mm = fM(wezly, wartosci, uL1(i), A);% to daje M
%         E(i+1) = fE(n, t(i+1));
%         poch = di1dt(i1(i),i2(i),uC(i),E(i),mm);
%         i1(i+1) = i1(i) + (h * poch);
%         i2(i+1) = i2(i) + (h * di2dt(i1(i),i2(i),uC(i),E(i),mm));
%         uC(i+1) = uC(i) + (h * duCdt(i1(i)));
%         uR2(i+1) = i2(i+1)*R2;
%         uL1(i+1) = abs(L1 * poch);
%         i = i + 1;
%     end
%     uC_2 = uC;
%     i1_2 = i1;
%     i2_2 = i2;
%     E_2 = E;
%     r1 = R1;
%     r2 = R2;
% end
% 
% function [M] = fM(wezly, wartosci, uL, A)%w zmienne M
%     global sp;
%     if (sp == 1)
%         M = Lagrange(wezly, wartosci, uL);
%     else
%         M = 0.8;
%     end
% end
% 
% function [y] = Lagrange(wezly, wartosci, x)%klasyczna funkcja licząca współczynniki met.Lagrange'a
%     wynik = 0;
%     for i=1:length(wezly)
%         iloczyn = 1;
%         for j=1:length(wezly)
%             if (j == i)
%                 continue;
%             end
%             iloczyn = iloczyn * ((x-wezly(j)) / (wezly(i)-wezly(j)));
%         end
%         wynik = wynik + (iloczyn * wartosci(i));
%     end
%     y = wynik;
% end

function [uC_2, i1_2, i2_2, E_2, r1, r2] = ObliczNumerycznie(n)
    global h;
    global tmax;
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

    t = 0:h:tmax;
    E = zeros(1, length(t));
    i1 = zeros(1, length(t));
    i2 = zeros(1, length(t));
    uC = zeros(1, length(t));
    i = 1;
    E(1) = fE(n, 0); %funkcja fE liczy E(t) dla roznych (n) podpunktow tego zadania (na dole)
    while (t(i) < tmax)
        t(i+1) = t(i) + h;
        E(i+1) = fE(n, t(i+1));%E(t) rysuje sie dobrze
        i1(i+1) = i1(i) + (h * di1dt(i1(i),i2(i),uC(i),E(i)));
        i2(i+1) = i2(i) + (h * di2dt(i1(i),i2(i),uC(i),E(i)));
        uC(i+1) = uC(i) + (h * duCdt(i1(i)));
        i = i + 1;
    end
    uC_2 = uC;
    i1_2 = i1;
    i2_2 = i2;
    E_2 = E;
    r1 = R1;
    r2 = R2;
end

function [t, y] = LiczMacierze(n)
    global tmax;
    global h;
    if (n(1) == 5)
        [t,y] = ode45(@odefun, [0 tmax], [0 0 0]);
    else
        t = [1 1];
        y = [1 1];
    end
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

function [u, i1, i2, energia] = ObliczAnal(n)
    global h;
    global tmax;
    t = 0:h:tmax;
    energia = 0;
    R1 = 0.1;
    R2 = 10;
    i1Anal = zeros(1,length(t));
    i2Anal = zeros(1,length(t));
    uAnal = zeros(1,length(t));
    for i=1:length(t)
        if (n == 4)
        %dla e=120sin(2pi50t)
        uAnal(i) = -0.00605859*exp(-2.0779*t(i)) + 0.449743*exp(-0.023*t(i))*sin(0.8184*t(i)) + i*(5.02926*10^-22*sin(314.759*t(i)) - 3.63769*10^-19 * exp(-2.0779*t(i))) - 0.00120892*sin(314.759*t(i)) + 0.00605907*exp(-0.023*t(i)) * cos(0.8184*t(i)) - 4.75862*10^-7 * cos(314.759*t(i));
        i1Anal(i) = 0.00645155*exp(-2.0799*t(i)) - 0.00770643*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.0000760998*sin(314.759*t(i)) + 0.183808*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.19026*cos(314.759*t(i)) + 0;
        i2Anal(i) = 0.0268642*exp(-2.0779*t(i)) - 0.0106901*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.000205402*sin(314.759*t(i)) + 0.00357606*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.0304402*cos(314.759*t(i)) + 0;
        elseif (n == 3)        
        %dla e=120sin(2pi5t)
        uAnal(i) = -0.105971*exp(-2.0779*t(i)) + 7.90587*exp(-0.023*t(i))*sin(0.8184*t(i)) - 0.212882*sin(31.416*t(i)) + 0.106809*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.000837577*cos(31.416*t(i)) + i*(1.46367*10^-18 * sin(31.416*t(i)) - 2.77556*10^-17 * cos(31.416*t(i)));
        i1Anal(i) = 0.110099*exp(-2.0779*t(i)) - 0.134624*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.0131567*sin(31.416*t(i)) + 3.23385*exp(-0.023*t(i))*cos(0.8184*t(i)) - 3.34395*cos(31.416*t(i)) + 0;
        i2Anal(i) = -0.469884*exp(-2.0779*t(i)) + 0.18792*exp(-0.023*t(i))*sin(0.8184*t(i)) - 0.0360202*sin(31.416*t(i)) - 0.0628553*exp(-0.023*t(i))*cos(0.8184*t(i)) + 0.532739*cos(31.416*t(i)) + i*(-5.55112*10^-17 * sin(31.416*t(i)) - 1.04083*10^-17 * cos(31.416*t(i)));
        elseif (n == 5)
        %dla e=sin(t)
        uAnal(i) = -0.0029943*exp(-2.0779*t(i)) + 3.52233*exp(-0.023*t(i))*sin(0.8184*t(i)) - 2.87857*sin(t(i)) + 0.449056*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.446062*cos(t(i)) + 0;
        i1Anal(i) = 0.00311093*exp(-2.0779*t(i)) - 0.224261*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.223031*sin(t(i)) + 1.43617*exp(-0.023*t(i))*cos(0.8184*t(i)) - 1.43928*cos(t(i)) + 0;
        i2Anal(i) = 0.0132769*exp(-2.0779*t(i)) - 0.0870447*exp(-0.023*t(i))*sin(0.8184*t(i)) + 0.0992511*sin(t(i)) + 0.0185062*exp(-0.023*t(i))*cos(0.8184*t(i)) - 0.0317831*cos(t(i)) + 0;
        i1fun = @(x) (0.00311093.*exp(-2.0779.*x) - 0.224261.*exp(-0.023.*x).*sin(0.8184.*x) + 0.223031.*sin(x) + 1.43617.*exp(-0.023.*x).*cos(0.8184.*x) - 1.43928.*cos(x) + 0).^2;
        i2fun = @(x) (0.0132769.*exp(-2.0779.*x) - 0.0870447.*exp(-0.023.*x).*sin(0.8184.*x) + 0.0992511.*sin(x) + 0.0185062.*exp(-0.023.*x).*cos(0.8184.*x) - 0.0317831.*cos(x) + 0).^2;
        else
        uAnal(i) = 1000;
        i1Anal(i) = 1000;
        i2Anal(i) = 1000;
        end
    end

    u = uAnal;
    i1 = i1Anal;
    i2 = i2Anal;
    if (n == 5)
        energia = R1*integral(i1fun,0,tmax) + R2*integral(i2fun,0,tmax);
    end
end
%funkcja do liczenia pradow
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
        case 6
            E = 1;
        otherwise
            E = 1000;
    end
end
%funkcja do liczenia pradow

function [result] = Czesc1_2(n)
    R1 = 0.1;
    R2 = 10;
    C = 0.5;
    L1 = 3;
    L2 = 5;
    M = 0.8;
    global h;
    h = 0.0001;
    global tmax;
    tmax = 30;
    [i1, i2, uC, E, t] = MidPointEuler(n);

    figure(1);
    hold on;
    plot(t, E, 'DisplayName', 'E(t)');
    plot(t, uC, 'DisplayName', 'uC(t)');
    yline(0);
    xlabel('t');
    ylabel('u');
    title('u(t) - wariant metody Eulera');
    legend;
    figure(2);
    hold on;
    plot(t, i1, 'DisplayName', 'i1(t)');
    plot(t, i2, 'DisplayName', 'i2(t)');
    xlabel('t');
    ylabel('i');
    title('i(t) - wariant metody Eulera');
    legend;
    hold off;
    result = 'wykonano';
end

function [i1_2, i2_2, uC_2, E_2, t_2] = MidPointEuler(n)
    R1 = 0.1;
    R2 = 10;
    C = 0.5;
    L1 = 3;
    L2 = 5;
    M = 0.8;    
    global h;
    global tmax;
    t = 0:h:tmax;
    di1dt = @(i1,i2,uC,E) ((1/((L1/M)-(M/L2)))*(((-R1/M)*i1)+((R2/L2)*i2)-((1/M)*uC)+((1/M)*E)));
    di2dt = @(i1,i2,uC,E) ((1/((M/L1) - (L2/M)))*(((-R1/L1)*i1)+((R2/M)*i2)-((1/L1)*uC)+((1/L1)*E)));
    duCdt = @(i1) ((1/C)*i1);
    E = zeros(1,length(t));
    i1 = zeros(1,length(t));
    i2 = zeros(1, length(t));
    uC = zeros(1, length(t));
    i = 1;
    E(1) = fE(n, t(1));
    while (t(i) < 30)
        E(i+1) = fE(n, t(i+1));
        i1c = di1dt(i1(i),i2(i),uC(i),E(i));
        i1(i+1) = i1(i) + (h * di1dt(i1(i)+(h/2)*i1c,i2(i)+(h/2)*i1c, uC(i)+(h/2)*i1c, E(i)+(h/2)*i1c));
        i2c = di2dt(i1(i),i2(i),uC(i),E(i));
        i2(i+1) = i2(i) + (h * di2dt(i1(i)+(h/2)*i2c,i2(i)+(h/2)*i2c,uC(i)+(h/2)*i2c,E(i)+(h/2)*i2c));
        uC2 = duCdt(i1(i));
        uC(i+1) = uC(i) + (h * duCdt(i1(i)+(h/2)*uC2));
        i = i + 1;
    end
    i1_2 = i1;
    i2_2 = i2;
    uC_2 = uC;
    E_2 = E;
    t_2 = t
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

function sth=Q2()

%{
    Rx is the Runge function
    function spline3 is the Spline interpolation function in second-order boundary condition
    function lpline is the Lagrange polynomial interpolation
    function pline is the Polynomial interpolation
    function npline4 is the Newton polynomial interpolation
%}


    format long;

    n = 20;

    Rx = @(x) 1./(1+25*x.^2);
    x = linspace(-1, 1, n);
    y = Rx(x);
   
    x_val = -1.0:0.001:1.0;
    S_val = spline3(x, y, 0, 0, x_val);
    L_val = lpline(x, y, x_val);
    P_val = pline(x, y, x_val);

    fg = figure(1);

    plot(x_val, P_val, 'c', x_val, S_val, 'r', ...
        x_val, L_val, 'b--', ...
        x, y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    % axis([0, 1, 0, 1]);
    grid on;
    set(fg, 'Position', [0 0 1000 600]);
    title(['n = ' int2str(n)]);
    % xlabel('X');ylabel('Y');

    
    
    function S_val = spline3(x, y, f0, fn, x_val)
        format long;
        
        k = length(x);

        h = zeros(1, k-1);
        % tmp: primary difference quotient
        tmp = zeros(1, k-1);
        for i=1:(k-1)
            h(i) = (x(i+1)-x(i));
            tmp(i) = (y(i+1)-y(i))/h(i);
        end

        D = zeros(k, 1);
        for i=2:(k-1)
            D(i) = 6*(tmp(i)-tmp(i-1))/(h(i-1)+h(i));
        end
        D(1) = f0;
        D(k) = fn;

        miu = zeros(1, k-1);
        lmd = zeros(1, k-1);
        for i=1:(k-2)
            miu(i) = h(i)/(h(i)+h(i+1));
            lmd(i+1) = 1-miu(i);
        end

        T = diag(lmd, 1)+diag(miu, -1)+eye(k)*2;
        M = T\D;
        S_val = zeros(size(x_val));
        for i=1:(k-1)
            S_tmp = M(i)*(x(i+1)-x_val).^3/6/h(i)+M(i+1)*(x_val-x(i)).^3/6/h(i)...
                +(y(i)-M(i)*h(i)^2/6)*(x(i+1)-x_val)/h(i)+(y(i+1)-M(i+1)*h(i)^2/6)*(x_val-x(i))/h(i);
            if i==k-1
                S_val = S_val+S_tmp.*(x(i)<=x_val & x_val<=x(i+1));
            else
                S_val = S_val+S_tmp.*(x(i)<=x_val & x_val<x(i+1));
            end
        end
    end

    function L_val = lpline(x, y, x_val)
        k = length(x);
        d = length(x_val);
        
        x = x(:);
        y = y(:);
        x_val = x_val(:);
         
        a = ones(1, k-1);
        b = ones(d, 1);
        
        L_val = 0;
        for i=1:k
            x_tmp = x([1:i-1 i+1:k]);
            L_val = L_val+y(i)*prod((x_val*a-b*x_tmp')./(b*(x(i)*a-x_tmp')), 2);
        end
    end

    function P_val = pline(x, y, x_val)
        v = vander(x);
        a = v\y(:);
        P_val = polyval(a, x_val);
    end

end






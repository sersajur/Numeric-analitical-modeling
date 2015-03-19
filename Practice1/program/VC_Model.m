function [] = VC_Model(filename, years, type)
    if nargin < 3
        type = 'HIP';
    end
    if nargin < 2
        years = 1;
    end
    
    function P = model(a)
        if strcmp(type,'HIP')
            P = a(1) * I + a(2) * H - a(3) * sqrt(I .* H) + a(4);
        else
            P = a(1) * I + a(2) * H + a(3) * G +...
                a(4) * sqrt(I .* H) + a(5) * sqrt(I .* G) + a(6) * sqrt(G .* H) +...
                a(7) * nthroot(I.*H.*G, 3) + a(8);
        end
        %P = a(1) * I.^2 + a(2) * H.^2 + a(3) * I .* H + a(4) * I + a(5) *H + a(6);
    end

    function [a, value] = fit_parameters(P)
        f = @(a) norm(model(a) - P);
        if strcmp(type,'HIP')
            start_point = [0 0 0 0];
            lb = [0 0 -Inf -Inf];
            ub = [1 1 Inf Inf];
        else
            start_point = [0 0 0 0 0 0 0 0];
            lb = [0 0 0 -Inf -Inf -Inf -Inf -Inf];
            ub = [1 1 1 Inf Inf Inf Inf Inf];
        end
        options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
        [a, value] = fmincon(f, start_point,[],[],[],[],lb,ub,[],options);
    end

    function P = p(B, cost, N, Ef)
        P = B ./ N * Ef / cost;
    end

    function [] = draw(P, I, H, color)
        S = 10;
        scatter3(I, H, P, S, color);
        xlabel('Influenza, I');
        ylabel('Bronchitis, H');
        zlabel('Pneumonia, P');
    end

    function [I, P, H, G] = readFile(filename, h)
        data = dlmread(filename);
        I = data(:,1:h); I = I(:);
        P = data(:,h+1:2*h); P = P(:);
        H = data(:,2*h+1:3*h); H = H(:);
        if ~strcmp(type,'HIP')
            G = data(:,3*h+1:4*h); G = G(:);
        else
            G = 0;
        end
    end

    function P = modelVC(a, B, N, cost, Ef)
        gamma = p(B, cost, N, Ef);
        if strcmp(type,'HIP')
            P = a(1) .* I .* (1-gamma) - a(3) .* sqrt(I .* H) * (1-sqrt(gamma));
        else
            P = (1 - gamma) * a(1) * I +...
                (sqrt(gamma) - 1) * (a(4) * sqrt(I.*H) + a(5) * sqrt(I.*G)) +...
                (1 - nthroot(gamma, 3)) * a(7) * nthroot(I.*H.*G, 3);
        end
    end

    % read data from file
    [I, P, H, G] = readFile(filename, years);
    % feat parameters
    a = fit_parameters(P);
    % compute error
    P_exp = model(a);
    eps = abs(P - P_exp);
    % display results
    display(['Model parameters a = ' num2str(a)]);
    display(['Maximum error = ' num2str(max(eps))]);
    display(['Mean error = ' num2str(mean(eps))]);
    % vaccination
    P_vc = modelVC(a, 1, 2, 1, 0.85);
    z = P_vc ./ P;
    display('People saved:');
    disp(z);
    display(['Maximum people saved = ' num2str(max(z))]);
    display(['Mean people saved = ' num2str(mean(z))]);
end


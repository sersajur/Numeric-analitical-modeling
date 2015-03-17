function [] = VC_Model()

    function P = model(a)
        P = a(1) * I + a(2) * H - a(3) * sqrt(I .* H) + a(4);
    end

    function a = fit_parameters(P, I, H)
        f = @(a) sum((model(a) - P) .^ 2);
        options = optimset('MaxFunEvals', 10000);
        [a,minv] = fminsearch(f, [0 0 0 0], options);
    end

    function P = p(B, cost, N, Ef)
        P = B / cost ./ N * Ef;
    end

    function [] = draw(P, I, H, color)
        S = 10;
        scatter3(I, H, P, S, color);
        xlabel('Influenza');
        ylabel('H');
        zlabel('Pneumonia');
    end

    function data_vec = readFile(filename)
        data_vec = dlmread(filename);
        data_vec = data_vec(:);
    end

    function P = modelVC(a, B, cost, Ef)
        gamma = p(B, cost, N, Ef);
        P = a(1) .* I .* (1-gamma) - a(3) .* sqrt(I .* H) * (1-sqrt(gamma));
    end
    
    I = readFile('I.csv');
    H = readFile('H.csv');
    P = readFile('P.csv');    

    a = fit_parameters(P,I,H);
    P_exp = model(a);
    draw(P, I, H, 'red');
    hold on
    draw(P_exp, I, H, 'green');
    hold off
    eps = abs(P-P_exp);
    display(max(eps));
    display(mean(eps));
end


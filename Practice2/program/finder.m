function [y u] = finder(nets, x, k)
    u = zeros(length(nets),k);
    y = zeros(length(nets),1);
    for i = 1:length(nets)
        fprintf('Processing net %d...\n', i);
        net = nets{i};
        f = @(u) abs(net([x u]'));
        [u(i,:), y(i)] = fmincon(f,zeros(1,k),[],[],[],[],...
        [0  0   0  0   0  0   0  0   0 0   0],...
        [20 Inf 20 Inf 20 Inf 20 Inf 1 Inf Inf]);
        fprintf('End processing net %d.\n', i);
    end
    [y,I] = min(y);
    u = u(I,:);
end
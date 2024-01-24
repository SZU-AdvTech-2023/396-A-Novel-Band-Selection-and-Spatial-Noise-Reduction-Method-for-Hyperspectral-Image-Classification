function [order_band] = NMFW_RANK(data)
    % data is an N*B matrix, where N is the number of pixels and B is the number of bands
    % order_band is the sorting result based on information content

    % Calculate information content
    meanVector = mean(data);
    covarianceMatrix = cov(data);
    inverseCovarianceMatrix = inv(covarianceMatrix);
    k = 1 ./ sum((((data - meanVector) * inverseCovarianceMatrix) .* (data - meanVector)), 2);  % Eq7
    weights = k' .* (inverseCovarianceMatrix * (data - meanVector)');  % Eq1

    % Calculate absolute values of weights
    weights = abs(weights);

    % Calculate rho (mean of absolute weights)
    rho = mean(weights, 2);  % Eq10

    % Sort bands based on rho in descending order
    [~, order_band] = sort(rho, 'descend');
    % [~, order_band] = sort(rho);  % Uncomment this line if you want ascending order
end

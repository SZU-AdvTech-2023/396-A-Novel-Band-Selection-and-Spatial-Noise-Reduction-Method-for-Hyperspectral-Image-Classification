function [data, M, m] = scaleData(data, M, m)
    % Scale input data between -1 and 1.
    %
    % INPUT
    %   data: Input data to rescale
    %   M: Maximum value of the output data
    %   m: Minimum value of the output data
    %
    % OUTPUT
    %   data: Rescaled data
    [Nb_s, Nb_b] = size(data);

    if nargin == 1
        M = max(data, [], 1);
        m = min(data, [], 1);
    end

    data = 2 * (data - repmat(m, Nb_s, 1)) ./ repmat(M - m, Nb_s, 1) - 1;
end

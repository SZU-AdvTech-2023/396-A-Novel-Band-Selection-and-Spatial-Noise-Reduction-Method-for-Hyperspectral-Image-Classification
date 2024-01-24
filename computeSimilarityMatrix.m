function similarityMatrix = computeSimilarityMatrix(X)
    % Normalization, which was not mentioned in the original paper
    % for i = 1 : L
    %     X(:, i) = X(:, i) / norm(X(:, i));
    % end

    % Compute squared Euclidean distance using L2_distance function
    similarityMatrix = sqrt(L2_distance(X, X));
end

% Compute squared Euclidean distance
% ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
function distanceMatrix = L2_distance(a, b)
    % a, b: two matrices, each column is a data
    % distanceMatrix: distance matrix of a and b
    if (size(a, 1) == 1)
        a = [a; zeros(1, size(a, 2))];
        b = [b; zeros(1, size(b, 2))];
    end

    aa = sum(a .* a);
    bb = sum(b .* b);
    ab = a' * b;
    distanceMatrix = repmat(aa', [1, size(bb, 2)]) + repmat(bb, [size(aa, 2), 1]) - 2 * ab;

    distanceMatrix = real(distanceMatrix);
    distanceMatrix = max(distanceMatrix, 0);
end

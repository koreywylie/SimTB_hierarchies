
function r = drchrnd(a, n)
%% Dirichlet random var. w/ param. vector alpha, size n
p = length(a);
r = gamrnd(repmat(a, n, 1), 1, n, p); % gamma random var., from stats
r = r ./ repmat(sum(r, 2), 1, p);
end
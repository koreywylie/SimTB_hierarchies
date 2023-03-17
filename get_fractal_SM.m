%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SM = get_fractal_SM(V, D, L, equal_size, verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates matrix of binary spatial maps,
%  w/ self-similar fractal structure,
%    by parcellating V total voxels,
%    into D groups,
%    over L total levels.
%  If equal_size is false,
%   group size drawn from a Dirichlet pdf w/ alpha=(D,D,..D)
%
% 3/17/2023 --kw--
%   -written by Korey Wylie, copyleft,
%      https://github.com/koreywylie/SimTB_hierarchies
%

if nargin < 2
  D = 2;
end
D = max(1,D);

Lmax = 1;
if D > 1
  Lmax = floor(log(V) / log(D));  % max. number of levels
end
if nargin >= 3
  L = min(Lmax, L);
else
  L = Lmax;
end
if nargin < 4
  equal_size = true;
  verbose = false;
end
if nargin < 5
  verbose = false;
end
if verbose
  fprintf('\nParcellating V = %d voxels, with degree D = %d (i.e., d-way split) over L = %d levels',V,D,L)
  if ~equal_size
    fprintf('\n%s', '...with randomly-sized blocks at each level, drawn from a Dirichlet distribution')
  end
end

K = sum(D.^(1:L));
SM = zeros(K, V);

l = 1;  % level index
k = 0;  % row index
kl_inds = [];  % running index of level indices for rows
D0 = 0; % index of initial row in prev. level
VB = V;  % voxels in blocks, total length of all blocks

if equal_size
  B_d = ones(1, D) * floor(VB / D);  % lengths of current evenly-sized blocks
else
  B = floor(VB / D);
  % B_d = drchrnd(ones(1, D), 1) * VB;  % alpha param. constant, roughly uniform dist. of block sizes at all levels
  B_d = drchrnd(ones(1, D) * D, 1) * VB;  % alpha param. proportional to degree at all levels
  % B_d = drchrnd(ones(1, D) * B, 1) * VB;  % alpha param. proportional to block size, increasing spread for smaller blocks
  B_d = round(B_d);
end

R = VB - sum(B_d);  % sum will be off by at most +/-D
if R < 0  % shrink all blocks by excess voxel
  B_d = B_d - 1;
  B_d = max(B_d, 0);
elseif R > D  % buffer all blocks w/ extra voxel
  B_d = B_d + 1;  
end
if VB ~= sum(B_d)  % if needed, pad early blocks w/ extra voxel
  R = VB - sum(B_d);
  for d=1:R
    B_d(d) = B_d(d) + 1;
  end
end

b1 = 0;  % start indice of block within level
b2 = 0;  % final indice of block within level
for d=1:D
  k = k + 1;
  kl_inds = [kl_inds, l];
  b1 = b2 + 1;
  b2 = b2 + B_d(d);
  SM(k, b1:b2) = ones(1, B_d(d));
end


if L > 1
  for l=2:L
    b1 = 0;
    b2 = 0;
    DLprev = sum(kl_inds == (l - 1));  % number of rows for prev. level
    D0 = sum(kl_inds < (l - 1));  % index of initial row in  prev. level
    
    for d0=1:DLprev
      D0 = D0 + 1;  % index of prev. block split in current level
      VB = sum(SM(D0,:) == 1);  % length of prev. block to be parcellated
      
      if equal_size
        B_d = ones(1, D) * floor(VB / D);  % lengths of current evenly-sized blocks
      else
        B = floor(VB / D);
        % B_d = drchrnd(ones(1, D), 1) * VB;  % alpha param. constant, roughly uniform dist. of block sizes at all levels
        B_d = drchrnd(ones(1, D) * D, 1) * VB;  % alpha param. proportional to degree at all levels
        % B_d = drchrnd(ones(1, D) * B, 1) * VB;  % alpha param. proportional to block size, increasing spread for smaller blocks
        B_d = round(B_d);
      end
      
      R = VB - sum(B_d);   % remainder of current block to distribute
      if R < 0
        B_d = B_d - 1;
        B_d = max(B_d, 0);
      elseif R > D
        B_d = B_d + 1;
      end
      if VB ~= sum(B_d)  % if needed, pad early blocks w/ extra voxel
        R = VB - sum(B_d);
        for d=1:R
          B_d(d) = B_d(d) + 1;
        end
      end
      
      for d=1:D
        k = k + 1;
        kl_inds = [kl_inds, l];
        b1 = b2 + 1;
        b2 = b2 + B_d(d);
        SM(k, b1:b2) = ones(1, B_d(d));
      end
    end
  end
end

% if ismember(0, sum(SM,2))
%   k0 = (0 ~= sum(SM, 2));
%   SM = SM(k0,:);
% end

if verbose
  fprintf('\n...resulting in %d total components', size(SM,1))
end
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = drchrnd(a, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet random var. w/ param. vector alpha, sample size n %
  p = length(a);  % prior sample sizes for Bayesian analysis
  r = gamrnd(repmat(a, n, 1), 1, n, p); % gamma random var., from stats package
  r = r ./ repmat(sum(r, 2), 1, p); % normalize sum of groups to 1
end

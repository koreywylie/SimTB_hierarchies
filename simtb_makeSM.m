
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SM = simtb_makeSM(sP, sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates component Spatial Maps (SMs) as nested hierarchy of networks,
%   for use with SimTB fMRI simulation toolbox
% Replaces fn. of same name in SimTB package,
%   must receive priority in matlab path
% Requires:
%   get_fractal_SM()
% Recommends:
%   simtb_createmask() modified version of fn., creates inclusive mask of all voxels
%

% 3/17/2023 --kw-- 
%   Instructions:
%      Script replaces original simtb_makeSM.m in simtb,
%        must receive priority over above in matlab path.

nC = sP.nC;
nV = sP.nV;
% SM_source_ID = sP.SM_source_ID;
% SM_translate_x = sP.SM_translate_x;
% SM_translate_y = sP.SM_translate_y;
% SM_theta = sP.SM_theta;
SM_present = sP.SM_present;
% SM_spread  = sP.SM_spread;
% %% initialize SM
% SM = zeros(nC, nV*nV);

if isfield(sP, 'HIER_groupSM_file') && ~isempty(sP.HIER_groupSM_file) && exist(sP.HIER_groupSM_file, 'file')
  load(sP.HIER_groupSM_file, 'SM');  % use shared group spatial map
  nC = size(SM, 1);
else
  D = 2;
  if isfield(sP, 'HIER_degree') && ~isempty(sP.HIER_degree)
    D = sP.HIER_degree;  % number of parcellations/subdivisions at each level of hierarchy
  end
  L = 3;
  if isfield(sP, 'HIER_levels') && ~isempty(sP.HIER_levels)
    L = sP.HIER_levels;     % number of levels in hierarchy
  end
  eq_size = true;
  if isfield(sP, 'HIER_equal_size') && ~isempty(sP.HIER_equal_size)
    eq_size = sP.HIER_equal_size;  % equal sized-parcellations if true, randomly-sized elsewise
  end
  
  SM = get_fractal_SM(nV*nV, D, L, eq_size);
  
  if isfield(sP, 'HIER_groupSM_file') && ~isempty(sP.HIER_groupSM_file)
    save(sP.HIER_groupSM_file, 'SM');
  end
end


mask = simtb_createmask(sP);  % 10/6/2021 --kw-- Recommend modified version of fn., returning matrix of all ones
mask = reshape(mask,1,nV*nV);

for c=1:nC
  Temp = SM(c,:);
  SM(c,:) = mask.*(SM_present(sub,c)*reshape(Temp,1,nV*nV) + 0.005*randn(1, nV*nV));
  % Temp =  simtb_generateSM(SM_source_ID(c), nV, SM_translate_x(sub,c), SM_translate_y(sub,c), SM_theta(sub,c), SM_spread(sub,c));
  %SM(c,:) = mask.*(SM_present(sub,c)*reshape(Temp,1,nV*nV));
  % add a bit of Gaussian noise so SMs will never be flat or completely identical
  % SM(c,:) = mask.*(SM_present(sub,c)*reshape(Temp,1,nV*nV) + 0.005*randn(1, nV*nV));    
  clear Temp     
end

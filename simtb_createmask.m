%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MASK = simtb_createmask(sP, saveMASKasNIFTI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates inclusive mask of all voxels,
%   to simplify spatial structure when creating hierarchy of networks.
%
% 3/17/2023 --kw-- 
%   Instructions:
%      Script replaces original simtb_createmask.m in simtb,
%        must receive priority over above in matlab path.


if nargin < 2
    saveMASKasNIFTI = 0;
end

nV = sP.nV;
% arg1 = linspace(-1,1,nV);
% [x,y] = meshgrid(arg1,arg1);
% r = sqrt(x.^2 + y.^2);

MASK = ones(nV,nV);
% MASK(r>1) = 0;

if saveMASKasNIFTI
    if exist('spm.m', 'file') 
        % reshape the data
        MASK = reshape(MASK, nV,nV,1,1);
        mfilename = simtb_makefilename(sP, 'MASK');
        if sP.D_motion_FLAG && (any(sP.D_motion_deviates(:)>0))
            % motion has been implemented for at least one subject
            warning('Mask is not appropriate if motion has been implemented')
        end
        fprintf('Writing MASK to %s\n', mfilename);
        simtb_saveasnii(MASK, mfilename);
    else
        fprintf('In order to save in nifti format, SPM must be on the search path.\n');
    end
end

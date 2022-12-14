
function InfoInitialPoint = InitialPoint(ReadInfo, point, coeff)
% InfoInitialPoint = InitialPoint(ReadInfo, point, coeff)
% 
% DESCRIPTION
%
%
% REQUIRED INPUTS
%   ReadInfo: an structure generated by 'ReadLayers' function
%
% OPTIONAL INPUTS
%   point: 
%   coeff: 
%
% OUTPUTS
%   InfoInitialPoint: a structure containing:
%       -idx: 
%       -SortNormDistance:
%       -coeff:    
%%
    Rows = ReadInfo.Dimensions(1);
    NumLayers = ReadInfo.Dimensions(2);
    NormalizedClimVar = ReadInfo.NormalizedClimVar;
    Distance = zeros(1, Rows);
    
    if nargin < 2
        point = rand(NumLayers, 1);
    end
    if nargin < 3
        coeff = rand(NumLayers, 1);
    end

    coeff = coeff/sum(coeff);
    point = point.*coeff;
    
    NormalizedClimVar = NormalizedClimVar.*coeff;     
    
    for i = 1 : Rows
        Distance(i) = norm(point - NormalizedClimVar(:, i))...
                      * (2 - corr2(point, NormalizedClimVar(:, i)));
    end
    
    NormDistance = 1 - normalize(Distance, 2, 'range');
    [SortNormDistance, idx] = sort(NormDistance, 2, 'descend');
    
    % OUTPUT STORAGE    
    InfoInitialPoint.idx = idx;
    InfoInitialPoint.SortNormDistance = SortNormDistance;
    InfoInitialPoint.coeff = coeff;

end
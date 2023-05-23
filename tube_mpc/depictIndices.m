%% depictIndices
%
% Function depictIndices plot parition and depicts indices of the selected
% critical region.
%
% depictIndices( ExplicitController, indices2depict, figureHandle )
%
% INPUTS:
%   ExplicitController  - mandatory: explicit controller eith constructed
%                         parition (explicit spultion map)
%   indices2depict      - optional: indices (vector of integers) to be depicted
%   figureHandle        - optional: figure handle
%
%   Example #1:
%   depictIndices(ExplicitController)
%
%   Example #2:
%   depictIndices(ExplicitController, [1 : 3 : 9])
%
%   Example #3:
%   figure, hold on
%   eMPC.partition.Set(1).plot
%   eMPC.partition.Set(3).plot
%   eMPC.partition.Set(5).plot
%   depictIndices(ExplicitController, [1,3,5], gcf)
function [] = depictIndices( obj, indices2depict, figureHandle )

% Check input Explicit Controller
if ( obj.isExplicit == 0 )
    error('MPTplus: explicit controller required to depict indices of its partition!')
end

% Problem size
Nr = obj.nr;
Nx = obj.nx;

% Check indices to be depicted
if ( nargin == 1 )
    indices2depict = [ 1 : Nr];
    % If no "figureHandle" is set, then assign current figure by "gcf"
    figureHandle = gcf;

elseif( nargin == 2 )
    % Check input set of indices to depict - vector of indeices is equired
    if( isequal(class(indices2depict),'double') == 0 )
        error(sprintf('MPTplus: depictIndices: indices2depict must be vector of integers!'))
    end
    % If input "indices2depict" is empty, then assign all indices
    if ( isempty(indices2depict) == 1 )
        indices2depict = [ 1 : Nr];
    end
    % If no "figureHandle" is set, then assign current figure by "gcf"
    figureHandle = gcf;

elseif( nargin == 3 )
    % Check "figureHandle"
    if ( isequal(class(figureHandle),'matlab.ui.Figure') == 0 )
        error(sprintf('MPTplus: depictIndices: figureHandle must be a figure handle (e.g.: "gcf")!'))
    end

else
    error(sprintf('MPTplus: depictIndices: Unexpected number of inputs: %d!',nargin))
end

% Call figure
figure(figureHandle)

% Depict indices of critical regions
for r = indices2depict
    xCheby = obj.partition.Set(r).chebyCenter;
    xChebyCenter = xCheby.x;
    if( Nx == 1 )
        text(xChebyCenter(1), 0, num2str(r) );
    elseif( Nx == 2 )
        text(xChebyCenter(1), xChebyCenter(2), num2str(r) );
    elseif( Nx == 2 )
        text(xChebyCenter(1), xChebyCenter(2), xChebyCenter(3), num2str(r) );
    else
        error(sprintf('MPTplus: Unable to depict indices for %d-D partitions!',obj.nx))
    end
end
end % function
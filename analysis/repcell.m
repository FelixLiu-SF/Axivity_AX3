function [cellout] = repcell(sizein, reparray)
% [cellout] = repcell(sizein, reparray)

if(iscell(reparray))
    
    cellout = cell(sizein);
    for ix = 1:size(cellout,1)
        for jx = 1:size(cellout,2)
            cellout{ix,jx} = reparray;
        end
    end
    
    
else

    tmpmat = repmat(reparray,sizein);

    cellout = mat2cell(tmpmat,ones(sizein),size(reparray,2));
    
end


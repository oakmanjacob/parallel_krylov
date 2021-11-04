function [L,U,perm,rperm,A] = my_ILUTP(A,droptol,lfil)
% Input:        A:          nxn MATLAB sparse matrix
%               dropTol:    Threshold for dropping
%               lfil:       Maximum number of nonzero elements allowed in
%                           each row of L and U
% Output:       L:          MATLAB sparse, lower triangular matrix
%               U:          MATLAB sparse, upper cleartriangular matrix
%               perm:       Permutation array representing the old
%                           locations of columns
%               rperm:      Inverse of the permutation array representing
%                           the new locations of columns (note: this array
%                           is necessary since the algorithm works in the 
%                           new ordering) 
%               A:          MATLAB sparse, column permuted A
% Modified from routine in SPARSKIT, Yousef Saad Sep 8, 1993 - Latest 
% revision, August 1996
% Author: Arielle (Grim-McNally) Car, Eric de Sturler 2013
% Latest revision: Arielle Carr June 2017
n = size(A,1);

% Set protection against zero pivots
pivTol = 1.e-12;
pivMod = 1.e-3;

perm = 1:n; 
rperm = 1:n;

[~,~,valA] = find(A');
nnzA = length(valA);
LU_nnz = nnzA+n + n*lfil;
At = transpose(A);

% Allocate space for MSR storage of LU (as a single matrix)
LUval = zeros(1,LU_nnz); % Stores the values of LU (LUval(1:n) 
                         % stores the diagonal elements of the matrix LU, 
                         % LUval(n+2:end) stores (by row) the off-diagonal 
                         % elements of LU)
LUindex = zeros(1,LU_nnz); % Stores position information for values in LU 
                           % (LUindex(n+2:end) stores the column indices
                           % of the off-diaganal elements; LUindex(1:n+1)
                           % stores pointers to the start of each row 
                           % for the off-diagonals)

% Pointers to the beginning of each row of U in LUval and LUindex
Urows = zeros(n,1);

Ustart = n + 2;
LUindex(1) = Ustart;

for row = 1:n
    
    % Initialize work arrays for current row - Refer to "Iterative Methods 
    % for Sparse Linear Systems 2nd Ed" (Yousef Saad 2003), Chapter 10.4.3 
    % for a more detailed explanataion of these data structures.
    w = zeros(1,n); % Contains the nonzero values of the current row 
    jw = zeros(1,n); % Points to the column numbers of the nonzero values   
                     % in the current row 
    rjw = zeros(1,n); % Vector of length n, which is nonzero only in 
                      % positions there is a nonzero value 
                      % in the current row and is zero elsewhere
                      % (This will help to determine if an 
                      % element is fill-in or not when performing the 
                      % elimination)
    
    % Compute the norm of the current row
    currentRow = At(:,row);
    dropNorm = norm(currentRow);
    if dropNorm == 0
        fprintf('Matrix is singular. (Zero row) \n');
        dropNorm = 1.e-6;
    end
    
    if mod(row,1000) == 0
        row
    end
    
    % Determine upper (right of diagonal) and lower (left of diagonal)
    % portions of current row
    [cols,~,vals] = find(currentRow);
    newcol = rperm(cols);
    diagcoef = perm(row);
    leftrow = find(newcol < row);
    lenL = length(leftrow);
    rightrow = find(newcol > row);
    lenU = length(rightrow) + 1;
    
    % Store diagonal coefficient,location in work arrays
    w(row) = A(row,diagcoef);
    jw(row) = row;
    rjw(row) = row;
    
    % Store lower part (all values to the left of diagonal coefficient) in
    % work arrays
    if lenL > 0
        jw(1:lenL) = newcol(leftrow);
        rjw(newcol(leftrow)) = 1:lenL;
        w(1:lenL) = vals(leftrow);
    end
    
    % Store upper part (all values to the right of diagonal coefficient) in
    % work arrays
    if lenU > 1
        jw(row+1:row+lenU-1) = newcol(rightrow);
        rjw(newcol(rightrow)) = row+1:row+lenU-1;
        w(row+1:row+lenU-1) = vals(rightrow);
    end
    
    % Eliminate previous rows
    nextw = 1;
    len = 0;
    while nextw <= lenL
        % Find the smallest column index among jw(nextw:lenl)
        minFind = jw(nextw:lenL);
        [mincol,idx] = min(minFind);
        if idx > 1 % If the smallest column is not (immediately) next,
                   % exchange w, rjw, jw 
            nextw_col = jw(nextw);
            jw(nextw) = jw(nextw+idx-1);
            jw(nextw+idx-1) = nextw_col;
            
            rjw(mincol) = nextw;
            rjw(nextw_col) = nextw+idx-1;
            
            nextw_val = w(nextw);
            w(nextw) = w(nextw+idx-1);
            w(nextw+idx-1) = nextw_val;
        end        
        nextCol = mincol;
        
        % Zero element out
        rjw(nextCol) = 0;
        
        multiplier = w(nextw)/LUval(nextCol);
        if abs(multiplier) > droptol % Apply dropping to L
            for k = Urows(nextCol):LUindex(nextCol+1)-1;
                MultALU = LUval(k);
                subtract = multiplier*MultALU;
                loc = rperm(LUindex(k)); % Determine the new column index
                fill = rjw(loc);
                
                if loc >= row % Update upper part 
                    if fill == 0 % The element is fill-in
                        lenU = lenU + 1;
                        i = row + lenU - 1;
                        jw(i) = loc;
                        rjw(loc) = i;
                        w(i) = -subtract;
                    else % It is not fill-in
                        w(fill) = w(fill) - subtract;
                    end
                elseif loc < row % Update the lower part
                    if fill == 0 % The element is fill-in
                        lenL = lenL + 1;
                        jw(lenL) = loc;
                        rjw(loc) = lenL;
                        w(lenL) = -subtract;
                    else % It is not fill-in
                        w(fill) = w(fill) - subtract;
                    end
                end
            end
            
            % Store pivot
            len = len + 1;
            w(len) = multiplier;
            jw(len) = nextCol;
            
        end
        
        nextw = nextw + 1;
        
    end
    
    % Enforce maximum fill on L
    lenL = len;
    len = min(lenL, lfil-1);
    absW = abs(w(1:lenL));
    holdW = w(1:lenL);
    [~,idxL] = sort(absW,'descend');
    w(1:lenL) = holdW(idxL);
    %     [w(1:lenL),idxL] = sort(w(1:lenL),'descend');
    jw(1:lenL) = jw(idxL);
    
    % Store elements of L (in old ordering)
    LUval(Ustart:((Ustart+len)-1)) = w(1:len);
    LUindex(Ustart:((Ustart+len)-1)) = perm(jw(1:len));
    
    % Save pointer to beginning of U(row,:)
    Ustart = Ustart+len;
    Urows(row) = Ustart;
    
    % Apply drop tolerance to U
    len = 0;
    for k = 1:lenU-1
        if abs(w(row+k)) >= droptol*dropNorm
            len = len + 1;
            w(row+len) = w(row+k);
            jw(row+len) = jw(row+k);
        end
    end
    
    % Enforce maximum fill on U
    lenU = len + 1;
    len = min(lenU,lfil);
    uppIdx = lenU-1;
    if len <= uppIdx
        absW = abs(w(row+1:(row+uppIdx)));
        holdW = w(row+1:(row+uppIdx));
        [x,idxU] = sort(absW,'descend');
        %         [x,idxU] = sort(w(row+1:(row+uppIdx)),'descend');
        jwHold = jw(row+1:(row+lenU-1));
        jx = jwHold(idxU);
        %         w(row+1:row+uppIdx) = x;
        w(row+1:row+uppIdx) = holdW(idxU);
        jw(row+1:row+uppIdx) = jx;
    end

%     % Update L
%     lenL = len;
%     len = min(lenL, lfil-1);
%     [w(1:lenL),jw(1:lenL)] = qSplit(w(1:lenL), jw(1:lenL), lenL, len);
%     %[w,jw] = qSplit(w, jw, len, lenl);
%     % Store L part (in original coordinates) in modified and compressed
%     for k = 1:len
%         LUval(Ustart) = w(k);
%         LUindex(Ustart) = perm(jw(k));
%         Ustart = Ustart + 1;
%     end
%     
%     % Save pointer to beginning of row ii of U
%     Urows(row) = Ustart;
%     
%     % Update U matrix (first apply dropping strategy)
%     len = 0;
%     for k = 1:lenU-1
%         if abs(w(row+k)) >= droptol*dropNorm
%             len = len + 1;
%             w(row+len) = w(row+k);
%             jw(row+len) = jw(row+k);
%         end
%     end
%     
%     lenU = len + 1;
%     len = min(lenU,lfil);
%     [x, jx] = qSplit(w(row+1:(row+lenU-1)),jw(row+1:(row+lenU-1)), lenU - 1, len);
%     w(row+1:row+lenU-1) = x;
%     jw(row+1:row+lenU-1) = jx;
    
    % Determine the next pivot
    pivIdx = row;
    piv = abs(w(pivIdx));
    for k = row+1:row+lenU-1
        pivComp = abs(w(k));
        if pivComp > piv
            pivIdx = k;
            piv = pivComp;
        end
    end
    
    % Exchange w, perm, rperm:
    tmp = w(row);
    w(row) = w(pivIdx);
    w(pivIdx) = tmp;
    
    loc = jw(pivIdx);
    i = perm(row);
    perm(row) = perm(loc);
    perm(loc) = i;
    
    rperm(perm(row)) = row;
    rperm(perm(loc)) = loc;
    
    % Store elements of U (in old ordering)
    Uend = length(row+1:row+len-1);
    if Uend ~= 0
        LUindex(Ustart:(Ustart+Uend-1)) = perm(jw(row+1:row+len-1));
        LUval(Ustart:Ustart+Uend-1) = w(row+1:row+len-1);
    end
    
    
    % Store diagonal element of U 
    if abs(w(row)) < pivTol
        w(row) = pivMod;
    end
    LUval(row) = w(row);
    
    % Point to the next row in U
    Ustart = Ustart+Uend;
    LUindex(row+1) = Ustart;
end

% Permute columns of LU (using new ordering)
LUindex(LUindex(1): LUindex(n+1) - 1) = ...
    rperm(LUindex(LUindex(1): LUindex(n+1) - 1));

% Convert to Compressed Sparse Row (from SPARSKIT function msrcsr.f, last 
% modified May 29, 1994 by R. Bramley)
LUvalHold = LUval;
rowIdx(1) = 1; % Row indices (for CSR format)
idx = 1;
for row = 1:n
    addDiag = 0;
    idxDiag = idx + (LUindex(row+1) - LUindex(row));
    for k = LUindex(row):LUindex(row+1)-1
        loc = LUindex(k);
        if loc < row || addDiag
            val(idx) = LUval(k);
            colIdx(idx) = loc;
            idx = idx + 1;
        else
            idxDiag = idx;
            idx = idx + 1;
            addDiag = 1;
            val(idx) = LUval(k);
            % Check for zero elements along the diagonal
            if val(idx) == 0
                fprintf('Matrix is singular. (Zero row) \n');
            end
            colIdx(idx) = loc;
            idx = idx + 1;
        end
    end
    val(idxDiag) = LUvalHold(row);
    colIdx(idxDiag) = row;
    if ~addDiag
        idx = idx + 1;
    end
    rowIdx(row+1) = idx;    
end

% Determine row indices (for COO format to convert to sparse matrix)
for i = 1:n
    rowNewIdx(rowIdx(i):rowIdx(i+1)-1) = i;
end

% Permute columns of A
A = A(:,perm(:));

% Convert LU to Matlab sparse matrix
LU = sparse(rowNewIdx,colIdx,val);

% Separate L and U using tril, triu so that Matlab recognizes them as 
% lower and upper, respectively, triangular matrices
L = speye(n) + tril(LU,-1);
U = triu(LU);
end

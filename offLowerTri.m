function lowTri=offLowerTri(theMat)
%%% returns the concatenated elements in the upper triangular
lowTri=theMat(triu(ones(size(theMat)))==0);
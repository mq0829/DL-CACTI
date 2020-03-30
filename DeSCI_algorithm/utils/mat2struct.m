function [sdata] = mat2struct(mdata)
%MAT2STRUCT Convert a matrix to a structure along the last dimension of the
%matrix.
nlast = size(mdata,ndims(mdata)); % size of the last dimension
odims = repmat({':'},1,ndims(mdata)-1);
sdata(1:nlast) = struct('cdata',zeros(size(mdata(odims{:},1)),'like',mdata));
for ilast = 1:nlast
    sdata(ilast).cdata = mdata(odims{:},ilast);
end

end


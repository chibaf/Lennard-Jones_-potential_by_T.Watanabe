function retval = readLJ(varargin)
datfile = varargin{1};
try
    fid = fopen(datfile,'r');
catch
    error('data file not found!');
end
i=0;
retval=[];
while feof(fid) == 0
  a = fgetl(fid);
  if (i>0)
    retval=[retval;str2num(a)];
  endif
  i=i+1;
end
fclose(fid);
endfunction

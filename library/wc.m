%  wc -l
%  Number of lines

function count = wc(fname)
fname = strrep(fname,' ','');

% try
    [s,out]=system(['wc -l ',fname]);
    out = strrep(strrep(out,fname,''),' ','');
    count = str2num(out);
%     warning(lasterr)
% catch
%     fh = fopen(fname, 'rt');
%     assert(fh ~= -1, 'Could not read: %s', fname);
%     x = onCleanup(@() fclose(fh));
%     count = 0;
%     while ~feof(fh)
%         count = count + sum( fread( fh, 16384, 'char' ) == char(10) );
%     end
% end

function r3_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end

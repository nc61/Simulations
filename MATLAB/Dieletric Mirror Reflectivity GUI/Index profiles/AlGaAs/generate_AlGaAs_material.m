lambda_m = 1e-6:1e-9:6e-6;
x = [0, 0.32, 0.7];

for ind = 1:length(x)
    n = index2(lambda_m, x(ind));
    fname = "AlGaAs(x="+num2str(x(ind))+").txt";
    
    fid = fopen(fname, 'w+');
    k = 0;
    format long
    for jnd = 1:length(lambda_m)
        line_to_write = sprintf('%s\t%s\t%s\r\n', num2str(lambda_m(jnd)*1e6), n(jnd), num2str(k));
        fwrite(fid, line_to_write);
    end
    fclose(fid);
end

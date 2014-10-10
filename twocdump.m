% twocdump - dump two color results into a CSV file
function twocdump(f,filename)
fd=fopen(filename,'w');
if fd<0
  error('Unable to open file: %s',filename);
end

fprintf(fd,'Index,Desc,N,N(gated),Median(GFP),Median(mCherry),mu,sigma,cilow,cihigh\n');
for i=1:length(f)
  ff=f(i);
  hdr=ff.hdr;
  sel=ff.P(ff.usegatenum,:);
  
  fprintf(fd,'%d,"%s"',i,hdr.cells);
  fprintf(fd,',%d,%d',size(ff.data,1),sum(sel));
  fprintf(fd,',%.0f',median(ff.gfp(sel)));
  fprintf(fd,',%.0f',median(ff.cherry(sel)));
  fprintf(fd,',%.3f',ff.mu);
  fprintf(fd,',%.3f',ff.sigma);
  fprintf(fd,',%.3f,%.3f',ff.muci80);
  fprintf(fd,'\n');
end
fclose(fd);

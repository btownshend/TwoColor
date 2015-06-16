in=fopen('12way-macros.xml','r');
out=fopen('12way-out.xml','w');
minx=30;
maxx=262144;
miny=minx;
maxy=maxx;
%elow=10^(3.82172-5.35478);   % Bin 1-2 edge
%ehigh=10^(5.47689-4.9286);
elow=10^(3.88548-5.41854);   % Bin 1-2 edge
ehigh=10^(5.418540-4.87025);
elow=2;
ehigh=30;
ratios=logspace(log10(elow),log10(ehigh),11);
fprintf('Ratios: %s\n', sprintf('%.2f ',ratios));
while true
  line=fgetl(in);
  if line==-1
    break;
  end
  line=strrep(line,'##Max.X##',sprintf('%f',log10(maxx)));
  line=strrep(line,'##Min.X##',sprintf('%f',log10(minx)));
  line=strrep(line,'##Max.Y##',sprintf('%f',log10(maxy)));
  line=strrep(line,'##Min.Y##',sprintf('%f',log10(miny)));
  for i=1:length(ratios)
    r=ratios(i);
    prefix=sprintf('##E%d.%d',i,i+1);
    if r<1
      line=strrep(line,[prefix,'.low.X##'],sprintf('%f',log10(miny/r)));
      line=strrep(line,[prefix,'.low.Y##'],sprintf('%f',log10(miny)));
      line=strrep(line,[prefix,'.high.X##'],sprintf('%f',log10(maxx)));
      line=strrep(line,[prefix,'.high.Y##'],sprintf('%f',log10(maxx*r)));
    else
      line=strrep(line,[prefix,'.low.X##'],sprintf('%f',log10(minx)));
      line=strrep(line,[prefix,'.low.Y##'],sprintf('%f',log10(minx*r)));
      line=strrep(line,[prefix,'.high.X##'],sprintf('%f',log10(maxy/r)));
      line=strrep(line,[prefix,'.high.Y##'],sprintf('%f',log10(maxy)));
    end
  end
  fprintf(out,'%s\n',line);
end
fclose(in);
fclose(out);

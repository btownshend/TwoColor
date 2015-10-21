% Create a BD FACS settings file with a set of ratio gates setup
% This will create the given number of bins with the lowest and highest bins as
%  triangular regions enclosing everything below or above the given ratio bounds
% The remaining pins are strips (actually 4-side polygons) that are equally space% in terms of log(ratio)
function bdfacs_setratios(nbins,ratiolow,ratiohigh,outfile,varargin)
defaults=struct('minx',30,'miny',30,'maxx',262144,'maxy',262144);
args=processargs(defaults,varargin);

macrofilename=sprintf('%dway-macros.xml',nbins);
in=fopen(macrofilename,'r');
if in<0
  error('Unable to open macro file: %s\n', macrofilename);
end
out=fopen(outfile,'w');
ratios=logspace(log10(ratiolow),log10(ratiohigh),11);
fprintf('Ratios: %s\n', sprintf('%.2f ',ratios));
while true
  line=fgetl(in);
  if line==-1
    break;
  end
  line=strrep(line,'##Max.X##',sprintf('%f',log10(args.maxx)));
  line=strrep(line,'##Min.X##',sprintf('%f',log10(args.minx)));
  line=strrep(line,'##Max.Y##',sprintf('%f',log10(args.maxy)));
  line=strrep(line,'##Min.Y##',sprintf('%f',log10(args.miny)));
  for i=1:length(ratios)
    r=ratios(i);
    prefix=sprintf('##E%d.%d',i,i+1);
    if r<1
      line=strrep(line,[prefix,'.low.X##'],sprintf('%f',log10(args.miny/r)));
      line=strrep(line,[prefix,'.low.Y##'],sprintf('%f',log10(args.miny)));
      line=strrep(line,[prefix,'.high.X##'],sprintf('%f',log10(args.maxx)));
      line=strrep(line,[prefix,'.high.Y##'],sprintf('%f',log10(args.maxx*r)));
    else
      line=strrep(line,[prefix,'.low.X##'],sprintf('%f',log10(args.minx)));
      line=strrep(line,[prefix,'.low.Y##'],sprintf('%f',log10(args.minx*r)));
      line=strrep(line,[prefix,'.high.X##'],sprintf('%f',log10(args.maxy/r)));
      line=strrep(line,[prefix,'.high.Y##'],sprintf('%f',log10(args.maxy)));
    end
  end
  fprintf(out,'%s\n',line);
end
fclose(in);
fclose(out);

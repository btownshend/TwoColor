% Create a BD FACS settings file with a set of ratio gates setup
% This will create the given number of bins with the lowest and highest bins as
%  triangular regions enclosing everything below or above the given ratio bounds
% The remaining pins are strips (actually 4-side polygons) that are equally space% in terms of log(ratio)
function ratios=bdfacs_setpolys(nbins,ratiolow,ratiohigh,outfile,varargin)
defaults=struct('minx',30,'miny',30,'maxx',262144,'maxy',262144,'doplot',true,'gapfrac',.001);
args=processargs(defaults,varargin);

macrofilename=sprintf('bd-%.0fpolys.xml',nbins);
in=fopen(macrofilename,'r');
if in<0
  error('Unable to open macro file: %s\n', macrofilename);
end
out=fopen(outfile,'w');
ratios=logspace(log10(ratiolow),log10(ratiohigh),11);
fprintf('Ratios: %s\n', sprintf('%.2f ',ratios));
% Build polygon coords
% The polygon for bin i is a strip or triangle bounded by the (minx,miny)-(maxx,maxy) rectangle, and the diagonal lines:
%    y=ratio(i)*x (except for the last bin)
%    y=ratio(i-1)*x  (except for bin 1)
polys={};
if args.doplot
  setfig('bdfacs_setpolys'); clf;
  loglog([args.minx,args.minx,args.maxx,args.maxx,args.minx],[args.miny,args.maxy,args.maxy,args.miny,args.miny],':');
  hold on;
end
stretch=10.^-(diff(log10(ratios(1:2)))*args.gapfrac/2);
fprintf('Scaling bin edges in by %fx\n', stretch);
for i=1:nbins
  if i==1
    r=ratios(1)*stretch;
    if args.miny<r*args.minx
      p=[args.minx,args.miny];
      p(end+1,:)=args.minx*[1,r];
    else
      p=args.miny*[1/r,1];
    end
    if args.maxy>r*args.maxx
      p(end+1,:)=args.maxx*[1,r];
    else
      p(end+1,:)=args.maxy*[1/r,1];
      p(end+1,:)=[args.maxx,args.maxy];
    end
    p(end+1,:)=[args.maxx,args.miny];
  elseif i==nbins
    r=ratios(nbins-1)/stretch;
    if args.miny<r*args.minx
      p=args.minx*[1,r];
    else
      p=[args.minx,args.miny];
      p(end+1,:)=args.miny*[1/r,1];
    end
    if args.maxy>r*args.maxx
      p(end+1,:)=args.maxx*[1,r];
      p(end+1,:)=[args.maxx,args.maxy];
    else
      p(end+1,:)=args.maxy*[1/r,1];
    end
    p(end+1,:)=[args.minx,args.maxy];
  else
    r1=ratios(i-1)/stretch;
    r2=ratios(i)*stretch;
    if args.miny<r1*args.minx
      p=args.minx*[1,r1];
    else
      p=args.miny*[1/r1,1];
    end
    if args.maxy>r1*args.maxx
      p(end+1,:)=args.maxx*[1,r1];
    else
      p(end+1,:)=args.maxy*[1/r1,1];
    end
    if (args.maxy>r1*args.maxx) ~= (args.maxy>r2*args.maxx)
      p(end+1,:)=[args.maxx,args.maxy];
    end
    if args.maxy>r2*args.maxx
      p(end+1,:)=args.maxx*[1,r2];
    else
      p(end+1,:)=args.maxy*[1/r2,1];
    end
    if args.miny<r2*args.minx
      p(end+1,:)=args.minx*[1,r2];
    else
      p(end+1,:)=args.miny*[1/r2,1];
    end
    if (args.miny<r1*args.minx) ~= (args.miny<r2*args.minx)
      p(end+1,:)=[args.minx,args.miny];
    end
  end
  str='';
  for j=1:1:size(p,1)
    str=[str,sprintf('<point x="%f" y="%f" />',p(j,:))];
  end
  polys{i}=str;
  if args.doplot
    loglog(p([1:end,1],1),p([1:end,1],2));
    pause(1);
  end
  fprintf('Bin%d: %s\n', i, polys{i});
end

while true
  line=fgetl(in);
  if line==-1
    break;
  end
  for i=1:nbins
    line=strrep(line,sprintf('##Bin%.0f##',i),polys{i});
  end
  fprintf(out,'%s\n',line);
end

fclose(in);
fclose(out);

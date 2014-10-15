% densplot - FACS data density plot
% Brent Townshend 10/2012
% Usage: densplot(x,y,bins,range,dolog)
%  x,y: data
%  bins: bin edges to annotate on plot
%  range: [minx, maxx, miny, maxy]  for plotting
%  dolog: true to use log data
function [z,rng]=densplot(x,y,bins,range,dolog)
if nargin<5
  dolog=0;
end
if isempty(x)
  fprintf('densplot: Nothing to plot\n');
  return;
end
if dolog
  if mean(x<=0 | y<=0)>0.01
    fprintf('densplot: Ignoring %.1f%% of points that are <= 0\n',mean(x<=0 | y<=0)*100);
  end
  sel=x>0 & y>0;
  x=x(sel);
  y=y(sel);
end
if nargin < 4 || length(range)==0
  xs=sort(x);
  ys=sort(y);
  range=[xs(floor(length(xs)*.005)+1),xs(floor(length(xs)*.995)+1),ys(floor(length(ys)*.005)+1),ys(floor(length(ys)*.995)+1)];
end
rng=range;  % Return value
if dolog
  x=log(x);
  y=log(y);
  range=log(range);
end
if nargin<3 || length(bins)==0
  bins=max(100,round(sqrt(length(x)/5)))*[1,1];
end
x=(x-range(1))/(range(2)-range(1));
y=(y-range(3))/(range(4)-range(3));

bx=floor(x*bins(1));
by=floor(y*bins(2));
bx=max(0,min(bins(1)+1,bx))+1;
by=max(0,min(bins(2)+1,by))+1;
z=zeros(bins+3);
for i=1:length(bx)
  z(by(i),bx(i))=z(by(i),bx(i))+1;
end
gx=(([0,1:bins(1)+2]-1)/bins(1))*(range(2)-range(1))+range(1);
gy=(([0,1:bins(2)+2]-1)/bins(2))*(range(4)-range(3))+range(3);
[mx,my]=meshgrid(gx,gy);
if dolog
  pcolor(exp(mx),exp(my),z);
  set(gca,'YScale','log');
  set(gca,'XScale','log');
else
  pcolor(mx,my,z);
end
caxis([0,max(max(z(2:end-2,2:end-2)))]);
shading flat
colorbar




% densplot - FACS data density plot
% Brent Townshend 10/2012
% Usage: densplot(x,y,bins,range,dolog)
%  x,y: data
%  bins: [bx,by] number of bins in each direction
%  range: [minx, maxx, miny, maxy]  for plotting
%  dolog: true to use log data, if it is a 2-element vector, control log setting for x,y independently
function [z,rng]=densplot(x,y,bins,range,dolog,dologicle,scaling)
if nargin<5 || isempty(dolog)
  dolog=[0,0];
end
if nargin<6 || isempty(dologicle)
  dologicle=[0,0];
end

if isempty(x)
  fprintf('densplot: Nothing to plot\n');
  return;
end
if length(dolog)==1
  dolog=[dolog,dolog];
end
if length(dologicle)==1
  dologicle=[dologicle,dologicle];
end
x=x(:);y=y(:);
sel=isfinite(x)&isfinite(y);

if any(dolog&dologicle)
  fprintf('densplot: Both log and logicle specified for the same axis - using logicle');
  dolog=dolog&~dologicle;
end
if dolog(1)
  if mean(x<=0)
    fprintf('densplot: Ignoring %.1f%% of points that have x <= 0\n',mean(x<=0)*100);
  end
  sel=sel & x>0;
end

if dolog(2)
  if mean(y<=0)>0.01
    fprintf('densplot: Ignoring %.1f%% of points that have y <= 0\n',mean(y<=0)*100);
  end
  sel=sel & y>0;
end

x=x(sel);
y=y(sel);

if nargin < 4 || length(range)==0
  xs=sort(x);
  ys=sort(y);
  range=[xs(floor(length(xs)*.005)+1),xs(floor(length(xs)*.995)+1),ys(floor(length(ys)*.005)+1),ys(floor(length(ys)*.995)+1)];
end
rng=range;  % Return value
if dolog(1)
  x=log10(x);
  range(1:2)=log10(range(1:2));
end
if dolog(2)
  y=log10(y);
  range(3:4)=log10(range(3:4));
end
if dologicle(1)
  if nargin<7
    lx=Logicle(prctile(x(x<0),5));   % From Park et al paper
  else
    lx=Logicle(scaling(1));
  end
  x=lx.map(x);
  range(1:2)=[0,lx.M];
end
if dologicle(2)
  if nargin<7
    ly=Logicle(prctile(y(y<0),5));   % From Park et al paper
  else
    ly=Logicle(scaling(2));
  end
  y=ly.map(y);
  range(3:4)=[0,ly.M];
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
if dolog(1)
  mx=10.^mx;
end
if dolog(2)
  my=10.^my;
end
pcolor(mx,my,z);
if dolog(1)
  set(gca,'XScale','log');
end
if dolog(2)
  set(gca,'YScale','log');
end
if dolog(1)||dolog(2)
  logticks(dolog(1),dolog(2));
end
% Setup logicle ticks
if dologicle(1)
  tickval=sort(unique([ceil(lx.r),0,100,1000,1e4,1e5,lx.T]));
  tickpos=lx.map(tickval);
  c=axis;
  sel=tickpos>=c(1) & tickpos<=c(2);
  tickpos=tickpos(sel); tickval=tickval(sel);
  set(gca,'XTick',tickpos);
  set(gca,'XTickLabel',num2cell(tickval));
end
if dologicle(2)
  tickval=sort(unique([ceil(ly.r),0,100,1000,1e4,1e5,ly.T]));
  tickpos=ly.map(tickval);
  c=axis;
  sel=tickpos>=c(3) & tickpos<=c(4);
  tickpos=tickpos(sel); tickval=tickval(sel);
  set(gca,'YTick',tickpos);
  set(gca,'YTickLabel',num2cell(tickval));
end

allz=z(2:end-2,2:end-2);
maxcnt=prctile(allz(:),99.5);
caxis([0,maxcnt]);
shading flat
colorbar
colormap('jet');
cmap=get(gcf,'Colormap');
cmapsize=max(round(maxcnt),100);

%  fprintf('Resetting color map to contain %d entries instead of %d\n', cmapsize,size(cmap,1));
cmap=jet(cmapsize);
cmap(1,:)=1;  % Make 0 counts white
set(gcf,'Colormap',cmap);






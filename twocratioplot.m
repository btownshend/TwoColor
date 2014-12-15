% twocolor ratio plot
% Brent Townshend 10/2012
% Usage:  twocratioplot(f,gates,options)
% Gates are setup using Gates class
% Options are:
%  bins: bins for doing multi-way sorting
%  sortbins: generate this number of bins automatically
%  firstbin: ratio for lowest bin edge
%  lastbin: ratio for highest bin edge
%  normalizeratios: true to normalize the ratios
function twocratioplot(f,varargin)
defaults=struct('bins',[],...
                'sortbins',8,...
                'firstbin',[],...
                'lastbin',[],...
                'maingates',{{}},...
                'ratiorange',[1e-2,20],...
                'normalizeratios',false,...
                'title','');
args=processargs(defaults,varargin);

% Plot the GFP/mCherry ratios

edges=[];
if ~isempty(args.bins)
  edges=args.bins;
elseif ~isempty(args.maingates)
  for i=1:length(args.maingates)
    edges(i*2-1)=args.maingates{i}(1);
    edges(i*2)=args.maingates{i}(2);
  end
elseif args.sortbins>0 && sum(f(1).P(f(1).usegatenum,:))>100
  if isempty(args.firstbin) || isempty(args.lastbin)
    sortratio=sort(f(1).ratio(f(1).P(f(1).usegatenum,:)));
    midratio=exp(mean(log(sortratio(round([.05,.95]*length(sortratio))))));
    fprintf('Mid ratio = %g\n', midratio);
  end
  if isempty(args.firstbin)
    low=sortratio(sortratio<midratio);
    args.firstbin=low(round(0.5*length(low)));
    fprintf('First bin at %g\n',args.firstbin);
  end
  if isempty(args.lastbin)
    high=sortratio(sortratio>midratio);
    args.lastbin=high(round(0.90*length(high)));
    fprintf('Last bin at %g\n',args.lastbin);
  end
  edges=args.firstbin*(args.lastbin/args.firstbin).^((0:(args.sortbins-2))/(args.sortbins-2));
end


setfig([f(1).fname,' Bin edges']);clf;
legh=[];
legv={};
ti='';
pluscolors=winter(length(f));
minuscolors=spring(length(f));
for i=1:length(f)
  x=f(i).data;
  hdr=f(i).hdr;
  P=f(i).P;
  ratio=f(i).ratio(P(f(i).usegatenum,:));
  if length(ratio)<200
    fprintf('Only %d data for %s; skipping\n',length(ratio), f(i).hdr.cells);
    continue;
  end
  second=0;
  if ~isempty(strfind(hdr.cells,'+'))
    col=pluscolors(i,:);
  else
    col=minuscolors(i,:);
  end
  if args.normalizeratios
    legh(end+1)=pdfplot(ratio/median(ratio),'-',1,max(10,round(length(ratio)/200)));
  else
    legh(end+1)=pdfplot(ratio,'-',1,max(10,round(length(ratio)/200)));
  end
  set(legh(end),'Color',col);
  hold on;
  legv{end+1}=sprintf('%s (mu=%.3f,sigma=%.2f,N=%d)',hdr.cells,f(i).mu,f(i).sigma,f(i).N);
  ti=[ti,', ',hdr.cells];
end
ti=ti(3:end);  % Remove leading ', '
if ~isempty(edges)
  f(i).edges=edges;
  fprintf('Edges for %s are [%s]\nPcts: ',hdr.cells,sprintf('%.2g ',edges));
  cum=0;
  c=axis;
  for k=1:length(edges)
    plot([edges(k),edges(k)],c(3:4),':');
    fprintf('%.1f%% ',(sum(ratio<edges(k))-cum)/length(ratio)*100);
    cum=sum(ratio<edges(k));
  end
  fprintf('%.1f%%\n',(length(ratio)-cum)/length(ratio)*100);
else
  f(i).edges=nan;
end

c=axis;
c(1:2)=args.ratiorange;
axis(c);
xlabel('GFP/mCherry');
if ~isempty(args.title)
  ti=args.title;
end
title(ti,'Interpreter','none');
legend(legh,legv,'Interpreter','none','Location','NorthWest');

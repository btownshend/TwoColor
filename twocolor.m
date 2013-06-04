% twocolor
% Do 2-color analysis of a group of FCS files
% Brent Townshend 10/2012
% Usage:  twocolor({file1,file2,...},gates,options)
% Gates are setup using Gates class
% Options are:
%  usegate: name of gate to use
%  bins: bins for doing multi-way sorting
%  sortbins: generate this number of bins automatically
%  firstbin: ratio for lowest bin edge
%  lastbin: ratio for highest bin edge
%  maingates: gates for bins stored in Gates
%  normalizeratios: true to normalize the ratios
%  maxevents: maximum events per file
%  desc: cell array of descriptions for each file for plot labeling 
%    if desc contains the substring '-theo' or '+theo' those traces are distinguished
function f=twocolor(fnames,gates,varargin)
  if nargin<2
    error('Usage: twocolor(fnames,gates,[options])\n');
  end
  defaults=struct('usegate','',...
                  'bins',[],...
                  'sortbins',8,...
                  'firstbin',[],...
                  'lastbin',[],...
                  'maingates',{{}},...
                  'ratiorange',[1e-2,10],...
                  'normalizeratios',false,...
                  'maxevents',1000000,...
                  'desc',{{}});
  args=processargs(defaults,varargin);

  fprintf('Loading data from %d files\n', length(fnames));

  % Load all data into f and construct gates
  f=struct([]);
  
  
  for i=1:length(fnames)
    fprintf('Loading %s...',fnames{i});
    [unscdata,f(i).hdr,f(i).data]=fca_readfcs(fnames{i});
    fprintf('read %d events\n', size(f(i).data,1));
    if size(f(i).data,1)>args.maxevents
      fprintf('Only keeping %d/%d events\n', args.maxevents,size(f(i).data,1));
      f(i).data=f(i).data(1:args.maxevents,:);
    end
    f(i).fname=fnames{i};
    if ~isempty(args.desc)
      f(i).hdr.cells=args.desc{i};
    elseif isempty(f(i).hdr.cells)
      f(i).hdr.cells=f(i).fname;
    end
    % Channel mapping for devices used
    fsca=findchannel(f(i).hdr.par,'FSC-A',{},1);
    ssca=findchannel(f(i).hdr.par,'SSC-A',{},1);
    sscw=findchannel(f(i).hdr.par,'SSC-W',{},0);
    ssch=findchannel(f(i).hdr.par,'SSC-H',{},1);
    fsch=findchannel(f(i).hdr.par,'FSC-H',{},1);
    gfp=findchannel(f(i).hdr.par,'GFP',{'GFP-A','B1-A','FITC-A'},1);
    cherry=findchannel(f(i).hdr.par,'mCherry-A',{'Y2-A'},1);
    dapi=findchannel(f(i).hdr.par,'DAPI-A',{'V1-A'},0);  % Optional 

    % Copy channels to explicit fields in f(i)
    channels={'fsca','fsch','ssca','sscw','ssch','gfp','cherry','dapi'};
    for j=1:length(channels)
      channel=channels{j};
      f(i).(channel)=f(i).data(:,eval(channel));
    end
    if isempty(f(i).dapi)
      % Kludge - fake DAPI as 300.0 so gates work
      f(i).dapi=300*ones(size(f(i).fsca));
    end
    
    f(i).ratio=f(i).gfp./f(i).cherry;


    % Form gates
    f(i).P=gates.applyall(f(i));
    gates.printstats(f(i));
  end

  usegatenum=gates.lookup(args.usegate);
  fprintf('Using gate %s (%d)\n', args.usegate,usegatenum);
  
  % Print data summary
  fprintf(' i           Condition                       Count      GFP   resid  cherry   ratio   resid       m       b   resid GFP(600)\n');
  for i=1:length(f)
    hdr=f(i).hdr;
    P=f(i).P;
    pf=polyfit(f(i).cherry(P(usegatenum,:)),f(i).gfp(P(usegatenum,:)),1);
    presid=f(i).gfp(P(usegatenum,:))-polyval(pf,f(i).cherry(P(usegatenum,:)));
    m=pf(1);b=pf(2);
    mgfp=median(f(i).gfp(P(usegatenum,:)));
    mcherry=median(f(i).cherry(P(usegatenum,:)));
    ratio=mgfp/mcherry;
    rresid=f(i).gfp(P(usegatenum,:))-ratio*f(i).cherry(P(usegatenum,:));
    fprintf('%2d %40s %6d  %7.1f %7.1f %7.1f %7.3f %7.1f %7.3f %7.1f %7.1f %7.1f\n',i,hdr.cells,sum(P(usegatenum,:)),mgfp,std(f(i).gfp(P(usegatenum,:))),mcherry,ratio,std(rresid),m,b,std(presid),m*600+b);
  end

  if 0
    % Plot the DAPI live/dead
    figure;
    h=pdfplot(f(15).data(:,dapi),'r',1);
    hold on;
    h(2)=pdfplot(f(16).data(:,dapi),'g',1);
    c=axis;
    plot([args.dapithresh,args.dapithresh],[c(3),c(4)],':');
    xlabel('DAPI');
    legend(h,{f(15).hdr.cells,f(16).hdr.cells});
    title(sprintf('DAPI Live/Dead Theshold (%.0f)',args.dapithresh));
  end

  % Plot all the density of GFP vs. mCherry
  figure;
  lastcells='';

  slow=1e10; shigh=1;
  glow=1e10; ghigh=1;
  for i=1:length(f)
    sel=f(i).P(usegatenum,:);
    sel=f(i).cherry>0 & f(i).gfp>0;
    sch=sort(f(i).cherry(sel));slow=min(slow,10^floor(log10(sch(round(length(sch)*.02)))));shigh=max(shigh,10^ceil(log10(sch(round(length(sch)*.98)))));
    sgfp=sort(f(i).gfp(sel));glow=min(glow,10^floor(log10(sgfp(round(length(sgfp)*.02)))));ghigh=max(ghigh,10^ceil(log10(sgfp(round(length(sgfp)*.98)))));
  end
  
  for i=1:length(f)
    sel=f(i).P(usegatenum,:);
    sel=f(i).cherry>0 & f(i).gfp>0;
    gfpval=exp(mean(log(f(i).gfp(sel))));

    val=sum(f(i).gfp(sel))/sum(f(i).cherry(sel));
    subplot(ceil(length(f)/2),2,i);
    densplot(f(i).cherry(sel),f(i).gfp(sel),[],[slow shigh glow ghigh],1);
    xlabel('mCherry');
    ylabel('GFP');
    info=sprintf('%16s G=%4.0f G/C=%4.3g ',f(i).hdr.cells,gfpval,val);
    title(info,'Interpreter','none');
    disp([fnames{i},': ',info])
    lastcells=f(i).hdr.cells;
  end


  % Plot the GFP/mCherry ratios
  legh=[];
  legv={};
  ti='';

  figure;

  edges=[];
  if ~isempty(args.bins)
    edges=args.bins;
  elseif ~isempty(args.maingates)
    for i=1:length(args.maingates)
      edges(i*2-1)=args.maingates{i}(1);
      edges(i*2)=args.maingates{i}(2);
    end
  elseif args.sortbins>0
    sortratio=sort(ratio);
    endpct=.02;
    if isempty(args.firstbin)
      args.firstbin=sortratio(ceil(endpct*length(sortratio)));
    end
    if isempty(args.lastbin)
      args.lastbin=sortratio(round((1-endpct)*length(sortratio)));
    end
    fprintf('median=%.2f\n',median(ratio));

    edges=args.firstbin*(args.lastbin/args.firstbin).^((0:(args.sortbins-2))/(args.sortbins-2));
  end
  
  cols='bcmyrg';
  for i=1:length(f)
    x=f(i).data;
    hdr=f(i).hdr;
    P=f(i).P;
    ratio=f(i).ratio(P(usegatenum,:));
    if length(ratio)<100
      continue;
    end
    second=0;
    if ~isempty([strfind(hdr.cells,'+theo'),strfind(hdr.cells,'+ theo'),strfind(f(i).fname,'+theo')])
      col='g';
    elseif ~isempty([strfind(hdr.cells,'-theo'),strfind(hdr.cells,'- theo'),strfind(f(i).fname,'-theo')])
      col='r';
    else
      col=cols(mod(i-1,length(cols))+1);
    end
    if args.normalizeratios
      legh(end+1)=pdfplot(ratio/median(ratio),col,1,round(length(ratio)/200));
    else
      legh(end+1)=pdfplot(ratio,col,1,round(length(ratio)/200));
    end
    hold on;
    % Calculate std of log(normal) distribution with extreme points removed
    stdev=std(log(ratio((ratio/median(ratio))>0.25&(ratio/median(ratio))<4)));
    legv{end+1}=sprintf('%s (mu=%.2f,sigma=%.2f,N=%d)',hdr.cells,median(ratio),exp(stdev),length(ratio));
    ti=[ti,', ',hdr.cells];
  end
  ti=ti(3:end);  % Remove leading ', '
  if ~isempty(edges)
    fprintf('Edges for %s are [%s],  Pcts: ',hdr.cells,sprintf('%.2g ',edges));
    cum=0;
    c=axis
    for k=1:length(edges)
      plot([edges(k),edges(k)],c(3:4),':');
      fprintf('%.1f%% ',(sum(ratio<edges(k))-cum)/length(ratio)*100);
      cum=sum(ratio<edges(k));
    end
    fprintf('%.1f%%\n',(length(ratio)-cum)/length(ratio)*100);
  end

  c=axis;
  c(1:2)=args.ratiorange;
  axis(c);
  xlabel('GFP/mCherry');
  title(ti,'Interpreter','none');
  legend(legh,legv,'Interpreter','none');
end

% Locate which channel number corresponds to a particular field
function channel=findchannel(par,name,altnames,reqd)
  if nargin<4
    reqd=false;
  end
  channel=[];
  altnames={name,altnames{:}};
  for j=1:length(altnames)
    channel=find(strcmp({par.name},altnames{j}));
    if ~isempty(channel)
      %      fprintf('Using channel %s for %s\n', altnames{j},name);
      break;
    end
  end
  if isempty(channel) && reqd
    error('Required channel %s not detected\n',name);
  end
end

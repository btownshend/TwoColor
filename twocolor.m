% twocolor
% Do 2-color analysis of a group of FCS files
% Brent Townshend 10/2012
% Usage:  twocolor({file1,file2,...},gates,options)
% Gates are setup using Gates class
% Options are:
%  usegate: name of gate to use
%  maxevents: maximum events per file
%  desc: cell array of descriptions for each file for plot labeling 
%    if desc contains the substring '-theo' or '+theo' those traces are distinguished
function f=twocolor(fnames,gates,varargin)
  if nargin<2
    error('Usage: twocolor(fnames,gates,[options])\n');
  end
  defaults=struct('usegate','',...
                  'ratiorange',[1e-2,20],...
                  'maxevents',1000000,...
                  'compensation',0,...
                  'desc',{{}});
  args=processargs(defaults,varargin);

  if isstruct(fnames)
    fprintf('Using already loaded data\n');
    f=fnames;
    fnames={f.fname};
  else
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
      fsca=findchannel(f(i).hdr.par,'FSC-A',{'*FSC'},1);
      ssca=findchannel(f(i).hdr.par,'SSC-A',{'*SSC'},1);
      sscw=findchannel(f(i).hdr.par,'SSC-W',{},0);
      ssch=findchannel(f(i).hdr.par,'SSC-H',{},0);
      fsch=findchannel(f(i).hdr.par,'FSC-H',{},0);
      gfp=findchannel(f(i).hdr.par,'GFP',{'GFP-A','B1-A','FITC-A','525/50 [488]'},1);
      cherry=findchannel(f(i).hdr.par,'mCherry-A',{'Y2-A','610/20 [561]'},1);
      dapi=findchannel(f(i).hdr.par,'DAPI-A',{'V1-A','460/50 [405]'},0);  % Optional 
      triggerPW=findchannel(f(i).hdr.par,'Trigger Pulse Width',{},0);
      % Copy channels to explicit fields in f(i)
      channels={'fsca','fsch','ssca','sscw','ssch','gfp','cherry','dapi','triggerPW'};
      for j=1:length(channels)
        channel=channels{j};
        if ~isempty(eval(channel))
          f(i).(channel)=f(i).data(:,eval(channel));
        end
      end
      if ~isfield(f(i),'dapi') || isempty(f(i).dapi)
        % Kludge - fake DAPI as 300.0 so gates work
        f(i).dapi=300*ones(size(f(i).fsca));
      end
      
      % Apply compensation
      f(i).gfp=f(i).gfp-args.compensation*f(i).cherry;
      
      f(i).ratio=f(i).gfp./f(i).cherry;


      % Form gates
      f(i).P=gates.applyall(f(i));
      gates.printstats(f(i));
    end
  end
  
  usegatenum=gates.lookup(args.usegate);
  f.usegatenum=usegatenum;
  fprintf('Using gate %s (%d)\n', args.usegate,usegatenum);
  
  % Print data summary
  fprintf(' i           Condition                       Count      GFP   resid  cherry   ratio  stdev   resid       m       b   resid GFP(600)\n');
  for i=1:length(f)
    hdr=f(i).hdr;
    P=f(i).P;
    pf=polyfit(f(i).cherry(P(usegatenum,:)),f(i).gfp(P(usegatenum,:)),1);
    presid=f(i).gfp(P(usegatenum,:))-polyval(pf,f(i).cherry(P(usegatenum,:)));
    m=pf(1);b=pf(2);
    mgfp=median(f(i).gfp(P(usegatenum,:)));
    mcherry=median(f(i).cherry(P(usegatenum,:)));
    ratio=f(i).gfp(P(usegatenum,:))./f(i).cherry(P(usegatenum,:));
    stdev=std(log(ratio((ratio/median(ratio))>0.25&(ratio/median(ratio))<4)));
    rresid=f(i).gfp(P(usegatenum,:))-median(ratio)*f(i).cherry(P(usegatenum,:));
    fprintf('%2d %40s %6d  %7.1f %7.1f %7.1f %7.3f %7.3f %7.1f %7.3f %7.1f %7.1f %7.1f\n',i,hdr.cells,sum(P(usegatenum,:)),mgfp,std(f(i).gfp(P(usegatenum,:))),mcherry,median(ratio),exp(stdev),std(rresid),m,b,std(presid),m*600+b);
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

  % Calculate ratios
  for i=1:length(f)
    ratio=f(i).ratio(f(i).P(usegatenum,:));
    f(i).N=length(ratio);
    if length(ratio)<100
      f(i).mu=nan;
      f(i).sigma=nan;
    else
      % Calculate std of log(normal) distribution with extreme points removed
      stdev=std(log(ratio((ratio/median(ratio))>0.25&(ratio/median(ratio))<4)));
      f(i).mu=median(ratio);
      f(i).sigma=exp(stdev);
    end
  end
  
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

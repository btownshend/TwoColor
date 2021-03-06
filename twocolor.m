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
%  cherryonly: two color data structure for mCherry-only under same conditions as test sample; used for compensation
function f=twocolor(fnames,gates,varargin)
  if nargin<2
    error('Usage: twocolor(fnames,gates,[options])\n');
  end
  defaults=struct('usegate','',...
                  'ratiorange',[1e-2,20],...
                  'maxevents',1000000,...
                  'comp',[],...,
  		  'basename',{{}},...,
                  'desc',{{}});
  args=processargs(defaults,varargin);
  if isempty(args.usegate)
    error('Must specify ''usegate'' on function call');
  end
  if isstruct(fnames)
    fprintf('Using already loaded data\n');
    f=fnames;
    fnames={f.fname};
  else
    fprintf('Loading data from %d files\n', length(fnames));

    % Load all data into f and construct gates
    f=struct([]);
    
    
    for i=1:length(fnames)
      fileinfo=dir(fnames{i});
      if isempty(fileinfo) || fileinfo.bytes==0
        fprintf('File %s does not exist or is empty\n', fnames{i});
        fileinfo
        return;
      end
      fprintf('Loading %s...',fnames{i});
      try
        [unscdata,f(i).hdr,f(i).data]=fca_readfcs(fnames{i});
      catch me
        fprintf('Error while reading %s:  only analyzing first %d files\n', fnames{i},i-1);
        fnames=fnames(1:i-1);
        break;
      end
      if ~isfield(f(i).hdr,'cells')
        fprintf('*** Bad data file: %s\n', fnames{i});
        f=f(1:i-1);
        continue;
      end
      
      fprintf('read %d events (%s)\n', size(f(i).data,1),f(i).hdr.cells);
      if size(f(i).data,1)>args.maxevents
        fprintf('Only keeping %d/%d events\n', args.maxevents,size(f(i).data,1));
        f(i).data=f(i).data(1:args.maxevents,:);
      end
      f(i).fname=fnames{i};
      if ~isempty(args.basename)
        f(i).basename=args.basename{i};
      end
      if ~isempty(args.desc)
        f(i).hdr.cells=args.desc{i};
      elseif isempty(f(i).hdr.cells) 
        if ~isempty(f(i).hdr.tube)
          f(i).hdr.cells=f(i).hdr.tube;
        else
          f(i).hdr.cells=f(i).fname;
        end
      end
      % Channel mapping for devices used
      fsca=findchannel(f(i).hdr.par,'FSC-A',{'*FSC'},1);
      ssca=findchannel(f(i).hdr.par,'SSC-A',{'*SSC'},1);
      sscw=findchannel(f(i).hdr.par,'SSC-W',{},0);
      ssch=findchannel(f(i).hdr.par,'SSC-H',{},0);
      fsch=findchannel(f(i).hdr.par,'FSC-H',{},0);
      fscw=findchannel(f(i).hdr.par,'FSC-W',{},0);
      gfp=findchannel(f(i).hdr.par,'GFP',{'GFP-A','B1-A','FITC-A','525/50 [488]'},1);
      cherry=findchannel(f(i).hdr.par,'mCherry-A',{'Y2-A','610/20 [561]','PE-Texas Red-A'},1);
      dapi=findchannel(f(i).hdr.par,'DAPI-A',{'V1-A','460/50 [405]','Pacific Blue-A'},0);  % Optional 
      triggerPW=findchannel(f(i).hdr.par,'Trigger Pulse Width',{},0);
      % Copy channels to explicit fields in f(i)
      channels={'fsca','fsch','fscw','ssca','sscw','ssch','gfp','cherry','dapi','triggerPW'};
      for j=1:length(channels)
        channel=channels{j};
        if ~isempty(eval(channel))
          f(i).(channel)=f(i).data(:,eval(channel));
        end
      end
      if ~isfield(f(i),'dapi') || isempty(f(i).dapi)
        % Kludge - fake DAPI as 300.0 so gates work
        f(i).dapi=300*ones(size(f(i).fsca));
        fprintf('Faking DAPI channel\n');
      end
      
      % Apply compensation
      if ~isempty(args.comp) 
        f(i)=args.comp.apply(f(i));
        f(i).comp=args.comp;
      end

      % Computed pratio - for plasmid counting
      f(i).pratio=f(i).cherry./f(i).ssca;
      
      % GFP/mCherry ratio
      f(i).ratio=f(i).gfp./f(i).cherry;

      % Fake FSCW if it is not in file
      if isfield(f(i),'fsca') && isfield(f(i),'fsch') && ~isfield(f(i),'fscw')
        f(i).fscw=f(i).fsca./f(i).fsch;
        fprintf('Faking FSCW channel\n');
      end

      % Form gates
      f(i).gates=gates;
      f(i).P=gates.applyall(f(i));
      gates.printstats(f(i));
    end
  end
  
  usegatenum=gates.lookup(args.usegate);
  for i=1:length(f)
    f(i).usegatenum=usegatenum;
  end
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
    fc=calcmu(f(i));
    % Copy (dissimilar) struct
    fnames=fieldnames(fc);
    for fi=1:length(fnames)
      f(i).(fnames{fi})=fc.(fnames{fi});
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

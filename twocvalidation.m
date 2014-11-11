% Analyze FCS files for a two-color validation experiment
% filename is a CSV file (with first line header) containing the following columns:
% index, basename, desc, fname, ligand, ligand conc, ignore, comment
% 
% Any entries with ignore non-empty will not be used in summary data
function result=twocvalidation(varargin)
  defaults=struct('filename','layout.csv','doplot',true);
  args=processargs(defaults,varargin);
  fd=fopen(args.filename,'r');
  l=textscan(fd,'%d %s %s %s %s %s %s %s','Delimiter',',','ReturnOnError',true,'HeaderLines',1);
  fclose(fd);
  if length(l{8})==length(l{7})-1
    % Lose last field if missing final newline
    l{8}{end+1}='';
  end
  layout=struct('basename',l{2},'desc',l{3},'cond',l{5},'fname',l{4},'note',l{8},'ignore',l{7});
  for i=1:length(layout)
    layout(i).fname=['data/',layout(i).fname];
  end

  %layout=layout(1:14);   % For testing

  % Setup gates>
  g=Gates;
  g.add('All Events','x.fsca>0');
  g.addrange('mCherry','cherry',1,log10([2000,250000]),'All Events');
  g.addrange('gfp','gfp',1,log10([100,250000]),'mCherry');
  g.addrange('oneplasmid','pratio',0,[0.2,0.4],'gfp');
  g.addrange('multiplasmid','pratio',0,[0.4,2.0],'gfp');

  % Print all gates
  g.print;

  % Process all FCS files
  fprintf('Loading %d files\n', length(layout));
  f=twocolor({layout.fname},g,'usegate','mCherry','desc',{layout.desc},'basename',{layout.basename});
  twocdump(f,'results.csv');

  basenames=unique({layout.basename});
  if args.doplot
    for i=1:length(basenames)
      fprintf('Plotting %s...',basenames{i});
      sel=strcmp({layout.basename},basenames{i});
      twocdensplot(f(sel),basenames{i});
      setfig([basenames{i},'-gates']);clf;
      subplot(211);
      fsel=find(sel);
      leg={};
      for j=1:length(fsel)
        g.plot(f(fsel(j)),'mCherry','');
        hold on;
        leg{j}=f(fsel(j)).hdr.cells;
      end
      title('mCherry');
      handles=get(gca,'Children');
      legend(handles(3:3:end),leg,'Location','Best');
      subplot(212);
      for j=1:length(fsel)
        g.plot(f(fsel(j)),'gfp','');
        hold on;
      end
      title('GFP (after gating for mCherry)');
      suptitle(sprintf('Gating for %s',basenames{i}));
      twocratioplot(f(sel),'sortbins',0);
      fprintf('\n');
    end
  end
  
  % Add gating to allow us to check mCherry for only high GFP cells
  g.addrange('high GFP','gfp',1,log10([2000,250000]),'All Events');
  g.addrange('low mCherry(high GFP)','cherry',1,[0,log10(2000)],'high GFP');
  for i=1:length(f)
    sel1=g.apply(f(i),'low mCherry(high GFP)');
    sel2=g.apply(f(i),'high GFP');
    dropoutfrac=sum(sel1)/sum(sel2);
    fprintf('%3d %20.20s: mCherry Dropouts=%4.1f%% of %d events with significant GFP\n', i, f(i).hdr.cells, dropoutfrac*100,sum(sel2));
  end

  % Composite analysis
  summary=[];
  for i=1:length(basenames)
    s=struct('basename',basenames{i},'minus',[],'minusci',{{}},'plus',[],'plusci',{{}},'target','','sample',{{}},'indices',[]);
    fsel=find(strcmp({layout.basename},basenames{i}));
    s.indices=fsel;
    for j=1:length(fsel)
      if ~isempty(layout(fsel(j)).ignore)
        fprintf('Skipping %s since it is marked as bad (note: %s)\n', layout(fsel(j)).desc, layout(fsel(j)).note);
        continue;
      end
      ppos=strfind(layout(fsel(j)).desc,'+');
      if isempty(ppos)
        sample=layout(fsel(j)).desc;
      else
        sample=layout(fsel(j)).desc(1:ppos-1);
      end
      loc=find(strcmp(s.sample,sample),1,'last');
      if isempty(loc)
        s.sample{end+1}=sample;
        loc=length(s.sample);
        s.minus(loc)=nan;
        s.plus(loc)=nan;
      end
      target=layout(fsel(j)).cond;
      if ~isempty(target)
        if ~isempty(s.target) && ~strcmp(s.target,target)
          fprintf('%s has multiple targets: %s and %s\n', sample, s.target, target);
          continue;
        end
        s.target=target;
      end
      if isempty(s.target)
        if ~isnan(s.minus(loc))
          fprintf('Warning: Multiple entries for %s -target\n', sample);
          s.sample{end+1}=sample;
          loc=length(s.sample);
          s.plus(loc)=nan;
        end
        s.minus(loc)=f(fsel(j)).mu;
        s.minusci{loc}=f(fsel(j)).muci80;
      else
        if ~isnan(s.plus(loc))
          fprintf('Warning: Multiple entries for %s +target\n', sample);
          s.sample{end+1}=sample;
          loc=length(s.sample);
          s.minus(loc)=nan;
        end
        s.plus(loc)=f(fsel(j)).mu;
        s.plusci{loc}=f(fsel(j)).muci80;
      end
    end
    summary=[summary,s];
  end
  % Strip out the raw data, other than ratio, first
  fsummary=rmfield(f,{'data','fsca','fsch','ssca','ssch','gfp','cherry','dapi','pratio','fscw','P'});
  for i=1:length(fsummary)
    fsummary(i).ratio=single(fsummary(i).ratio);
  end
  result=struct('layout',layout,'f',fsummary,'gates',g,'summary',summary);
  save('fsummary.mat','-struct','result');
  result.f=f;   % Keep full data as return value from function

  % Make a validation plot
  
  twocvalidationplot(summary);
end

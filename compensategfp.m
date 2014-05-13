% Perform compensation of GFP signal
% Uses 'fcherryonly' arg for mCherry-only to figure out compensation
% Compensation is linear function of SSCA and Cherry, which is subtracted from GFP values
% Stores compensation fit in f.gfpcompfit, new gfp in f.gfpcomp, new ratio in f.ratiocomp

function f=compensategfp(f,fcherryonly,varargin)
defaults=struct('maxgfp',2000,'doplot',true);
args=processargs(defaults,varargin);

% Check Auto fluor
if ~strcmp(fcherryonly.desc,'mCherry-only')
  fprintf('Warning: 2nd arg is "%s," expected "mCherry-only"\n', fcherryonly.desc);
end

sel=fcherryonly.gfp<args.maxgfp;
A=[fcherryonly.ssca(sel),fcherryonly.cherry(sel)];
gfp=fcherryonly.gfp(sel);
A(:,end+1)=1;
fit=A\gfp;
rsqd=1-sum((gfp-A*fit).^2)/sum(gfp.^2);
fprintf('GFP(comp)=GFP-(%.4f*SSCA + %.4f*Cherry + %.0f) RSqd=%.2f\n', fit, rsqd);

% Compute compensated gfp
for i=1:length(f)
  A=[f(i).ssca,f(i).cherry];
  A(:,end+1)=1;
  f(i).gfpcompfit=fit;
  f(i).gfpcompvars={'ssca','cherry'};
  f(i).gfpcomp=f(i).gfp-A*fit;
  f(i).ratiocomp=f(i).gfpcomp./f(i).cherry;
end

if args.doplot
  ti=sprintf('Autofluor. %d:%s %s', j, fcherryonly.desc, fcherryonly.sorter);
  setfig(ti);clf;
  chrange=prctile(fcherryonly.cherry,[0,98]);
  grange=prctile([fcherryonly.gfp;fcherryonly.gfpcomp],[1,95]);
  sscrange=prctile(fcherryonly.ssca,[2,98]);

  for rpt=1:2
    if rpt==1
      gfp=fcherryonly.gfp;
    else
      gfp=fcherryonly.gfpcomp;
    end

    subplot(2,2,rpt);
    densplot(fcherryonly.cherry,gfp,[],[chrange,grange],false);
    xlabel('Cherry');
    ylabel('GFP');
    if rpt==1
      title('No autofl. subtraction');
    else
      title('With autofl. subtraction');
    end
    
    subplot(2,2,rpt+2);
    densplot(fcherryonly.ssca,gfp,[],[sscrange,grange],false);
    xlabel('SSC-A');
    ylabel('GFP');
  end
  suptitle(ti);
end

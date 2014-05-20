% Perform compensation of GFP signal
% Uses 'fcherryonly' arg for mCherry-only to figure out compensation
% Compensation is linear function of SSCA and Cherry, which is subtracted from GFP values
% Stores compensation fit in f.gfpcompfit, new gfp in f.gfpcomp, new ratio in f.ratiocomp

function f=compensategfp(f,fcherryonly,varargin)
defaults=struct('maxgfp',2000,'doplot',true);
args=processargs(defaults,varargin);

% Check Auto fluor
if ~strcmp(fcherryonly.desc,'mCherry-only') && ~strcmp(fcherryonly.desc,'mCherry only')
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
f=[f,fcherryonly];
for i=1:length(f)
  A=[f(i).ssca,f(i).cherry];
  A(:,end+1)=1;
  f(i).uncomp=struct('gfp',f(i).gfp,'ratio',f(i).ratio,'mu',f(i).mu,'sigma',f(i).sigma,'muci80',f(i).muci80);
  f(i).gfpcompfit=fit;
  f(i).gfpcompvars={'ssca','cherry'};
  f(i).gfp=f(i).uncomp.gfp-A*fit;
  f(i).ratio=f(i).gfp./f(i).cherry;
  f(i)=calcmu(f(i));
end
fcherryonly=f(end);
f=f(1:end-1);

if args.doplot
  if ~isfield(fcherryonly,'sorter')
    fcherryonly.sorter='Unknown Sorter';
  end
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

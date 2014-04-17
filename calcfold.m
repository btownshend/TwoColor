% Calculate fold change between two data
function [medianfold,ci80,map]=calcfold(fminus,fplus,varargin)
  defaults=struct('doplot',false);
  args=processargs(defaults,varargin);

  % Calculate as a function of Cherry

  edges=prctile([fminus.cherry;fplus.cherry],0:1:100);
  cminus=[];gminus=[];cplus=[];gplus=[];
  for k=1:length(edges)-1
    selminus=fminus.cherry>=edges(k) & fminus.cherry<edges(k+1) & fminus.P(fminus.usegatenum,:)';
    cminus(k)=median(fminus.cherry(selminus));
    gminus(k)=median(fminus.gfp(selminus));
    ratiominus(k)=median(fminus.gfp(selminus)./fminus.cherry(selminus));
    selplus=fplus.cherry>=edges(k) & fplus.cherry<edges(k+1) & fplus.P(fplus.usegatenum,:)';
    cplus(k)=median(fplus.cherry(selplus));
    gplus(k)=median(fplus.gfp(selplus));
    ratioplus(k)=median(fplus.gfp(selplus)./fplus.cherry(selplus));
  end
  fold=interp1(cplus(isfinite(cplus)),ratioplus(isfinite(cplus)),cminus)./ratiominus;
  medianfold=median(fold(isfinite(fold)));
  ci80=prctile(fold(isfinite(fold)),[10,90]);
  map=struct('cherry',cminus,'mu1',ratiominus,'mu2',ratioplus,'fold',fold);

  
  if args.doplot
    fprintf('Fold change:\n');
    if isfield(fminus,'oldmu') && isfield(fplus,'oldmu') && ~isempty(fminus.oldmu) && ~isempty(fplus.oldmu)
      fprintf('\tUsing ratio of old mu values: %.2f\n', fplus.oldmu/fminus.oldmu);
    end
    fprintf('\tUsing ratio of mu values: %.2f\n', fplus.mu/fminus.mu);
    fprintf('\tUsing median of fold at each cherry level: %.2f [%.2f,%.2f]\n', medianfold,ci80);
    xrange=[100,300000];
    yrange=[.9*min([ratiominus,ratioplus,fminus.mu,fplus.mu]),1.1*max([ratiominus,ratioplus,fminus.mu,fplus.mu])];
    if yrange(2)/yrange(1) < 10
      yrange=[1,10]/sqrt(10)*mean(yrange);
    end
    ti='calcfold';
    setfig(ti);clf;
    plotyaxis={'GFP/mCherry','GFP/mCherry','Fold'};
    plotnames={'-target','+target','Fold Change'};
    for r=1:3
      subplot(3,1,r);
      c=axis;
      if r==1
        loglog(cminus,ratiominus,'.');
        hold on;
        plot(xrange,fminus.mu*[1,1],':');
        yrange=[0.01,10];
      elseif r==2
        loglog(cplus,ratioplus,'.');
        hold on;
        plot(xrange,fplus.mu*[1,1],':');
        yrange=[0.01,10];
      else
        loglog((cminus+cplus)/2,fold,'.');
        hold on;
        loglog(xrange,fplus.mu/fminus.mu*[1,1],':')
        yrange=[.5,20];
      end
      axis([xrange(1),xrange(2),yrange(1),yrange(2)]);
      if (r==3)
        xlabel('cherry');
      end
      ylabel(plotyaxis{r});
      if r==1
        title(sprintf('%s -theo: mu=%.2f [%.2f,%.2f]',fminus.desc,fminus.mu,prctile(ratiominus,[10,90])));
      elseif r==2
        title(sprintf('%s +theo: mu=%.2f [%.2f,%.2f]',fplus.desc,fplus.mu,prctile(ratioplus,[10,90])));
      else
        title(sprintf('%s: Fold=%.2f [%.2f,%.2f]',fminus.desc,medianfold,ci80));
      end
    end
    suptitle(sprintf('%s@%.1fmM / %s@%.1fmM', fplus.desc, fplus.theoconc, fminus.desc, fminus.theoconc));
  end
end

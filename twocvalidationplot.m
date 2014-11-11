% Composite analysis from a summary generated during a vyb analysis
function twocvalidationplot(s)
for i=1:length(s)
  s(i).minusmean=nanmean(s(i).minus);
  s(i).plusmean=nanmean(s(i).plus);
  if ~isempty(s(i).minus)
    s(i).ciminus=[min(s(i).minus),max(s(i).minus)];
  else
    s(i).ciminus=[nan,nan];
  end
  if ~isempty(s(i).plus)
    s(i).ciplus=[min(s(i).plus),max(s(i).plus)];
  else
    s(i).ciplus=[nan,nan];
  end
end

% Bar plot
tgts=unique({s.target});
for i=1:length(tgts)
  tgt=tgts{i};
  if isempty(tgt)
    continue;
  end
  fsel1=find(strcmp({s.target},tgt));
  fsel2=find(strcmp({s.basename},'sTRSV')|strcmp({s.basename},'sTRSVctl'));
  fsel=unique([fsel2,fsel1],'stable');
  % Get rid of +target data for a different target
  pmean=[s(fsel).plusmean];
  reject=~strcmp({s(fsel).target},tgt);
  fprintf('Rejecting %d samples whose target does not match %s\n', sum(reject),tgt);
  keyboard
  pmean(reject)=nan;
  % Find sTRSV+tgt to do background subtraction
  strsvSel=find(strcmp({s.basename},'sTRSV') & strcmp({s.target},tgt));
  if ~isempty(strsvSel)
    nonspecific=s(strsvSel).plusmean-s(strsvSel).minusmean;
    fprintf('Compensating for non-specific effects of %s: subtracting %.2f from all plus measurements\n', tgt, nonspecific);
  else
    nonspecific=0;
  end
  setfig(['twocvalidation-',tgt]);clf;
  for k=1:2
    foldbar({s(fsel).basename},[s(fsel).minusmean],pmean-nonspecific,'ciminus',reshape([s(fsel).ciminus],2,[])','ciplus',reshape([s(fsel).ciplus]-nonspecific,2,[])');
    legend off
    title(tgt);
    ps=get(gca,'Position');
    width=length(fsel)*0.4/ps(3);
    pubfigure(sprintf('fig_fold%s',tgt),gcf,width,3.5);
  end
end

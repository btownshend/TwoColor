% Composite analysis from a summary generated during a vyb analysis
% Input is a structure array with:
%  basename - basename 
%  ratio(c,n) - mu
%  ci(c,n,2) - confidence intervals
%  conds{c} - conditions ('' is -target)
%  samples{n} - sample names
function twocvalidationplot(s,varargin)
defaults=struct('publish',false,'figname','twocvalidation','yrange',[]);
args=processargs(defaults,varargin);
allconds={''};
for i=1:length(s)
  allconds=union(allconds,s(i).conds);
end
% Form global matrices mean(i,c), ci(i,c,2) for allconds{c},s(i).basename
mean=nan(length(s),length(allconds));
ci=nan(length(s),length(allconds),2);
for i=1:length(s)
  for c=1:length(s(i).conds)
    cc=find(strcmp(allconds,s(i).conds{c}));
    if isempty(s(i).ratio)
      mean(i,cc)=nan;
    else
      mean(i,cc)=nanmean(s(i).ratio(c,:));
    end
    if isempty(s(i).ratio)
      ci(i,cc,1:2)=nan;
    else
      ci(i,cc,:)=[min(s(i).ratio(c,:)),max(s(i).ratio(c,:))];
    end
  end
end
strsvSel=find(strcmp(upper({s.basename}),'STRSV'));

for c=2:length(allconds)
  % Find sTRSV+tgt to do background subtraction
  nonspecific=0;
  if ~isempty(strsvSel)
    nonspecific=mean(strsvSel,c)-mean(strsvSel,1);
    if isfinite(nonspecific)
      fprintf('Compensating for non-specific effects of %s: subtracting %.2f from all plus measurements\n', allconds{c}, nonspecific);
    else
      nonspecific=0;
    end
  end

  setfig([args.figname,'-',allconds{c}]);clf;
  
  % Only plot samples that are not single-target for a different target
  sel=(isfinite(mean(:,1)) & sum(isfinite(mean),2)==1)|isfinite(mean(:,c));
  % Always keep strsv*
  sel=sel|strncmp(upper({s.basename}),'STRSV',5)';

  
  % execute the figure drawing twice so that the figure will get resized by pubfigure before the final annotations are added
  for k=1:2
    foldbar({s(sel).basename},mean(sel,1),mean(sel,c)-nonspecific,'ciminus',reshape(ci(sel,1,:),sum(sel),[])-nonspecific,'ciplus',reshape(ci(sel,c,:),sum(sel),[])-nonspecific,'yrange',args.yrange);
    if isfield(s,'section')
      % Draw section dividers
      hold on;
      sections=[s(sel).section];
      cax=axis;
      for i=2:length(sections)
        if sections(i)~=sections(i-1)
          plot((i-0.5)*[1,1],cax(3:4),':k');
        end
      end
    end
    %    legend off
    if args.publish
      ps=get(gca,'Position');
      width=length(s)*0.2/ps(3);
      pubfigure(sprintf('fig_fold%s',allconds{c}),gcf,width,3.5);
    else
      title([args.figname,'-',allconds{c}]);
    end
  end
end

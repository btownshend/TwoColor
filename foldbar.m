% Foldbar - barchar of fold changes
function foldbar(labels,minus,plus,varargin)
defaults=struct('ciminus',[],'ciplus',[],'control',[]);
args=processargs(defaults,varargin);
if ~isempty(args.control)
  minus=minus/args.control;
  plus=plus/args.control;
end
%clf;
% Don't clear, instead just delete any annnotations
% This way, we can draw on a figure after pubfigure() has been called and changed the 'Position' of the axes resulting in 
% moving any annotations
delete(findall(gcf,'Type','doubleendarrowshape'))
delete(findall(gcf,'Type','textboxshape'))
all=[minus(:);plus(:)];
bar([minus(:),plus(:)]);
set(gca,'YScale','log');
set(gca,'XTick',1:length(minus));
set(gca,'XTickLabel',labels);
set(gca,'XTickLabelRotation',45);
set(gca,'box','off');
legend('-target','+target','Location','Best');
hold on;
mpos=-.15; ppos=.15;
if ~isempty(args.ciminus)
  if size(args.ciminus,1)~=length(minus) || size(args.ciminus,2)~=2
    error('ciminus must be %d x 2', length(minus));
  end
  errorbar((1:length(minus))+mpos,minus,minus(:)-args.ciminus(:,1),args.ciminus(:,2)-minus(:),'.k');
end
if ~isempty(args.ciplus)
  if size(args.ciplus,1)~=length(plus) || size(args.ciplus,2)~=2
    error('ciplus must be %d x 2', length(plus));
  end
  errorbar((1:length(plus))+ppos,plus,plus(:)-args.ciplus(:,1),args.ciplus(:,2)-plus(:),'.k');
end
ax=axis;
ax(1)=0.5; ax(2)=length(minus)+0.5;
axis(ax);

twidth=0.05; theight=0.05;
for i=1:length(minus)
  if isnan(plus(i)) || isnan(minus(i))
    continue;
  end
  p1=pos2normalized([i+mpos*1.5,plus(i)]);
  p2=pos2normalized([i+mpos*1.5,minus(i)]);
  fold=plus(i)/minus(i);
  if fold >=2.4
    annotation('doublearrow',[p1(1),p2(1)],[min(p1(2),p2(2))+.01,max(p1(2),p2(2))-.01],'LineWidth',2);
  end
  annotation('textbox',[p1(1)-twidth/2,max(p1(2),p2(2))+.01,twidth,theight],'String',sprintf('%.1f',plus(i)/minus(i)),'HorizontalAlignment','center','VerticalAlignment','bottom','FitBoxToText','on','LineStyle','none');
end
colormap([0.8 0 0; 0 0.5 0]);
if ~isempty(args.control)
  ylabel('Normalized \mu (GFP/mCherry)');
else
  ylabel('\mu (GFP/mCherry)');
end

function c=pos2normalized(p)
gcapos=get(gca,'Position');
ax=axis;
c(1)=(p(1)-ax(1))/(ax(2)-ax(1))*gcapos(3)+gcapos(1);
c(2)=(log(p(2))-log(ax(3)))/(log(ax(4))-log(ax(3)))*gcapos(4)+gcapos(2);
%fprintf('p=(%f,%f), c=(%f,%f)\n', p,c);
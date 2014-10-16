% Apply non-linear compensation and compute gfp/mcherry ratio
% With a non-linearity, there is no "global" ratio
% Instead reference everything to a reference level of cherry (15000 for the VYB)
function f=calcmu(f,varargin)
defaults=struct('doplot',false,...
                'nonlinear',false,...
                'offset',0,...
                'refcherry',15000);
args=processargs(defaults,varargin);

% Calculate tightened ratio
sel=f.P(f.usegatenum,:)';

ratio=f.gfp./f.cherry;
stdmu=median(ratio(sel));
predict=f.cherry*stdmu;
predstd=sqrt(mean((log(predict(sel))-log(f.gfp(sel))).^2));
fprintf('Standard mu = %.3f [%.3f,%.3f] predstd=%.3f\n', stdmu, prctile(ratio(sel),10),prctile(ratio(sel),90),predstd);

if isempty(args.offset)
  % Remove any offset in the GFP signal to make the raw data closer to a ratio
  sellow=sel&f.cherry<prctile(f.cherry(sel),10);
  pf=polyfit(f.cherry(sellow),f.gfp(sellow),1);
  args.offset=pf(2);
end

gfp=f.gfp-args.offset;

ratio=max(0,gfp./f.cherry);
mu=median(ratio(sel));
predict=f.cherry*mu+args.offset;
predstd=sqrt(mean((log(predict(sel))-log(f.gfp(sel))).^2));
fprintf('Offset %.0f,  mu = %.3f [%.3f,%.3f] predstd=%.3f\n', args.offset, mu, prctile(ratio(sel),10),prctile(ratio(sel),90),predstd);


% Calculate non-linear compensation
if args.nonlinear
  fitsel=sel&gfp>prctile(gfp(sel),50)&f.cherry>prctile(f.cherry(sel),50);
  pf=polyfit(log(f.cherry(fitsel)),log(gfp(fitsel)),1);
  compratio=max(0,gfp)./((f.cherry/args.refcherry).^pf(1)*args.refcherry);
  compratio=real(compratio);
  compratio(f.cherry<0)=nan;
  mu=median(compratio(sel));
  setfig('nlfit');clf;
  loglog(f.cherry(sel),gfp(sel),'.');
  hold on;
  loglog(f.cherry(fitsel),gfp(fitsel),'.g');
  ax=axis;
  plot(ax(1:2),exp(polyval(pf,log(ax(1:2)))),':r');
  predict=((f.cherry/args.refcherry).^pf(1)*args.refcherry)*mu+args.offset;
  predstd=sqrt(mean((log(predict(sel))-log(f.gfp(sel))).^2));
  fprintf('Corrected mu = %.3f [%.3f,%.3f] predstd=%.3f (exp=%.2f)\n', mu, prctile(compratio(sel),10),prctile(compratio(sel),90),predstd,pf(1));
end


%fprintf('Using reference cherry of %.0f\n', args.refcherry);
if args.nonlinear && args.doplot
  ratiosel=sel & ratio<f.mu*3 & ratio>f.mu/3;  % Plot over this range
  % Compensation plot
  setfig('calcmu-comp');clf;

  subplot(211);
  densplot(f.cherry(ratiosel),ratio(ratiosel),[],[],true);
  xlabel('Cherry');
  ylabel('Ratio');
  c=axis;
  if isfield(f,'desc')
    desc=f.desc;
  else
    desc=f.hdr.cells;
  end
  title(sprintf('%s: uncompensated, offset=%.0f',desc,args.offset));

  subplot(212);
  densplot(f.cherry(ratiosel),compratio(ratiosel),[],[],true);
  title(sprintf('%s: compensated with %.2f',desc,pf(1)));
  xlabel('Cherry');
  ylabel('Ratio');
  axis(c);

  setfig('calcmu-hist');clf;
  pdfplot(compratio(sel),[],1);
  hold on;
  c=axis;
  plot(mu*[1,1],c(3:4),':g');
  plot(2*mu*[1,1],c(3:4),':r');
  plot(0.5*mu*[1,1],c(3:4),':r');
  xlabel('Compensated Ratio');
  ylabel('Distribution');
  title(sprintf('Ratio Distribution (N=%d)',sum(sel)));

  setfig('chdropout-norm');clf;
  normplot(log10(compratio(ratiosel)));
  hold on;
  c=axis;
  plot(log10(mu)*[1,1],c(3:4),':g');
  plot(log10(2*mu)*[1,1],c(3:4),':r');
  xlabel('log10(ratio)');

end

f.mu=mu;

% Calculate std of log(normal) distribution with extreme points removed
stdev=std(log(ratio(sel&(ratio/f.mu)>0.25&(ratio/f.mu)<4)));
f.sigma=exp(stdev);
f.muci80=prctile(ratio(sel),[10,90]);

if args.nonlinear
  f.refcherry=args.refcherry;
  f.compexponent=pf(1);
  f.compratio=compratio;
end
f.offset=args.offset;




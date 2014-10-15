% Apply non-linear compensation and compute gfp/mcherry ratio
% With a non-linearity, there is no "global" ratio
% Instead reference everything to a reference level of cherry (15000 for the VYB)
function f=calcmu(f,varargin)
defaults=struct('doplot',false,...
                'nonlinear',false,...
                'refcherry',15000);
args=processargs(defaults,varargin);

% Calculate tightened ratio
sel=f.P(f.usegatenum,:)';

% Calculate non-linear compensation
if args.nonlinear
  pf=polyfit(log(f.cherry(sel)),log(f.gfp(sel)),1);
  ratio=(f.gfp./((f.cherry/args.refcherry).^pf(1)*args.refcherry));
else
  ratio=f.gfp./f.cherry;
end

mu=median(ratio(sel));

%fprintf('Using reference cherry of %.0f\n', args.refcherry);
if args.nonlinear && args.doplot
  ratiosel=sel & f.ratio<f.mu*3 & f.ratio>f.mu/3;  % Plot over this range
  % Compensation plot
  setfig('calcmu-comp');clf;
  subplot(211);
  densplot(f.cherry(ratiosel),f.ratio(ratiosel),[],[],true);
  xlabel('Cherry');
  ylabel('Ratio');
  c=axis;
  if isfield(f,'desc')
    desc=f.desc;
  else
    desc=f.hdr.cells;
  end
  title(sprintf('%s: uncompensated',desc));
  subplot(212);
  densplot(f.cherry(ratiosel),ratio(ratiosel),[],[],true);
  title(sprintf('%s: compensated with %.2f',desc,pf(1)));
  xlabel('Cherry');
  ylabel('Ratio');
  axis(c);

  setfig('calcmu-hist');clf;
  hist(log10(ratio(ratiosel)),200);
  hold on;
  c=axis;
  plot(log10(mu)*[1,1],c(3:4),':g');
  plot(log10(2*mu)*[1,1],c(3:4),':r');
  plot(log10(0.5*mu)*[1,1],c(3:4),':r');
  xlabel('log10(ratio)');
  ylabel('Distribution');
  title(sprintf('Ratio Distribution (N=%d)',sum(ratiosel)));

  setfig('chdropout-norm');clf;
  normplot(log10(ratio(sel)));
  hold on;
  c=axis;
  plot(log10(mu)*[1,1],c(3:4),':g');
  plot(log10(2*mu)*[1,1],c(3:4),':r');
  xlabel('log10(ratio)');

  fprintf('Standard mu = %.3f [%.3f,%.3f] std=%.3f\n', f.mu, prctile(f.ratio(sel),10),prctile(f.ratio(sel),90),std(log(f.ratio(sel))));
  fprintf('Corrected mu = %.3f [%.3f,%.3f] std=%.3f (exp=%.2f)\n', mu, prctile(ratio(sel),10),prctile(ratio(sel),90),std(log(ratio(sel))),pf(1));
end

f.mu=mu;
% Calculate std of log(normal) distribution with extreme points removed
stdev=std(log(ratio(sel&(ratio/f.mu)>0.25&(ratio/f.mu)<4)));
f.sigma=exp(stdev);
f.muci80=prctile(f.ratio(sel),[10,90]);

if args.nonlinear
  f.refcherry=args.refcherry;
  f.compexponent=pf(1);
  f.compratio=ratio;
end




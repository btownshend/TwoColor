function twocdensplot(f,ti)
% Plot all the density of GFP vs. mCherry
setfig([f(1).fname,' GFPvsCherry']); clf;
lastcells='';

slow=1e10; shigh=1;
glow=1e10; ghigh=1;
for i=1:length(f)
  sel=f(i).P(f(i).usegatenum,:);
  sel=f(i).cherry>0 & f(i).gfp>0;
  sch=sort(f(i).cherry(sel));slow=min(slow,10^floor(log10(sch(round(length(sch)*.02)))));shigh=max(shigh,10^ceil(log10(sch(round(length(sch)*.98)))));
  sgfp=sort(f(i).gfp(sel));glow=min(glow,10^floor(log10(sgfp(round(length(sgfp)*.02)))));ghigh=max(ghigh,10^ceil(log10(sgfp(round(length(sgfp)*.98)))));
end

for i=1:length(f)
  sel=f(i).P(f(i).usegatenum,:);
  %    sel=f(i).cherry>0 & f(i).gfp>0;
  gfpval=exp(mean(log(f(i).gfp(sel))));

  val=sum(f(i).gfp(sel))/sum(f(i).cherry(sel));
  subplot(ceil(length(f)/2),1+(length(f)>1),i);
  densplot(f(i).cherry(sel),f(i).gfp(sel),[],[slow shigh glow ghigh],1);
  xlabel('mCherry');
  ylabel('GFP');
  info=sprintf('%16s G=%4.0f G/C=%4.3g ',f(i).hdr.cells,gfpval,val);
  title(info,'Interpreter','none');
  lastcells=f(i).hdr.cells;
end
h=suptitle(ti);
set(h,'Interpreter','none');


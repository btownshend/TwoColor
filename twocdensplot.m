function twocdensplot(f,ti)
% Plot all the density of GFP vs. mCherry
setfig([f(1).fname,' ',ti,' GFPvsCherry']); clf;
lastcells='';

slow=1e10; shigh=1;
glow=1e10; ghigh=1;
for i=1:length(f)
  sel=f(i).P(f(i).usegatenum,:)';
  sel=sel & f(i).cherry>0 & f(i).gfp>0;
  sch=prctile(f(i).cherry(sel),[2,99]);slow=min(slow,sch(1)/1.02);shigh=max(shigh,sch(2)*1.02);
  sgfp=prctile(f(i).gfp(sel),[2,99]);glow=min(glow,sgfp(1)/1.02);ghigh=max(ghigh,sgfp(2)*1.02);
end

for i=1:length(f)
  sel=f(i).P(f(i).usegatenum,:);
  if sum(sel)==0
    fprintf('No data falls in gate for %s (%s) with gate %s\n', ti, f(i).hdr.cells, f(i).gates.g{f(i).usegatenum}.name);
    continue;
  end
  %    sel=f(i).cherry>0 & f(i).gfp>0;
  gfpval=exp(mean(log(f(i).gfp(sel))));
  chval=exp(mean(log(f(i).cherry(sel))));
  
  subplot(ceil(length(f)/2),1+(length(f)>1),i);
  densplot(f(i).cherry(sel),f(i).gfp(sel),[],[slow shigh glow ghigh],1);
  xlabel('mCherry');
  ylabel('GFP');
  info=sprintf('%16s N=%d G=%4.0f R=%.0f mu=%4.3g ',f(i).hdr.cells,sum(sel),gfpval,chval,f(i).mu);
  if isfield(f(i),'offset') && f(i).offset~=0
    info=[info,sprintf('offs=%.0f ',f(i).offset)];
  end
  c=axis;
  rng=[max(min(f(i).cherry(sel)),min(f(i).gfp(sel)-f(i).offset)/f(i).mu),min(max(f(i).cherry(sel)),max(f(i).gfp(sel)-f(i).offset)/f(i).mu)];
  r=exp(log(rng(1)):.01:log(rng(2)));
  hold on;
  plot(r,r*f(i).mu+f(i).offset,'y');
  title(info,'Interpreter','none');
  lastcells=f(i).hdr.cells;
end
h=suptitle(ti);
set(h,'Interpreter','none');


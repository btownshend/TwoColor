% Create new gates structure
g=Gates;

% Add a scatter gate
% Copied down from VYB settings
% Values in FCS file are 262.144 times larger than shown on VYB screen shot
g.add('scatter','inpolygon(x.fsca,x.ssca,262.144*[50,250,650,900,900,600,200,25,50],262.144*[40,45,150,300,600,700,350,75,40])');
% Add a mcherry gate
g.add('cherry','x.cherry>1000 & x.cherry<260000','scatter');

% Print out our gates
g.print

% Analyze a pair of files using gate 'cherry'
f=twocolor({'data-theo.fcs','data+theo.fcs'},g,'usegate','cherry');

% Compare bulk with cell level measures
for i=1:length(f)
  bulkcherry(i)=median(f(i).cherry);
  bulkgfp(i)=median(f(i).gfp);
  cellratio(i)=median(f(i).ratio);
end


  
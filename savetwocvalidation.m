% savetwocvalidation - save 2-color results in a .mat file
% Strips out the raw data, other than ratio, first
function savetwocvalidation(layout,f,g,summary,filename)
if nargin<5
  filename='fsummary.mat';
end

fsummary=rmfield(f,{'data','fsca','fsch','ssca','ssch','gfp','cherry','dapi','pratio','fscw','P'});
for i=1:length(fsummary)
  fsummary(i).ratio=single(fsummary(i).ratio);
end
s=struct('layout',layout,'f',fsummary,'gates',g,'summary',summary);
save(filename,'-struct','s');


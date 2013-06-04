% From XML data file
scalefactor=31457.16*2.5;   % From <scaling_factor> entries  == 262144 * 0.12 (set as scaling factor in s/w)
maingates={[-2429.92,8387.56],[8518.62,14166.29],[14203.37,24103.90],[24140.98,40889.00],[40913.05,69315.79],[69460.14,118108.38],[118223.61,199820.65],[199992.54,261823.59]};  % From <region> entries
g=Gates;
g.add('scatter','inpolygon(x.fsca,x.ssca,262.144*[50,250,650,900,900,600,200,25,50],262.144*[40,45,150,300,600,700,350,75,40])');
g.add('viable','inpolygon(x.dapi,x.ssca,262.144*[.01,2.8,40,.01,.01],262.144*[38,36,900,900,38])','scatter');   % By eye
g.add('cherry','x.cherry>200');

mg={};
for i=1:length(maingates)
  g.add(sprintf('n%d',i),sprintf('(x.gfp./x.cherry>%f) & (x.gfp./x.cherry<%f)',maingates{i}/scalefactor),'cherry');
  for j=1:2
    mg{i}(j)=maingates{i}(j)/scalefactor;
  end
end

% Print all gates
g.print;


% Process all FCS files
ff=dir('*.fcs');
f={};
targets={'2013-02-06','2013-02-075x'};
for i=1:length(targets)
  sel=arrayfun(@(x) ~isempty(strfind(x.name,targets{i})), ff);
  fprintf('Using files %s for %s\n', sprintf('%d ',find(sel)), targets{i});
  f{i}=twocolor({ff(sel).name},g,'usegate','cherry','maingates',mg,'maxevents',100000);
  densplot(f{i}(end).dapi,f{i}(end).ssca);
end

g.plot(f{end}(end),'cherry','cherry','gfp',1,1);

% Plot gates the same way they were done on the aria
g.plot(f{end}(3),'scatter','fsca','ssca',0,1);
g.plot(f{end}(end),'viable','dapi','ssca',1,1);

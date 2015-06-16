% Dump gates from a BDFACS (aria) .xml configuration file
% Save in a .m file that can be read to create gates using Gates module
function results=bdfacs_gates(infile)
if nargin<1
  error('Usage: bdfacs_gates(infile)');
end
%global x
%if ~exist('x') || isempty(x)
x=xml2struct(infile);
%end
if ~isfield(x,'bdfacs')
  fprintf('Make sure that it is the private copy of xml2struct that is on your path\n');
  if isfield(x,'Configuration')
    error('This appears to be an INFLUX file\n');
  end
  error('This does not seem to be a BDFACS file -- outer element, "bdfacs", missing\n');
end

if isfield(x.bdfacs,'worksheet_template')
  fprintf('Template file\n');
  wt=x.bdfacs.worksheet_template;
  fprintf('  Worksheet Template: %s\n', wt.Attributes.name);
  results.template=makegates(wt.gates);
  return
end

e=x.bdfacs.experiment;
fprintf('Experiment: %s (%s)\n', e.Attributes.name, e.date.Text);
aw=e.acquisition_worksheets;
fprintf(' %s\n', aw.Attributes.name);
for iwt=1:length(aw.worksheet_template)
  wt=aw.worksheet_template{iwt};
  fprintf('  Worksheet Template: %s\n', wt.Attributes.name);
  if strcmp(wt.Attributes.name,'Analysis')
    results.template=makegates(wt.gates);
  end
end

% Load gates from specific tubes
tubes=[];
for i=1:length(e.specimen)
  s=e.specimen{i};
  for j=1:length(s.tube)
    t=s.tube{j};
    fprintf('Loading gates from specimen %s, tube %s\n', s.Attributes.name, t.Attributes.name);
    g=makegates(t.gates);
    tubes=[tubes,struct('filename',t.data_filename.Text,'gates',g.gates,'xml',g.xml)];
  end
end
results.tubes=tubes;
return;

c=x.Configuration;
fprintf('Configuration: %s, Created: %s\n', c.Name.Text, c.CreationDate.Text);
ex=c.CytometerConfiguration.Experiment;
dg={};
for it=1:length(ex.Tube)+1
  if it<=length(ex.Tube)
    if iscell(ex.Tube)
      t=ex.Tube{it};
    else
      t=ex.Tube;
    end
    fname=sprintf('Tube %d: FCS=%s',it,t.FCSFile.Text);
  else
    t=ex;
    fname='Template';
  end
  dg{end+1}='';


  fprintf('%s\n',fname);
  setfig(fname); clf;
  gates=t.DataSource.EventSource.Gates;
  suptitle(strrep(fname,'\','\\'));
  for ig=1:length(gates.Gate)
    subplot(3,ceil(length(gates.Gate)/3),ig);
    g=gates.Gate{ig};
    gtype=g.GateType.Text;
    fprintf('\tGate %d, Name:%s, Type:%s', ig, g.Name.Text, gtype);
    if isfield(g,'RegionType')
      fprintf('(%s)', g.RegionType.Text);
    end
    if isfield(g,'Input')
      fprintf(', Input:%s',g.Input.Text);
      input=g.Input.Text;
    end
    if isfield(g,'ParamX')
      xname=fixname(g.ParamX.Text);
      fprintf(', X=%s (%s)',xname, g.ScaleX.Text);
      if strcmp(g.ScaleX.Text,'Log')
        xsel=sprintf('log10(x.%s)',xname);
      elseif strcmp(g.ScaleX.Text,'Linear')
        xsel=sprintf('x.%s',xname);
      else
        error('Unsupported scale type: %s\n', g.ScaleX.Text);
      end
    else
      xname='None';

      xsel='none';
    end
    if isfield(g,'ParamY')
      yname=fixname(g.ParamY.Text);
      fprintf(', Y=%s (%s)',yname, g.ScaleY.Text);
      if strcmp(g.ScaleY.Text,'Log')
        ysel=sprintf('log10(x.%s)',yname);
      elseif strcmp(g.ScaleY.Text,'Linear')
        ysel=sprintf('x.%s',yname);
      else
        error('Unsupported scale type: %s\n', g.ScaleY.Text);
      end
    else
      yname='None';
      ysel='none';
    end
    fprintf('\n');
    if isfield(g,'Vertex')
      fprintf('\t\tVertices: ');
      p=[];
      for iv=1:length(g.Vertex)
        v=g.Vertex{iv};
        fprintf('[%s,%s] ',v.X.Text,v.Y.Text);
        p(iv,1)=str2num(v.X.Text);
        p(iv,2)=str2num(v.Y.Text);
      end
      fprintf('\n');
      if strcmp(g.RegionType.Text,'Polygon')
        plot([p(:,1);p(1,1)],[p(:,2);p(1,2)]);
        xlist=sprintf('%f,',p(:,1));xlist=xlist(1:end-1);
        ylist=sprintf('%f,',p(:,2));ylist=ylist(1:end-1);
        dg{end}=[dg{end},sprintf('g.add(''%s'',''inpolygon(%s,%s,[%s],[%s])'',''%s'');\n',g.Name.Text,xsel,ysel, xlist,ylist,g.Input.Text)];
      elseif strcmp(g.RegionType.Text,'Quadrant')
        plot(p([1,2,1,3,1,4,1,5],1),p([1,2,1,3,1,4,1,5],2));
        line=p(4,:)-p(1,:); slope=line(2)/line(1);
        fprintf(' Slope=%.2f, Intercept=%.2f',slope,p(1,2)-slope*p(1,1));
      elseif strcmp(g.RegionType.Text,'Interval')
        plot(p(1,1)*[1,1],[0,1]);
        hold on;
        plot(p(2,1)*[1,1],[0,1]);
        dg{end}=[dg{end},sprintf('g.add(''%s'',''%s>=%g & %s<=%g'',''%s'');\n',g.Name.Text,xsel,p(1,1),xsel,p(2,1),g.Input.Text)];
      else
        fprintf(' Unsupported region type');
      end
    elseif strcmp(g.GateType.Text,'AllEvents')
      dg{end}=[dg{end},sprintf('g.add(''%s'',''x.fsca>0'');\n',g.Name.Text)];
    else
      fprintf(' Unsupported gate type: %s\n',g.GateType.Text);
    end
    xlabel(xsel);
    ylabel(ysel);
    title(sprintf('Gate %d - %s',ig,g.Name.Text));
  end
end
for gloc=1:length(dg)
  outfname=sprintf('%s_%d.m',strrep(x.Configuration.Name.Text,' ','_'),gloc);
  if exist(outfname,'file')
    error('File %s already exists\n',outfname);
  end
  fd=fopen(outfname,'w');
  fprintf(fd,'%s',dg{gloc});
  fclose(fd);
  fprintf('Saved gates to %s\n',outfname);
end
end

% Adjust names of parameters
function s=fixname(s)
if strcmp(s,'FSC-H')
  s='fsch';
elseif strcmp(s,'FSC-A')
  s='fsca';
elseif strcmp(s,'FSC-W')
  s='fscw';
elseif strcmp(s,'SSC-A')
  s='ssca';
elseif strcmp(s,'SSC-H')
  s='ssch';
elseif strcmp(s,'SSC-W')
  s='sscw';
elseif strcmp(s,'FITC-A')
  s='gfp';
elseif strcmp(s,'mCherry-A')
  s='cherry';
elseif strcmp(s,'Ratio: FITC-A/mCherry-A')
  s='ratio';
elseif strcmp(s,'Pacific Blue-A')
  s='dapi';
elseif strcmp(s,'PE-Texas Red-A')
  s='cherry';
end
end


function result=makegates(xml)
GG=Gates;
g=xml.gate;
if length(g)==1
  g={g};
end
for ig=1:length(g)
  gate=g{ig};
  fprintf('   Gate: %s, Type: %s', gate.Attributes.fullname, gate.Attributes.type);
  if isfield(gate,'parent')
    fprintf(', Parent: %s',gate.parent.Text);
  end
  fprintf('\n');
  if strcmp(gate.Attributes.type,'EventSource_Classifier')
    GG.add(gate.Attributes.fullname,'1');
  elseif strcmp(gate.Attributes.type,'Region_Classifier')
    r=gate.region;
    ra=r.Attributes;
    fprintf('    Region: %s, type:%s',ra.name,ra.type);
    if strcmp(ra.type,'POLYGON_REGION') || strcmp(ra.type,'RECTANGLE_REGION')
      islog=[strcmp(gate.is_x_parameter_log.Text,'true'),strcmp(gate.is_y_parameter_log.Text,'true')];
      isscaled=[strcmp(gate.is_x_parameter_scaled.Text,'true'),strcmp(gate.is_y_parameter_scaled.Text,'true')];
      scale=[str2num(gate.x_parameter_scale_value.Text),str2num(gate.y_parameter_scale_value.Text)];
      fprintf(', x:%s(log:%d), y:%s(log:%d)\n      [', ra.xparm, islog(1), ra.yparm,islog(2));
      xv=cellfun(@(z) str2num(z.Attributes.x), r.points.point);
      yv=cellfun(@(z) str2num(z.Attributes.y), r.points.point);
      for i=1:length(xv)
        fprintf('(%f,%f) ',xv(i),yv(i));
      end
      fprintf(']\n');
      GG.addpolygon(gate.Attributes.fullname,{fixname(ra.xparm),fixname(ra.yparm)},islog,[xv;yv]',gate.parent.Text,'isscaled',isscaled,'scale',scale);
    elseif strcmp(ra.type,'INTERVAL_REGION')
      islog=strcmp(gate.is_x_parameter_log.Text,'true');
      isscaled=strcmp(gate.is_x_parameter_scaled.Text,'true');
      scale=str2num(gate.x_parameter_scale_value.Text);
      xv=cellfun(@(z) str2num(z.Attributes.x), r.points.point);
      fprintf(', x:%s (log:%d)\n      [%s]\n', ra.xparm,islog,sprintf('%f ',xv));
      GG.addrange(gate.Attributes.fullname,fixname(ra.xparm),islog,xv',gate.parent.Text,'isscaled',isscaled,'scale',scale);
    else
      fprintf('(Unknown region type)\n');
    end
  elseif strcmp(gate.Attributes.type,'NOT_Classifier')
    fprintf('     NOT: %s\n', gate.input.Text);
    GG.addnot(gate.Attributes.fullname, gate.input.Text, gate.parent.Text);
  elseif strcmp(gate.Attributes.type,'RestOf_Classifier')
    fprintf('     Rest of: %s (ignored)\n', gate.input.Text);
  else
    fprintf('Unhandled gate type\n');
    keyboard
  end
end
result=struct('gates',GG','xml',xml);
end
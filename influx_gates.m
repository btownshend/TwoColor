% Dump gates from a Carmen .xml configuration file
% Save in a .m file that can be read to create gates using Gates module
function influx_gates(infile)
if nargin<1
  error('Usage: influx_gates(infile)');
end
if ~exist('x')
  x=xml2struct(infile);
end
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

% Adjust names of parameters
function s=fixname(s)
if strcmp(s,'Trigger Pulse Width')
  s='triggerPW';
elseif strcmp(s,'525/50 [488]')
  s='gfp';
elseif strcmp(s,'610/20 [561]')
  s='cherry';
elseif strcmp(s,'460/50 [405]')
  s='dapi';
elseif strcmp(s,'FSC')
  s='fsca';
elseif strcmp(s,'SSC')
  s='ssca';
end

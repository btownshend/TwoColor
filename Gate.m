% Class for single gate definitions
classdef Gate < handle;
  properties
    name;
    parent;
    expr;    % Gate as a matlab expression
    gatetype;  % 1=expr, 2=polygon, 3=range, 4=not
    vars;   % {2} Variable names (e.g. 'fsca', or 'dapi')
    islog;   % (2) Whether to take log10 of v{1},v{2} before applying polygon
    polygon;  % (:,2) Polygon coords
    range;    % (2) Range of values [low,high]
  end
  methods
    function obj=Gate(name,parent,vars,islog,polygon,isnot)
    % Construct a gate using either Gate(name,parent,expr), Gate(name,parent,vars,islog,polygon), or Gate(name,parent,vars,islog,range)
      obj.name=name;
      obj.parent=parent;
      if nargin==6
        obj.islog=0;
        obj.gatetype=4;
        obj.vars=vars;
      elseif nargin==3
        obj.expr=vars;
        obj.gatetype=1;
      else
        obj.vars=vars;
        obj.islog=islog;
        if length(vars)==1
          obj.range=polygon;
          obj.gatetype=3;
        else
          obj.polygon=polygon;
          obj.gatetype=2;
        end
      end
    end
    
    function print(obj)
      fprintf('%-20s %3d %s\n',obj.name,obj.parent,obj.desc());
    end
    
    function s=desc(obj)   
      if obj.gatetype==1
        s=sprintf('%s',obj.expr);
      elseif obj.gatetype==2 || obj.gatetype==3 || obj.gatetype==4
        for j=1:length(obj.vars)
          if obj.islog(j)
            v{j}=['log10(',obj.vars{j},')'];
          else
            v{j}=obj.vars{j};
          end
        end
        if obj.gatetype==2
          s=sprintf('POLYGON(%s,%s) [%s]',v{1},v{2},sprintf('(%f,%f) ',obj.polygon'));
        elseif obj.gatetype==3
          s=sprintf('RANGE(%s) [%f,%f]',v{1},obj.range);
        else
          s=sprintf('NOT(%s)',v{1});
        end
      end
    end
    

    function drawgate(obj) 
    % Draw gate onto current plot
      if obj.gatetype==1
        fprintf('Unable to draw "expr" gates\n');
      elseif obj.gatetype==2
        hold on;
        for i=1:2
          if obj.islog(i)
            v(:,i)=10.^obj.polygon(:,i);
          else
            v(:,i)=obj.polygon(:,i);
          end
        end
        plot(v([1:end,1],1),v([1:end,1],2),'r');
      elseif obj.gatetype==3
        c=axis;
        hold on;
        r=obj.range;
        if obj.islog
          r=10.^r;
        end
        plot([r(1),r(1)],c(3:4),'r:');
        plot([r(2),r(2)],c(3:4),'r:');
      elseif obj.gattetype==4
        fprintf('Ignoring draw of NOT gate\n');
      else
        error('Bad gatetype: %d\n', obj.gatetype);
      end
    end

    function range=getrange(obj)
    % get the range:  [minx,maxx,miny,maxy] of the gate
    % Handle log conversions
      range=[nan,nan,nan,nan];
      if obj.gatetype==2
        poly=obj.polygon;
        if obj.islog(1)
          poly(:,1)=10.^(poly(:,1));
        end
        if obj.islog(2)
          poly(:,2)=10.^(poly(:,2));
        end
        range(1)=min(poly(:,1));
        range(2)=max(poly(:,1));
        range(3)=min(poly(:,2));
        range(4)=max(poly(:,2));
      else
        error('Gate.getrange does not support gate type %d\n', obj.gatetype);
      end
    end
  end
end


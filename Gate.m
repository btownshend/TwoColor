% Class for single gate definitions
classdef Gate < handle;
  properties
    name;
    parent;
    expr;    % Gate as a matlab expression
    gatetype;  % 1=expr, 2=polygon, 3=range, 4=not
    vars;   % {2} Variable names (e.g. 'fsca', or 'dapi')
    islog;   % (2) Whether to take log10 of v{1},v{2} before applying polygon
    isscaled;  % (2) Whether data is logicle scaled (takes precedence over islog)
    scale; %  (2) Scale if logicle scaled
    logicle;   % Logicle map
    polygon;  % (:,2) Polygon coords
    range;    % (2) Range of values [low,high]
  end
  methods
    function obj=Gate(name,gatetype,parent,varargin)
    % Construct a gate 
      defaults=struct('isscaled',[],'scale',[],'islog',[],'polygon',[],'range',[],'vars',{{}},'expr','');
      args=processargs(defaults,varargin);

      obj.gatetype=gatetype;
      obj.name=name;
      obj.parent=parent;
      if isempty(args.isscaled)
        obj.isscaled=false(1,length(args.vars));
        obj.scale=zeros(1,length(args.vars));
      else
        obj.isscaled=args.isscaled;
        obj.scale=args.scale;
      end
      obj.logicle={[],[]};
      if isempty(args.islog)
        obj.islog=false(1,length(args.vars));
      else
        obj.islog=args.islog;
      end
      if gatetype==4
        obj.vars=args.vars;
      elseif gatetype==1
        obj.expr=args.expr;
      elseif gatetype==2
        obj.vars=args.vars;
        obj.polygon=args.polygon;
      elseif gatetype==3
        obj.vars=args.vars;
        obj.range=args.range;
      else
        error('Bad gatetype: %d\n', gatetype);
      end
    end
    
    function print(obj)
      fprintf('%-20s %3d %s\n',obj.name,obj.parent,obj.desc());
    end
    
    function s=desc(obj)   
      if obj.gatetype==1
        s=sprintf('%s',obj.expr);
      elseif obj.gatetype==4
        s=sprintf('NOT(%s)',obj.vars{1});
      elseif obj.gatetype==2 || obj.gatetype==3
        for j=1:length(obj.vars)
          if obj.isscaled(j)
            v{j}=['scaled(',obj.vars{j},',',num2str(obj.scale(j)),')'];
          elseif obj.islog(j)
            v{j}=['log10(',obj.vars{j},')'];
          else
            v{j}=obj.vars{j};
          end
        end
        if obj.gatetype==2
          s=sprintf('POLYGON(%s,%s) [%s]',v{1},v{2},sprintf('(%f,%f) ',obj.polygon'));
        elseif obj.gatetype==3
          s=sprintf('RANGE(%s) [%f,%f]',v{1},obj.range);
        end
      end
    end
    

    function drawgate(obj,popts) 
    % Draw gate onto current plot
    % TODO: this needs work:  currently "dumb" about whether axes are the same as the gate parameters and what the plotting scale is
      if nargin<2
        popts='r';
        if obj.gatetype==3
          popts=[popts,':'];
        end
      end
      gscale=822;  % Unsure if this s correct -- just a rough estimate (TODO: Check Aria)
      if obj.gatetype==1
        fprintf('Unable to draw "expr" gates\n');
      elseif obj.gatetype==2
        hold on;
        for i=1:2
          if obj.isscaled(i)
            % Don't need to map since the plot will have been drawn in logicle coords already
            %l=Logicle(obj.scale(i));
            %v(:,i)=l.unmap(obj.polygon(:,i));
            v(:,i)=obj.polygon(:,i)/gscale;
          elseif obj.islog(i)
            v(:,i)=10.^obj.polygon(:,i);
          else
            v(:,i)=obj.polygon(:,i);
          end
        end
        plot(v([1:end,1],1),v([1:end,1],2),popts);
      elseif obj.gatetype==3
        c=axis;
        hold on;
        r=obj.range;
        if obj.isscaled
          % We need to map the points to real values and then, if the plot is logicle, map back to plot coordinates
          % TODO: map points if the plot is logicle
          l=Logicle(obj.scale);
          r=l.unmap(r/gscale);
        elseif obj.islog
          r=10.^r;
        end
        % UGLY hack, assume y-axis if gate is on gfp
        if strcmp(obj.vars{1},'gfp')
          plot(c(1:2),[r(1),r(1)],popts);
          plot(c(1:1),[r(2),r(2)],popts);
        else
          plot([r(1),r(1)],c(3:4),popts);
          plot([r(2),r(2)],c(3:4),popts);
        end
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
        for i=1:2
          if obj.isscaled(i)
            l=Logicle(obj.scale(i));
            poly(:,i)=l.unmap(poly(:,i));
          elseif obj.islog(i)
            poly(:,i)=10.^(poly(:,i));
          end
        end
        range(1)=min(poly(:,1));
        range(2)=max(poly(:,1));
        range(3)=min(poly(:,2));
        range(4)=max(poly(:,2));
      else
        error('Gate.getrange does not support gate type %d\n', obj.gatetype);
      end
    end
    
    function l=getLogicle(obj,i)
      if ~obj.isscaled(i)
        error('Gate.getLogicle called for object that is not scaled\n');
      end
      if isempty(obj.logicle{i})
        obj.logicle{i}=Logicle(obj.scale(i));
      end
      l=obj.logicle{i};
    end
  end
end

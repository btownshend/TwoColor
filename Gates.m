% Class for FACS gates (tuned for two-color assays)
% Brent Townshend 10/2012
% Example:
%  g=Gates;
%  g.add('scatter','x.ssca>50000&x.ssca<200000&x.fsca>50000&x.fsca<200000');
%  g.add('s1','inpolygon(x.fsca,x.fsch,[6747.16,217976.80,255422.06,257342.33,69155.92,7707.29,3866.75],[19034.68,254139.43,260142.11,231129.18,58052.07,7029.33,18034.24])','scatter');
%  g.add('s2','(x.ssch>10^2.7188915) & (x.ssch<10^4.231259)','s1');
%  g.add('P8','inpolygon(log10(x.dapi),log10(x.ssca),[1.3906,3.8241,3.0596,2.9850,1.1575],[4.4185,4.5036,3.1845,1.6100,1.5888])','s2');
%  g.add('P3','x.gfp>1000');
%  g.print
classdef Gates < handle
  properties
    g;   % Cell array of gates
    hasbeenwarned;
  end
  
  methods
    function obj=Gates
      obj.g={};
      obj.hasbeenwarned=false;
    end
    
    function add(obj,name,expr,parent)
    % Add a gate.  expr can x.* fields for FCS data access
      if nargin<4
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,1,parent,'expr',expr);
    end
    
    function addpolygon(obj,name,vnames,islog,polygon,parent,varargin) 
      defaults=struct('isscaled',[0,0],'scale',[0,0]);
      args=processargs(defaults,varargin);
      
      assert(length(vnames)==2);
      assert(length(islog)==2);
      assert(size(polygon,2)==2);
      assert(size(polygon,1)>=3);
      if nargin<6
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,2,parent,'vars',vnames,'islog',islog,'polygon',polygon,'isscaled',args.isscaled,'scale',args.scale);
    end
    
    function addrange(obj,name,vname,islog,range,parent,varargin) 
      defaults=struct('isscaled',0,'scale',0);
      args=processargs(defaults,varargin);
      assert(length(islog)==1);
      assert(length(range)==2);
      if nargin<6
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,3,parent,'vars',{vname},'islog',islog,'range',range,'isscaled',args.isscaled,'scale',args.scale);
    end
    
    function addnot(obj,name,vname,parent)
      if nargin<4
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,4,parent,'vars',{vname});
    end
    
    function num=lookup(obj,name)
      for num=1:length(obj.g)
        if strcmp(obj.g{num}.name,name)
          return;
        end
      end
      error('Unable to find gate "%s"\n', name);
    end
    
    function sel=apply(obj,x,name,parents,nowarn)
    % Apply a gate
      if nargin<5
        nowarn=false;
      end
      gscale=822;  % Unsure if this s correct -- just a rough estimate (TODO: Check Aria)
      if isempty(name)
        % No gate
        sel=true(size(x.fsca,1),1);
        return;
      elseif isnumeric(name)
        gnum=name;
      else
        gnum=obj.lookup(name);
      end
      
      if isfinite(obj.g{gnum}.parent)
        if nargin>=4 && ~isempty(parents)
          sel=parents(obj.g{gnum}.parent,:)';
        else
          sel=obj.apply(x,obj.g{gnum}.parent,[],nowarn);
        end
      else
        sel=true(size(x.fsca,1),1);
      end
      if obj.g{gnum}.gatetype==1   % Expr gate
        sel=sel&eval(obj.g{gnum}.expr);
      elseif obj.g{gnum}.gatetype==2  % Polygon gate
        v=[];
        poly=obj.g{gnum}.polygon;
        for i=1:length(obj.g{gnum}.vars)
          v(:,i)=x.(obj.g{gnum}.vars{i});
          if obj.g{gnum}.isscaled(i)
            if ~obj.hasbeenwarned
              fprintf('Polygon gate has variable %s scaled -- gates.apply handling of this is untested\n',obj.g{gnum}.vars{i});
              obj.hasbeenwarned=true;
            end
            l=obj.g{gnum}.getLogicle(i);
            v(:,i)=l.map(v(:,i));
            poly=poly/gscale;
          elseif obj.g{gnum}.islog(i);
            fracneg=mean(v(:,i)<0);
            if fracneg>0.01 && ~nowarn
              fprintf('Gates.apply: Warning, ignoring %.1f%% of events that have negative values for %s\n', fracneg*100, obj.g{gnum}.vars{i});
            end
            v(v(:,i)>0,i)=log10(v(v(:,i)>0,i));
            v(v(:,i)<=0,i)=nan;
          end
        end
        % if strcmp(obj.g{gnum}.vars{2},'fscw')
        %   keyboard;
        % end
        inp=inpolygon(v(:,1),v(:,2),poly(:,1),poly(:,2));
        if (size(inp,1)~=1) ~= (size(sel,1)~=1)
          inp=inp';
        end
        sel=sel&inp;
      elseif obj.g{gnum}.gatetype==3 % Range gate
        v=x.(obj.g{gnum}.vars{1});
        if obj.g{gnum}.isscaled
          if ~obj.hasbeenwarned
            fprintf('Range gate on %s is scaled -- gates.apply handling of this is untested\n',obj.g{gnum}.vars{1});
              obj.hasbeenwarned=true;
          end
          l=obj.g{gnum}.getLogicle(1);
          v=l.map(v);
          sel=sel&v>=(obj.g{gnum}.range(1)/gscale)&v<=(obj.g{gnum}.range(2)/gscale);
        elseif obj.g{gnum}.islog
          sel=sel&v>=10.^obj.g{gnum}.range(1)&v<=10.^obj.g{gnum}.range(2);
        else
          sel=sel&v>=obj.g{gnum}.range(1)&v<=obj.g{gnum}.range(2);
        end
      elseif obj.g{gnum}.gatetype==4 % NOT gate
        ngate=obj.g{gnum}.vars{1};
        % Apply the referenced gate
        if nargin>=4 && ~isempty(parents)
          nsel=parents(obj.lookup(ngate),:)';
        else
          nsel=obj.apply(x,ngate,[],nowarn);
        end
        % And invert
        sel=sel&~nsel;
      end
    end

    function p=applyall(obj,x)
      p=false(length(obj.g),size(x.fsca,1));
      for i=1:length(obj.g);
        p(i,:)=obj.apply(x,i,p(1:i-1,:));
      end
    end
    
    function print(obj)
      for i=1:length(obj.g)
        fprintf('%2d ',i);
        obj.g{i}.print();
      end
    end

    function printstats(obj,x)
      p=applyall(obj,x);
      total=size(p,2);
      fprintf('Gating counts for %s\n', x.hdr.cells);
      fprintf(' # %30.30s Parent Events %%Parent  %%Total\n','Gate');
      for i=1:length(obj.g)
        if isfinite(obj.g{i}.parent)
          psize=sum(p(obj.g{i}.parent,:));
        else
          psize=total;
        end
        nm=obj.g{i}.name;
        if length(nm)>30
          nm=['..',nm(end-27:end)];
        end
        fprintf('%2d %30.30s %3d   %7d  %6.2f  %6.2f %s\n',i,nm,obj.g{i}.parent,sum(p(i,:)),sum(p(i,:))/psize*100,sum(p(i,:))/total*100,obj.g{i}.desc());
      end
    end

    function plot(obj, x, name,figname, dosubplots)
      if nargin<4
        figname=name;
      end
      if nargin<5
        dosubplots=false;
      end
      if isempty(x)
        return;
      end
      gnum=obj.lookup(name);
      gate=obj.g{gnum};
      parent=gate.parent;
      if isfinite(parent)
        psel=obj.apply(x,parent);
      else
        psel=true(size(x.fsca,1),1);
      end
      if sum(psel)<10
        fprintf('Gates.plot: only %d points; not plotting\n',sum(psel));
        return;
      end
      sel=obj.apply(x,gnum);
      if ~isempty(figname)
        setfig(figname);
        clf;
      end
      if gate.gatetype<=2
        v1=gate.vars{1};
        v2=gate.vars{2};
        v1log=gate.islog(1);
        v2log=gate.islog(2);
        grange=gate.getrange();
        range=[];
        if dosubplots
          subplot(211);
          [~,range]=densplot(x.(v1),x.(v2),[],[],[v1log,v2log]);
          % Adjust range to include gate
          range(1)=min(range(1),grange(1));
          range(2)=max(range(2),grange(2));
          range(3)=min(range(3),grange(3));
          range(4)=max(range(4),grange(4));
          [~,range]=densplot(x.(v1),x.(v2),[],range,[v1log,v2log],gate.isscaled,gate.scale);
          xlabel(v1);
          ylabel(v2);
          title('All events');
          gate.drawgate();
          subplot(212);
        end
        [~,range]=densplot(x.(v1)(psel),x.(v2)(psel),[],range,[v1log,v2log],gate.isscaled,gate.scale);
        % Adjust range to include gate
        range(1)=min(range(1),grange(1));
        range(2)=max(range(2),grange(2));
        range(3)=min(range(3),grange(3));
        range(4)=max(range(4),grange(4));
        densplot(x.(v1)(psel),x.(v2)(psel),[],range,[v1log,v2log],gate.isscaled,gate.scale);
        xlabel(v1);
        ylabel(v2);
        title('Included in parent gate');
        gate.drawgate();
        
        %      plot(x.(v1)(~psel),x.(v2)(~psel),'k.');
        %hold on;
        %plot(x.(v1)(psel&~sel),x.(v2)(psel&~sel),'r.');
        % plot(x.(v1)(psel&sel),x.(v2)(psel&sel),'g.');
        % if v1log
        %   set(gca,'XScale','log');
        % end
        % if v2log
        %   set(gca,'YScale','log');
        % end
      elseif gate.gatetype==3 % Range gate
        v1=gate.vars{1};
        if dosubplots
          subplot(211);
          if gate.isscaled
            l=gate.getLogicle(1);
            v=l.map(x.(v1));
            pdfplot(v);
          else
            pdfplot(x.(v1),[],gate.islog);
          end
          xlabel(v1);
          title('All events');
          gate.drawgate();
          subplot(212);
        end
        if gate.isscaled
          l=gate.getLogicle(1)
          v=l.map(x.(v1));
          pdfplot(v(psel));
        else
          pdfplot(x.(v1)(psel),[],gate.islog);
        end
        xlabel(v1);
        title('Included in parent gate');
        gate.drawgate();
      else
        fprintf('Gates.plot: Unhandled gate type: %d\n', gate.gatetype);
      end
      h=suptitle(sprintf('Gate %s applied to %s', name, x.hdr.cells));
      set(h,'Interpreter','none');
    end

  function plotgates(obj)
    ngate=length(obj.g);
    nrow=ceil(sqrt(ngate));
    ncol=ceil(ngate/nrow);
    for i=1:length(obj.g)
      subplot(nrow,ncol,i);
      g=obj.g{i};
      g.drawgate();
      if isfield(g,'v1')
        xlabel(g.v1);
      end
      if isfield(g,'v2')
        ylabel(g.v2);
      end
      h=title(sprintf('%s (%s)',g.name,g.parent ));
      set(h,'Interpreter','none');
    end
  end
  end
  
end

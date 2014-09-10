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
  end
  
  methods
    function obj=Gates
      obj.g={};
    end
    
    function add(obj,name,expr,parent)
    % Add a gate.  expr can x.* fields for FCS data access
      if nargin<4
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,parent,expr);
    end
    
    function addpolygon(obj,name,vnames,islog,polygon,parent) 
      assert(length(vnames)==2);
      assert(length(islog)==2);
      assert(size(polygon,2)==2);
      assert(size(polygon,1)>=3);
      if nargin<6
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,parent,vnames,islog,polygon);
    end
    
    function addrange(obj,name,vname,islog,range,parent) 
      assert(length(islog)==1);
      assert(length(range)==2);
      if nargin<6
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,parent,{vname},islog,range);
    end
    
    function addnot(obj,name,vname,parent)
      if nargin<6
        parent=nan;
      elseif ~isnumeric(parent)
        parent=obj.lookup(parent);
      end
      obj.g{end+1}=Gate(name,parent,{vname},[],[],1);
    end
    
    function num=lookup(obj,name)
      for num=1:length(obj.g)
        if strcmp(obj.g{num}.name,name)
          return;
        end
      end
      error('Unable to find gate "%s"\n', name);
    end
    
    function sel=apply(obj,x,name)
    % Apply a gate
      if isempty(name)
        % No gate
        sel=true(length(x),1);
        return;
      elseif isnumeric(name)
        gnum=name;
      else
        gnum=obj.lookup(name);
      end
      
      if isfinite(obj.g{gnum}.parent)
        sel=obj.apply(x,obj.g{gnum}.parent);
      else
        sel=true(length(x),1);
      end
      if obj.g{gnum}.gatetype==1   % Expr gate
        sel=sel&eval(obj.g{gnum}.expr);
      elseif obj.g{gnum}.gatetype==2  % Polygon gate
        v=[];
        for i=1:length(obj.g{gnum}.vars)
          v(:,i)=x.(obj.g{gnum}.vars{i});
          if obj.g{gnum}.islog(i);
            fracneg=mean(v(:,i)<0);
            if fracneg>0.01
              fprintf('Warning: Ignoring %.2f%% of data that is negative\n', fracneg*100);
            end
            v(v(:,i)>0,i)=log10(v(v(:,i)>0,i));
            v(v(:,i)<=0,i)=nan;
          end
        end
        sel=sel&inpolygon(v(:,1),v(:,2),obj.g{gnum}.polygon(:,1),obj.g{gnum}.polygon(:,2));
      elseif obj.g{gnum}.gatetype==3 % Range gate
        v=x.(obj.g{gnum}.vars{1});
        if obj.g{gnum}.islog
          fracneg=mean(v<0);
          if fracneg>0.01
            fprintf('Warning: Ignoring %.2f%% of data that is negative\n', fracneg*100);
          end
          v(v>0)=log10(v(v>0));
          v(v<=0)=nan;
        end
        sel=sel&v>=obj.g{gnum}.range(1)&v<=obj.g{gnum}.range(2);
      end
    end

    function p=applyall(obj,x)
      p=false(length(obj.g),size(x.data,1));
      for i=1:length(obj.g);
        p(i,:)=obj.apply(x,i);
      end
    end
    
    function print(obj)
      for i=1:length(obj.g)
        fprintf('%2d ',i);
        obj.g{i}.print();
      end
    end

    function printstats(obj,x)
      fprintf('Population Parent Events %%Parent  %%Total\n');
      p=applyall(obj,x);
      total=size(p,2);
      for i=1:length(obj.g)
        if isfinite(obj.g{i}.parent)
          psize=sum(p(obj.g{i}.parent,:));
        else
          psize=total;
        end
        fprintf('%10s %3d   %7d  %6.2f  %6.2f %s\n',obj.g{i}.name,obj.g{i}.parent,sum(p(i,:)),sum(p(i,:))/psize*100,sum(p(i,:))/total*100,obj.g{i}.desc());
      end
    end

    function plot(obj, x, name)
      if isempty(x)
        return;
      end
      gnum=obj.lookup(name);
      gate=obj.g{gnum};
      parent=gate.parent;
      if isfinite(parent)
        psel=obj.apply(x,parent);
      else
        psel=true(size(x.data,1),1);
      end
      sel=obj.apply(x,gnum);
      setfig(name);
      clf;
      if gate.gatetype<=2
        v1=gate.vars{1};
        v2=gate.vars{2};
        v1log=gate.islog(1);
        v2log=gate.islog(2);
        subplot(211);
        [~,range]=densplot(x.(v1),x.(v2),[],[],v1log||v2log);
        xlabel(v1);
        ylabel(v2);
        title('All events');
        gate.drawgate();

        subplot(212);
        densplot(x.(v1)(psel),x.(v2)(psel),[],range,v1log||v2log);
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
      else % Range gate
        v1=gate.vars{1};
        subplot(211);
        pdfplot(x.(v1),[],gate.islog);
        xlabel(v1);
        title('All events');
        gate.drawgate();

        subplot(212);
        pdfplot(x.(v1)(psel),[],gate.islog);
        xlabel(v1);
        title('Included in parent gate');
        gate.drawgate();
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

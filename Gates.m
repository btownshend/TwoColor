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
    name;
    parent;
    expr;
  end
  
  methods
    function obj=Gates
      obj.name={};
      obj.parent=[];
      obj.expr={};
    end
    
    function add(obj,name,expr,parent)
    % Add a gate.  expr can x.* fields for FCS data access
      if nargin<4
        parent=nan;
      elseif ~isnumeric(parent)
          parent=obj.lookup(parent);
      end
      obj.name{end+1}=name;
      obj.parent=[obj.parent,parent];
      obj.expr{end+1}=expr;
    end
    
    function num=lookup(obj,name)
      num=find(strcmp(obj.name,name),1);
      if isempty(num)
        error('Unable to find gate "%s"\n', name);
      end
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
      
      if isfinite(obj.parent(gnum))
        sel=obj.apply(x,obj.parent(gnum));
      else
        sel=true(length(x),1);
      end
      sel=sel&eval(obj.expr{gnum});
    end

    function p=applyall(obj,x)
      p=false(length(obj.name),size(x.data,1));
      for i=1:length(obj.name)
        p(i,:)=obj.apply(x,i);
      end
    end
    
    function print(obj)
      for i=1:length(obj.name)
        fprintf('%-10s %3d %s\n',obj.name{i},obj.parent(i),obj.expr{i});
      end
    end

    function printstats(obj,x)
      fprintf('Population Parent  Events  %%Parent  %%Total\n');
      p=applyall(obj,x);
      total=size(p,2);
      for i=1:length(obj.name)
        if isfinite(obj.parent(i))
          psize=sum(p(obj.parent(i),:));
        else
          psize=total;
        end
        fprintf('%10s %3d %7d  %6.2f  %6.2f %s\n',obj.name{i},obj.parent(i),sum(p(i,:)),sum(p(i,:))/psize*100,sum(p(i,:))/total*100,obj.expr{i});
      end
    end

    function plot(obj, x, name, v1, v2, v1log, v2log)
      if nargin<6
        v1log=false;
      end
      if nargin<7
        v2log=false;
      end
      gnum=obj.lookup(name);
      parent=obj.parent(gnum);
      if isfinite(parent)
        psel=obj.apply(x,parent);
      else
        psel=true(length(x.(v1)),1);
      end
      sel=obj.apply(x,gnum);
      setfig(name);
      clf;
      subplot(311);
      [~,range]=densplot(x.(v1),x.(v2),[],[],v1log||v2log);
      xlabel(v1);
      ylabel(v2);
      title('All events');
      subplot(312);
      densplot(x.(v1)(psel),x.(v2)(psel),[],range,v1log||v2log);
      xlabel(v1);
      ylabel(v2);
      title('Included in parent gate');
      subplot(313);
      densplot(x.(v1)(psel&sel),x.(v2)(psel&sel),[],range,v1log||v2log);
      xlabel(v1);
      ylabel(v2);
      title('Included in this gate');

      %      plot(x.(v1)(~psel),x.(v2)(~psel),'k.');
      %hold on;
      %plot(x.(v1)(psel&~sel),x.(v2)(psel&~sel),'r.');
      % plot(x.(v1)(psel&sel),x.(v2)(psel&sel),'g.');
      if v1log
        set(gca,'XScale','log');
      end
      if v2log
        set(gca,'YScale','log');
      end
      xlabel(v1);
      ylabel(v2);
      suptitle(sprintf('Gate %s applied to %s', name, x.hdr.cells));
    end
  end
end

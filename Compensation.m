% Calculate compensation values using a set of FCS files
% GFP compensation is linear function of SSCA and Cherry, which can be subtracted from GFP values
% Cherry compensation is linear function of SSCA and GFP, which can be subtracted from cherry values
% Stores compensation fit in f.gfpcompfit, f.cherrycompfit

classdef Compensation < handle
  properties
    comp;   % Cell array of compensation structures
  end
  
  methods

    function obj=Compensation()
      comp=[];
    end

    function computeGFPCherry(obj,cherryOnlyFiles,gfpOnlyFiles, varargin)
      defaults=struct('maxgfp',2000,'maxcherry',2000,'doplot',true);
      args=processargs(defaults,varargin);

      if ~isempty(cherryOnlyFiles)
        obj.compute(cherryOnlyFiles,'gfp',{'ssca','cherry'},args.maxgfp);
        fprintf('GFP(obj.comp)=GFP-(%.4f*SSCA + %.4f*Cherry + %.0f) RSqd=%.2f\n', obj.comp(end).fit, obj.comp(end).rsqd);
      end

      if ~isempty(gfpOnlyFiles)
        obj.compute(gfpOnlyFiles,'cherry',{'ssca','gfp'},args.maxcherry);
        fprintf('Cherry(comp)=Cherry-(%.4f*SSCA + %.4f*GFP + %.0f) RSqd=%.2f\n', obj.comp(end).fit, obj.comp(end).rsqd);
      end
      if args.doplot
        for i=1:length(cherryOnlyFiles)
          obj.plot(cherryOnlyFiles{i},'desc',sprintf('Cherry-only %d',i),'compsel',1);
        end
        for i=1:length(gfpOnlyFiles)
          obj.plot(gfpOnlyFiles{i},'desc',sprintf('GFP-only %d',i),'compsel',2);
        end
      end
    end

    function plot(obj,file,varargin)
      defaults=struct('desc',[],'compsel',[]);
      args=processargs(defaults,varargin);

      g=Gates;
      g.add('All Events','x.fsca>0');
      f=twocolor({file},g,'usegate','All Events');
      if isempty(args.desc)
        args.desc=f.hdr.cells;
      end
      if isempty(args.compsel)
        usecomps=obj.comp;
      else
        usecomps=obj.comp(args.compsel);
      end
      for k=1:length(usecomps)
        c=usecomps(k);

        fc=obj.apply(f);
        tgtrange=prctile(f.(c.tgt),[1,95]);
        tgtrange(2)=tgtrange(2)*1.05;
        
        ti=sprintf('Compensation of %s on %s', c.tgt, f.hdr.cells);
        setfig(ti);clf;
        for rpt=1:2
          if rpt==1
            ff=f;
          else
            ff=fc;
          end

          tgt=ff.(c.tgt);

          for i=1:length(c.compvars)
            subplot(length(c.compvars),2,rpt+(i-1)*2);
            compvar=ff.(c.compvars{i});
            srcrange=prctile(compvar,[2,98]);
            srcrange(2)=srcrange(2)*1.05;
            densplot(compvar,tgt,[],[srcrange,tgtrange],false);
            xlabel(c.compvars{i});
            ylabel(c.tgt);
            if rpt==1
              title('Uncompensated');
            else
              title('Compensated');
            end
          end
        end
        suptitle(sprintf('%s comp(%s)=[%s]',args.desc,c.tgt,sprintf('%.3f ', c.fit)));
      end
    end

    function compute(obj, files, tgt, compvars, maxtgt)
    % Compensate tgt as a function of compvars
    % files should be 'compvars'-only runs (i.e. no tgt present)
    % e.g. tgt='gfp', compvars={'cherry','ssca'} files={'cherryonly1.fcs','cherryonly2.fcs'}
    % Uses a selection to only keep cells with tgt<=maxtgt
      A=[];
      y=[];
      npts=0;  y=[];

      g=Gates;
      g.add('All Events','x.fsca>0');
      f=twocolor(files,g,'usegate','All Events');
      for i=1:length(f)
        sel=f(i).(tgt)<=maxtgt;
        fprintf('Using %.2f%% of data in %s that has %s <= %f\n', mean(sel)*100, files{i}, tgt, maxtgt);
        for j=1:length(compvars)
          s=(f(i).(compvars{j})<250000);  % Not clipped
          if mean(s)<0.99
            fprintf('Eliminating %.2f%% of data in %s that has %s clipped\n', (1-mean(s))*100, files{i}, compvars{j});
          end
          sel=sel&s;
        end
        n=length(f(i).(tgt)(sel));
        for j=1:length(compvars)
          A(npts+1:npts+n,j)=f(i).(compvars{j})(sel);
        end
        y(npts+1:npts+n,1)=f(i).(tgt)(sel);
      end
      A(:,end+1)=1;
      fit=A\y;
      rsqd=1-sum((y-A*fit).^2)/sum(y.^2);
      c=struct('tgt',tgt,'compvars',{compvars},'fit',fit,'rsqd',rsqd,'maxtgt',maxtgt);
      obj.comp=[obj.comp,c];
    end

    function fc=apply(obj,f)
      fc=f;
      for k=1:length(obj.comp)
        % Apply compensation
        A=[];
        for j=1:length(obj.comp(k).compvars)
          A(:,j)=f.(obj.comp(k).compvars{j});
        end
        A(:,end+1)=1;
        fc.(obj.comp(k).tgt)=f.(obj.comp(k).tgt)-A*obj.comp(k).fit;
      end
    end

  end % methods
end % classdef
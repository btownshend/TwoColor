% Logicle mapping based on Parks et al, 2006 paper
%
% S=Texp(-(m-w))(exp(x-w)-p^2*exp(-(x-w)/p)+p^2-1), for x>=w
% Usage:  s=logicle(x,r,M)
%  x - display scale position [0..M]
%  r - minimum (negative) value to display
%  M - number of decades to display (default 4.5)
%  s - measurement value (up to 262144), can be negative
classdef Logicle  < handle
  properties
    M;  % Number of decades of range
    T;  % Maximum signal value
    W;	% Width of mapping in decades
    p;  % p such that W=2*p*log10(p)/p+1
    r;  % Most negative value representable
    x,y;  % Map of values from x=0:W to y=logicle value
  end
  methods
    function obj=Logicle(r,M,T)
      if r>0
        fprintf('Logicle: assuming scaling value of %f is actually -%f\n', r, r);
        r=-r;
      end
      obj.r=r;
      if nargin<3
        obj.T=262144;
      else
        obj.T=T;
      end
      if nargin<2
        obj.M=4.5;
      else
        obj.M=M;
      end
      if r==0
        obj.p=0.001;
        obj.W=0;
      else
        Wdesired=(obj.M-log10(obj.T/abs(r)))/2;
        obj.p=obj.findp(Wdesired);
        obj.W=2*obj.p*log10(obj.p)/(obj.p+1);
      end
      fprintf('r=%.0f, T=%.0f, M=%f, W=%f, p=%f\n', r, obj.T, obj.M, obj.W, obj.p);
      obj.x=0:.005:obj.M;
      obj.y=unmap(obj,obj.x);
    end

    function y=unmap(obj,x)
      sp=obj.T*10.^(-(obj.M-obj.W))*(10.^(x-obj.W)-obj.p^2*10.^(-(x-obj.W)/obj.p)+obj.p^2-1);
      sM=-obj.T*10.^(-(obj.M-obj.W))*(10.^(-x+obj.W)-obj.p^2*10.^(-(-x+obj.W)/obj.p)+obj.p^2-1);
      y=sp;
      y(x<obj.W)=sM(x<obj.W);
    end

    function x=map(obj,y)
      x=interp1(obj.y,obj.x,y);
    end
  end
  methods(Static)
    function p=findp(W)
    % Invert W=2*p*log10(p)/p+1 equation
      p=fminsearch(@(p) ((2*p*log10(p)/(p+1))-W).^2, 1);
    end
  end
end
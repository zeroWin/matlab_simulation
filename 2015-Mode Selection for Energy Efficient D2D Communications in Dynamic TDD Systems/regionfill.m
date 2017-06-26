function varargout=regionfill(fun,domain,npts,varargin)
%REGIONFILL Filled region satisfied by inequalities
%   REGIONFILL(FUN,DOMAIN) fills the region defined by FUN within
%   the box bounded by DOMAIN. FUN is either a handle of M-function
%   with 2 argruments and a logical return, or a set of inequalities
%   concatenated by logical operators in expression form (some logical
%   eqations involved in FLOOR, CEIL, and FIX is ok as well). DOMAIN
%   is a 1-by-4 vector([xmin,xmax,ymin,ymax]) denoting the plot range.
%
%   REGIONFILL(...,npts) sets the number of sample points in each
%   direction of the bounding box specified by DOAMIN. Default value
%   is 28.
%
%   REGIONFILL(...,'PropertyName',propertyvalue,...) specifys
%   properties for this function. The property names can be reduce
%   to a string without those characters in the brackets.
%
%   ----------------------------------------------------------------------
%        Properties               Values
%   ----------------------------------------------------------------------
%      'MaxRecursion'        |   0(default),1,2,...
%                            |
%      'Segments'            |   TRUE(default) or FALSE
%                            |
%      'RealSol'             |   TRUE(default) or FALSE
%                            |
%      'FaceColor'           |   'r', 'g', 'b', 'w', 'y', 'c', 'm',
%                            |   [R,G,B] array(automaticlly)
%                            |
%      'MeshColor'           |   'r', 'g', 'b', 'w', 'y', 'c', 'm',
%                            |   [R,G,B], 'none'(defualt)
%                            |
%      'MeshStyle'           |   '.-', '--', ':', '-'(default), etc.
%                            |
%      'B(or)d(er)Color'     |   'r', 'g', 'b', 'w', 'y', 'c', 'm',
%                            |   'k' or [R,G,B] array(automaticlly)
%                            |
%      'B(or)d(er)LineStyle' |   '--','-.',':','-'(default) or 'r.-', etc.
%   ----------------------------------------------------------------------  
%
%   'MaxRecursion' specifys the maximun number of recurrions to find
%   the border lines. Output argument LINES is a 2-by-M array stores
%   contours(border lines) seperated by NaNs for setting 'Segments' to
%   FALSE(0), and line segments without order for TRUE(1). Set 'RealSol'
%   to TRUE(1) for comparing in real field and FAlSE(0) for comparing
%   in real part or imaginary part. 'FaceColor'(automaticlly set by
%   program if no assignment) specifys the face color of filled region.
%   'BdColor'(short for 'BorderColor') and 'BdLineStyle'(short for
%   'BorderLineStyle') specify the color and line style of the border
%   lines respectively. 'BdColor' is set to light gray as default.
%
%   [VERTICES,LINES]=REGIONFILL(...) return a filled region and its
%   contours. VERTICES, an M-by-2 array, stores the coordinates of
%   vertices of that contours. LINES is a 1-by-k  cell array with
%   in every cell a list of indices corresponding to one connect
%   contour. If the contour is closed then the last index value is
%   equal to first index value.
%
%   [VERTICES,LINES,FV] = REGIONFILL(...) also returns structure FV,
%   which contains the fields vertices, faces. These fields correspond
%   to the Vertices, Faces patch properties of the filled region.
%
%   Example 1
%       % logical inequalities input
%       f='x^sin(x)+y*sin(y)>1||y^2/8+x^2<9';
%       regionfill(f,[-12 12,-12 12]) 
%
%       f='sqrt((4*x.^3.*y-4*x.*y.^3).^2+...
%         (y.^4-6*x.^2.*y.^2+x.^4-1).^2)<2&&cos(y)<1';
%       regionfill(f,[-3 3 -3 3],50) % 50 sample points
%
%   Example 2
%       % special equations is also supported
%       g='floor(x)/ceil(x*y)==1&x^2+sin(y)>2';
%       regionfill(g,[-4 4,-4 4],80)
%
%   Example 3
%       f='imag(x^sin(x)+y^sin(y))>0&x^sin(x)+y^sin(y)>1';
%       regionfill(f,[-6,6,-6,6],50,'RealSol',1)
%
%       % following result is differ from that of above 
%       f='imag(x^sin(x)+y^sin(y))>0&x^sin(x)+y^sin(y)>1';
%       regionfill(f,[-6,6,-6,6],50,'RealSol',0)
%
%   See also REGIONSURF, SURFP, SURFS.

%   Author's e-mail: kimist@qq.com
%   $Revision: 1.0 $  $Date: 2013/11/16 $

if nargin<2
    error('MATLAB:regionfill:NeedMoreArgs','At least 2 inputs are required.')
end

if ~isnumeric(domain) || ~isreal(domain) ||  size(domain,2)~=4
    error('MATLAB:regionfill:InvalidArg','Plot range must be a 1-by-4 real vector.')
end

% default sets
seg=true;
rsol=true;
cdflt=[0.701961, 0.745098, 0.823529]; % default FaceColor
bdc=[0.3 0.3 0.3]; % default BorderColor
bdl='-';           % default BorderLineStyle
maxrec=0;          % default # of MaxRecursion

if nargin<3
    npts=28;
elseif ~isnumeric(npts)||floor(npts)~=npts||npts<0
    error('MATLAB:regionfill:NptsNotInteger','npts must be a positive integer.')
end

if ~isempty(varargin)
    if mod(length(varargin),2)~=0
        error('MATLAB:IncompletePropsPairs',...
            'Need input property value.')
    end
    
    for k=1:2:length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case {'bdcolor','bordercolor'}
                    bdc=varargin{k+1};
                    if ~isletter(bdc)&&(~isvector(bdc)||size(bdc,2)~=3)
                        error('MATLAB:regionfill:InvalidColSpec',...
                            ['''BorderColor'' property value must be ',...
                            'either a special color character or a',...
                            '1-by-3 RGB array.'])
                    end
                case {'bdlstyle','boaderlinestyle'}
                    [bdl,bdc,ignore,msg]=colstyle(varargin{k+1});
                    if ~isempty(msg)
                        error(msg)
                    end
                case 'facecolor'
                    cdflt=varargin{k+1};
                    if ~isletter(cdflt)&&(~isvector(cdflt)||size(cdflt,2)~=3)
                        error('MATLAB:regionfill:InvalidColSpec',...
                            ['''FaceColor'' property value must be ',...
                            'either a special color character or a',...
                            '1-by-3 RGB array.'])
                    end
                case 'segments'
                    seg=varargin{k+1};
                    if ~isscalar(seg)
                        error('MATLAB:regionfill:InvalidValue',...
                            '''Segments'' must be a scalar logical value.')
                    end
                case 'realsol'
                    rsol=varargin{k+1};
                    if ~isscalar(rsol)
                        error('MATLAB:regionfill:InvalidValue',...
                            '''RealSol'' must be a scalar logical value.')
                    end
                case {'meshcolor','meshstyle'}
                    warning('MATLAB:regionfill:NotImplOption',...
                        [varargin{k} 'property will be implemented in the future.\n',...
                        'You can set them manually in this script.'])     
                 case 'maxrecursion'
                     maxrec=varargin{k+1};
                     if ~isscalar(maxrec)||~isreal(maxrec)||...
                             maxrec<0||floor(maxrec)~=maxrec
                         error('MATLAB:regionfill:WrongPropValue',...
                             ['''MaxRec'' property value must be a real ',...
                             'non-negative scalar integer.'])
                     else
                         warning('MATLAB:regionfill:NotCompleted',...
                             '''MaxRecursion'' now has not been implemented completely.')
                     end
                otherwise
                    error('MATLAB:regionfill:InvalidPropertyName',...
                [varargin{k} ' is not a valid property name.'])
            end
        else
            error('MATLAB:regionfill:InvalidType',...
                'Property name must be a string.')
        end
    end
end

[f,g,msg]=parsefun(fun); %#
if ~isempty(msg)
    error(msg);
end
xmin=domain(1);
xmax=domain(2);
ymin=domain(3);
ymax=domain(4);
dx=(xmax-xmin)/(npts+1);
dy=(ymax-ymin)/(npts+1);
if dx<0||dy<0
    error('MATLAB:regionfill:InvalidValue',...
        'The value of DOMAIN is invalid.')
end
jumpPrecision=0.02*min([dx,dy]);
x=linspace(xmin,xmax,npts+2);
y=linspace(ymin,ymax,npts+2);
[xi,yi]=meshgrid(x,y);
ap=true(size(xi));
mp=false(size(xi));

if ~isempty(f{1})&&rsol
    for k=1:length(f)
        for j = 1:numel(f{k})
            v=f{k}{j}(xi,yi);
            good=imag(v)==0;
            ap = and(good,ap);
        end
        mp=or(mp,ap);
    end
    if ~all(ap(:))
        warning('MATLAB:ineqfill:ComplexIgnored',...
            'Complex value can''t be compared. Ignore the region producing complex value.')
    end
    z=g(xi,yi);
else
    g=strrep(fun,'||','|'); %#
    g=strrep(g,'&&','&');
    g=str2func(['@(x,y)' vectorize(g)]);
    z=g(xi,yi);
end
%  Cell
%  8----4
%  |    |
%  |    |
%  1----2
V0=z(1:end-1,1:end-1);
V1=z(1:end-1,2:end);
V2=z(2:end,2:end);
V3=z(2:end,1:end-1);
V=V0+V1*2+V2*4+V3*8;
[row,col]=find(V>0);
idx=row+(col-1)*size(z,1);
v=V(sub2ind(size(V),row,col));

% Warning: DO NOT use offset method(e.g. p2=p1+[dx,0],...)£¬or it may
% lead to a precision error less than 1e-16, which makes FINDJUMP
% returns some wrong results.

% Cell corner coordinates
p1=[xi(idx) yi(idx)];
p2=[xi(idx+size(z,1)),yi(idx+size(z,1))];
p3=[xi(idx+size(z,1)+1),yi(idx+size(z,1)+1)];
p4=[xi(idx+1),yi(idx+1)];

i2=find(v==2);
i4=find(v==4);
i6=find(v==6);
i8=find(v==8);
i9=find(v==9);
i10=find(v==10);
i11=find(v==11);
i12=find(v==12);
i13=find(v==13);
i14=find(v==14);
%---case 1-----
temp=p1(i2,:);
p1(i2,:)=p2(i2,:);
p2(i2,:)=p3(i2,:);
p3(i2,:)=p4(i2,:);
p4(i2,:)=temp;

temp=p1(i4,:);
p1(i4,:)=p3(i4,:);
p3(i4,:)=temp;
temp=p2(i4,:);
p2(i4,:)=p4(i4,:);
p4(i4,:)=temp;

temp=p4(i8,:);
p4(i8,:)=p3(i8,:);
p3(i8,:)=p2(i8,:);
p2(i8,:)=p1(i8,:);
p1(i8,:)=temp;

%----case 3------
temp=p1(i6,:);
p1(i6,:)=p2(i6,:);
p2(i6,:)=p3(i6,:);
p3(i6,:)=p4(i6,:);
p4(i6,:)=temp;

temp=p4(i9,:);
p4(i9,:)=p3(i9,:);
p3(i9,:)=p2(i9,:);
p2(i9,:)=p1(i9,:);
p1(i9,:)=temp;

temp=p1(i12,:);
p1(i12,:)=p3(i12,:);
p3(i12,:)=temp;
temp=p2(i12,:);
p2(i12,:)=p4(i12,:);
p4(i12,:)=temp;

%-----case7------
temp=p4(i11,:);
p4(i11,:)=p3(i11,:);
p3(i11,:)=p2(i11,:);
p2(i11,:)=p1(i11,:);
p1(i11,:)=temp;

temp=p1(i13,:);
p1(i13,:)=p3(i13,:);
p3(i13,:)=temp;
temp=p2(i13,:);
p2(i13,:)=p4(i13,:);
p4(i13,:)=temp;

temp=p1(i14,:);
p1(i14,:)=p2(i14,:);
p2(i14,:)=p3(i14,:);
p3(i14,:)=p4(i14,:);
p4(i14,:)=temp;

%--case5--
temp=p1(i10,:);
p1(i10,:)=p2(i10,:);
p2(i10,:)=p3(i10,:);
p3(i10,:)=p4(i10,:);
p4(i10,:)=temp;

case1=[find(v==1);i2;i4;i8];
case3=[find(v==3);i6;i9;i12];
case5=[find(v==5);i10];
case7=[find(v==7);i11;i13;i14];

x15=[p1(v==15,1)';p2(v==15,1)';p3(v==15,1)';p4(v==15,1)'];
y15=[p1(v==15,2)';p2(v==15,2)';p3(v==15,2)';p4(v==15,2)'];
cax=newplot;
patch('xdata',x15,'ydata',y15,'facecolor',cdflt,'edgecolor','none','parent',cax)
lx=[x15([1,2],:) x15([2,3],:) x15([3,4],:) x15([4,1],:)];
ly=[y15([1,2],:) y15([2,3],:) y15([3,4],:) y15([4,1],:)];

% process case1
touch1=findjump(g,p1(case1,:),p4(case1,:),jumpPrecision);
touch2=findjump(g,p1(case1,:),p2(case1,:),jumpPrecision);
x1=[touch1(:,1)'; p1(case1,1)'; touch2(:,1)'];
y1=[touch1(:,2)'; p1(case1,2)'; touch2(:,2)'];
patch('xdata',x1,'ydata',y1,'facecolor',cdflt,'edgecolor','none','parent',cax)
lx=[lx x1([1,2],:) x1([2,3],:) x1([3,1],:)];
ly=[ly y1([1,2],:) y1([2,3],:) y1([3,1],:)];
% process case3
touch1=findjump(g,p1(case3,:),p4(case3,:),jumpPrecision);
touch2=findjump(g,p2(case3,:),p3(case3,:),jumpPrecision);
x3=[touch1(:,1)'; p1(case3,1)'; p2(case3,1)'; touch2(:,1)'];
y3=[touch1(:,2)'; p1(case3,2)'; p2(case3,2)'; touch2(:,2)'];
patch('xdata',x3,'ydata',y3,'facecolor',cdflt,'edgecolor','none','parent',cax)
lx=[lx x3([1,2],:) x3([2,3],:) x3([3,4],:) x3([4,1],:)];
ly=[ly y3([1,2],:) y3([2,3],:) y3([3,4],:) y3([4,1],:)];

% process case7
touch1=findjump(g,p1(case7,:),p4(case7,:),jumpPrecision);
touch2=findjump(g,p3(case7,:),p4(case7,:),jumpPrecision);
x7=[touch1(:,1)'; p1(case7,1)'; p2(case7,1)'; p3(case7,1)'; touch2(:,1)'];
y7=[touch1(:,2)'; p1(case7,2)'; p2(case7,2)'; p3(case7,2)'; touch2(:,2)'];
patch('xdata',x7,'ydata',y7,'facecolor',cdflt,'edgecolor','none','parent',cax)
lx=[lx x7([1,2],:) x7([2,3],:) x7([3,4],:) x7([4,5],:) x7([5,1],:)];
ly=[ly y7([1,2],:) y7([2,3],:) y7([3,4],:) y7([4,5],:) y7([5,1],:)];

% process case5 (ambitious)
touch1=findjump(g,p1(case5,:),p4(case5,:),jumpPrecision);
touch2=findjump(g,p1(case5,:),p2(case5,:),jumpPrecision);
touch3=findjump(g,p3(case5,:),p2(case5,:),jumpPrecision);
touch4=findjump(g,p3(case5,:),p4(case5,:),jumpPrecision);

mid=(p1(case5,:)+p3(case5,:))/2;
asy=g(mid(:,1),mid(:,2))>=0;
if maxrec>0  % for special purpose
    touch5=findjump(g,mid,p4(case5,:),jumpPrecision);
    touch6=findjump(g,mid,p2(case5,:),jumpPrecision);
    x5=[touch1(asy,1)'; p1(case5(asy),1)'; touch2(asy,1)';...
        touch6(asy,1)'; touch3(asy,1)'; p3(case5(asy),1)';...
        touch4(asy,1)'; touch5(asy,1)'];
    y5=[touch1(asy,2)'; p1(case5(asy),2)'; touch2(asy,2)';...
        touch6(asy,2)'; touch3(asy,2)'; p3(case5(asy),2)';...
        touch4(asy,2)'; touch5(asy,2)'];
    patch('xdata',x5,'ydata',y5,'facecolor',cdflt,'edgecolor','none','parent',cax)
    lx=[lx x5([1,2],:) x5([2,3],:) x5([3,4],:) x5([4,5],:) x5([5,6],:) ...
        x5([6,7],:) x5([7,8],:) x5([8,1],:)];
    ly=[ly y5([1,2],:) y5([2,3],:) y5([3,4],:) y5([4,5],:) y5([5,6],:) ...
        y5([6,7],:) y5([7,8],:) y5([8,1],:)];
    sx=[x15([1,2,3],:) x15([3 4 1],:), x1, x3([1,2,3],:) x3([3 4 1],:),...
        x7([1,2,3],:) x7([1 3 5],:) x7([3 4 5],:),...
        x5([1 2 8],:) x5([8 2 6],:) x5([4 6 2],:) x5([3 4 2],:),...
        x5([5 6 4],:) x5([7 6 8],:)];
    sy=[y15([1,2,3],:) y15([3 4 1],:), y1, y3([1,2,3],:) y3([3 4 1],:),...
        y7([1,2,3],:) y7([1 3 5],:) y7([3 4 5],:),...
        y5([1 2 8],:) y5([8 2 6],:) y5([4 6 2],:) y5([3 4 2],:),...
        y5([5 6 4],:) y5([7 6 8],:)];
    asy=~asy;
    touch5=findjump(g,mid,p1(case5,:),jumpPrecision);
    touch6=findjump(g,mid,p3(case5,:),jumpPrecision);
    x5=[[touch1(asy,1)'; p1(case5(asy),1)'; touch2(asy,1)'; touch5(asy,1)'],...
        [touch3(asy,1)'; p3(case5(asy),1)'; touch4(asy,1)'; touch6(asy,1)']];
    y5=[[touch1(asy,2)'; p1(case5(asy),2)'; touch2(asy,2)'; touch5(asy,2)'],...
        [touch3(asy,2)'; p3(case5(asy),2)'; touch4(asy,2)'; touch6(asy,2)']];
    patch('xdata',x5,'ydata',y5,'facecolor',cdflt,'edgecolor','g','parent',cax)
    lx=[lx x5([1,2],:) x5([2,3],:) x5([3,4],:) x5([4,1],:)];
    ly=[ly y5([1,2],:) y5([2,3],:) y5([3,4],:) y5([4,1],:)];
    sx=[sx x5([1 2 4],:) x5([2 3 4],:)];
    sy=[sy y5([1 2 3],:) y5([2 3 4],:)];
else
% subcase 1
x5=[touch1(asy,1)'; p1(case5(asy),1)'; touch2(asy,1)'; touch3(asy,1)'; p3(case5(asy),1)'; touch4(asy,1)'];
y5=[touch1(asy,2)'; p1(case5(asy),2)'; touch2(asy,2)'; touch3(asy,2)'; p3(case5(asy),2)'; touch4(asy,2)'];
patch('xdata',x5,'ydata',y5,'facecolor',cdflt,'edgecolor','none','parent',cax)
lx=[lx x5([1,2],:) x5([2,3],:) x5([3,4],:) x5([4,5],:) x5([5,6],:) x5([6,1],:)];
ly=[ly y5([1,2],:) y5([2,3],:) y5([3,4],:) y5([4,5],:) y5([5,6],:) y5([6,1],:)];
sx=[x15([1,2,3],:) x15([3 4 1],:), x1, x3([1,2,3],:) x3([3 4 1],:),...
    x7([1,2,3],:) x7([1 3 5],:) x7([3 4 5],:),...
    x5([1 2 3],:) x5([1 3 4],:) x5([1 4 6],:) x5([4 5 6],:)];
sy=[y15([1,2,3],:) y15([3 4 1],:), y1, y3([1,2,3],:) y3([3 4 1],:),...
    y7([1,2,3],:) y7([1 3 5],:) y7([3 4 5],:),...
    y5([1 2 3],:) y5([1 3 4],:) y5([1 4 6],:) y5([4 5 6],:)];
% subcase 2
asy=~asy;
x5=[[touch1(asy,1)'; p1(case5(asy),1)'; touch2(asy,1)'],[touch3(asy,1)'; p3(case5(asy),1)'; touch4(asy,1)']];
y5=[[touch1(asy,2)'; p1(case5(asy),2)'; touch2(asy,2)'],[touch3(asy,2)'; p3(case5(asy),2)'; touch4(asy,2)']];
patch('xdata',x5,'ydata',y5,'facecolor',cdflt,'edgecolor','y','parent',cax)
xlabel('x'); ylabel('y');
lx=[lx x5([1,2],:) x5([2,3],:) x5([3,1],:)];
ly=[ly y5([1,2],:) y5([2,3],:) y5([3,1],:)];
sx=[sx x5];
sy=[sy y5];
end

if ~isempty(lx)
% cut off the 'tails' or make up the 'gap'(beyond 10 digits)
% to avoid possible computational errors.(If boundary line is not
% proper, please ger rid of the comment symbols in the following
% block of codes.)

% lx=floor(1e10*lx)/1e10;
% ly=floor(1e10*ly)/1e10;

[v,ignore,n]=unique([lx(:) ly(:)],'rows');
edge=reshape(n,size(lx));
[D,ignore,n]=uniquerow(sort(edge)');
line=D(n==1,:)';
 hold on
 plot(cax,[v(line(1,:),1)';v(line(2,:),1)'],...
     [v(line(1,:),2)'; v(line(2,:),2)'],...
     'color',bdc,'linestyle',bdl)
 
% You can set grid line style or color manually here.
%   plot(cax,[v(line1(1,:),1)';v(line1(2,:),1)'],...
%      [v(line1(1,:),2)'; v(line1(2,:),2)'],'--','Color','g')
 hold off
if ~seg  %#
 line=sortlines(line');
end
[v,ignore,n]=unique([sx(:) sy(:)],'rows');
fvs=struct('Vertices',v,'Faces',reshape(n,size(sx))');
else
line=cell(1,0); %#
v=zeros(0,3);
fvs=struct('Vertices',[],'Faces',[]);
end
if nargout<4&&nargout>0
varargout{1}=v;
varargout{2}=line;
varargout{3}=fvs;
elseif nargout>3
    error('MATLAB:regionfill:RedundantOutputs','Too many outputs.')
end
    


function jpts=findjump(phase,l,r,tol)
jpts=zeros(size(l));
for k=1:size(l,1)
    if sqrt(sum((l(k,:)-r(k,:)).^2))<tol
        jpts(k,:)=(l(k,:)+r(k,:))*0.5;
    elseif phase(l(k,1),l(k,2))==phase(0.5*(l(k,1)+r(k,1)),0.5*(l(k,2)+r(k,2)))
        jpts(k,:)=findjump(phase,0.5*(l(k,:)+r(k,:)),r(k,:),tol);
    else
        jpts(k,:)=findjump(phase,l(k,:),0.5*(l(k,:)+r(k,:)),tol);
    end
end

function [f,g,msg]=parsefun(expr)
msg=[];
f=[];
g=[];
if isa(expr,'function_handle')
    funi = func2str(expr);
    if ~isempty(strfind(funi,'('))||~isempty(strfind(funi,')'))||...
            ~isempty(strfind(funi,'['))||~isempty(strfind(funi,']'))||...
            ~isempty(strfind(funi,'>'))||~isempty(strfind(funi,'<'))||...
            ~isempty(strfind(funi,'='))
        msg='If FUN is a function handle, only a M-function is valid.';
        return;
    end
    g=expr;
elseif isa(expr,'char')
    pos=findstr(expr,'=');
    if ~isempty(pos)
        for k=1:length(pos)
            if pos(k)==1||pos(k)==length(expr)||...
                    expr(pos(k)+1)=='<'||expr(pos(k)+1)=='>'||...
                    expr(pos(k)-1)~='~'&&expr(pos(k)-1)~='<'&&expr(pos(k)-1)~='>'&&expr(pos(k)-1)~='='&&expr(pos(k)+1)~='='||...%##
                    (expr(pos(k)-1)~='~'&&expr(pos(k)-2)=='='||expr(pos(k)-2)=='<'||expr(pos(k)-2)=='>')&&expr(pos(k)-1)=='='%###
                
                msg='Invalid logical function expression. use ''=='', ''<='', or ''>='' instead.';
                return;
            end
        end
    end
    
    expr=strrep(expr,'&&','&');
    expr=strrep(expr,'||','|');
    C=textscan(expr,'%s','delimiter','|','MultipleDelimsAsOne', 1);
    C=C{1};
    f=cell(1,length(C));
    g='';
    for k=1:length(C)
        R=textscan(C{k},'%s','delimiter','&','MultipleDelimsAsOne', 1);
        R=R{1};
        a=[];
        for jr=1:length(R)
            if ~isempty(findstr(R{jr},'<='))
                ff=[strrep(R{jr},'<=','-(') ')'];
            elseif ~isempty(findstr(R{jr},'<'))
                ff=[strrep(R{jr},'<','-(') ')'];
            elseif ~isempty(findstr(R{jr},'>='))
                ff=['-(' strrep(R{jr},'>=',')+')];
            elseif ~isempty(findstr(R{jr},'>'))
                ff=['-(' strrep(R{jr},'>',')+')];
            elseif ~isempty(findstr(R{jr},'=='))
                ff=[strrep(R{jr},'==','-(') ')']; %#########
            elseif ~isempty(findstr(R{jr},'~='))
                ff=[strrep(R{jr},'~=','-(') ')']; %#########
            else
                msg='Expression may need logical operator.';
                return;
            end
            f{k}{jr}=str2func(vectorize(['@(x,y)' ff]));
            if isempty(findstr(R{jr},'imag'))&&isempty(findstr(R{jr},'real')) %###############
                if jr==1
                    a=[R{jr} '&abs(imag(' ff '))==0'];
                else % '&abs(imag(' ff '))==0' can be replaced by '&1e-14*abs(real(' ff '))'
                    a=[a '&' R{jr} '&abs(imag(' ff '))==0'];
                end
            else
                if jr==1
                    a=R{jr};
                else % '&abs(imag(' ff '))==0' can be replaced by '&1e-14*abs(real(' ff '))'
                    a=[a '&' R{jr}];
                end
            end
        end
        if ~isempty(ff)
            if k==1
                g=a;
            else
                g=[a '|' g];
            end
        end
    end
    g=str2func(['@(x,y)' vectorize(g)]);
else
    msg='Invalid logical function expression.';
end

function [D,nr,c] = uniquerow(R)
% UNIQUE  Extraction of unique rows out of matrix
%	[D,NR,C] = UNIQUE(R)  Returns matrix D containing unique rows of
%   input matrix R, vector NR (size(R,1) by 1) showing index into rows
%	of D for each row of R, vector C (size(D,1) by 1) containing number
%   of occurences (count) in R of each row of D.

%  Based on the program MUNIQUE.M by Richard Aufrichtig
%  (winner of M-file contest V3.0)

y = R*rand(size(R,2),1);
[y,i] = sort(y);
y = find([1; diff(y)]);

nr = zeros(size(R,1),1);
nr(y) = [i(1); diff(i(y))];
nr(i) = cumsum(nr);
y = sort(nr);
c = find(diff([y; length(nr)+1]));
c = [c(1); diff(c)];

y = zeros(size(nr));
y(nr) = ones(size(nr));
i = find(y);
y(i) = 1:length(i);
nr = y(nr);
D(nr,:) = R;

function Objects=sortlines(Lines)
% sort the order of Lines
% This part of code is stolen from ISOCONTOUR
% written by D.Kroon University of Twente (March 2011)
if isempty(Lines)
    Objects=cell(1,0);
    return
end
Obj=zeros(100,2); 
Obj(1,1)=1; nObjects=1; reverse=false;
for i=1:size(Lines,1)-1
    F=Lines(i,2);
    Q=i+1;
    R=find(any(Lines(Q:end,:)==F,2),1,'first');
    if(~isempty(R))
        R=R+i;
        TF=Lines(Q,:);
        if(Lines(R,1)==F)
            Lines(Q,:)=Lines(R,[1 2]);
        else
            Lines(Q,:)=Lines(R,[2 1]);
        end
        if(R~=Q)
            Lines(R,:)=TF;
        end
    else
        F=Lines(Obj(nObjects,1),1);
        R=find(any(Lines(Q:end,:)==F,2),1,'first');
        if(~isempty(R))
            reverse=true;
            Lines(Obj(nObjects,1):i,:)=Lines(i:-1:Obj(nObjects,1),[2 1]);
            R=R+i;
            TF=Lines(Q,:);
            if(Lines(R,1)==F),Lines(Q,:)=Lines(R,[1 2]); else Lines(Q,:)=Lines(R,[2 1]); end
            if(R~=Q), Lines(R,:)=TF; end
        else
            if(reverse)
                Lines(Obj(nObjects,1):i,:)=Lines(i:-1:Obj(nObjects,1),[2 1]);
                reverse=false;
            end
            Obj(nObjects,2)=i;
            nObjects=nObjects+1;
            Obj(nObjects,1)=i+1;
        end
    end
end
Obj(nObjects,2)=i+1;

% Object index list, to real connect object lines
Objects=cell(1,nObjects);
for i=1:nObjects
    % Determine if the line is closed
    if(Lines(Obj(i,1),1)==Lines(Obj(i,2),2)) 
        Objects{i}=[Lines(Obj(i,1):Obj(i,2),1);Lines(Obj(i,1),1)];
    else
        Objects{i}=Lines(Obj(i,1):Obj(i,2),1);
    end
end
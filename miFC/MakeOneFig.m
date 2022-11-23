
model = createpde;
geometryFromEdges(model,@lshapeg);
applyBoundaryCondition(model,'dirichlet', ...
                             'Edge',1:model.Geometry.NumEdges, ...
                             'u',0);
c = 1;
a = 0;
f = 1;
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);
generateMesh(model,'Hmax',0.05);
results = solvepde(model);

v = linspace(-1,1,101);
[X,Y] = meshgrid(v);
querypoints = [X(:),Y(:)]';
uintrp = interpolateSolution(results,querypoints);
% For plots
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\ColorBrewer');

uintrp = reshape(uintrp,size(X));
figure()
colormap(summer);
mesh(X,Y,uintrp)
xlabel('x')
ylabel('y')

%%
clear t 
t = 1:0.01:2;
s = sin(2*pi*t);

for jj=1:5
for ii=1:length(t)
   s(jj,ii)=s(1,ii)*(rand(1)+rand(1));
end
end

figure()
hold on
colormap(summer);
scatter(t,s,15,[0 0.4470 0.7410],'filled')
plot(t,sin(2*pi*t),'k')
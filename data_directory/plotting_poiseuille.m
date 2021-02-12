%% 0---------------0--------------------0--------------------0-----------------------0-----------0
%%% animation

%% N.B. passing a viscosity of 1 will result in Re=1, viscosity of 10 will result in Re=0.01.
dt = 0.01;       %% outputting interval
ts = 0;          %% start time 
tf = 2.0;         %% end time
visc = 1.0;      %% Viscosity
lambd = 1.0;    %% Relaxation time
bta = 0.1;       %% Viscosity ratio




ns=floor(ts/dt)+1;nf=floor(tf/dt); %% start and end frames
figure(3)
for n=ns:nf
hold off;poiseuille(dt*n,visc,lambd,bta);hold on;
fname=strcat('PART',num2str(n));A=load(fname);
plot(A(:,2),A(:,3),'rx','linewidth',2);axis([0 1 -0.1 3]); %3
set(gca,'Fontname','Times');set(gca,'Fontsize',14)
xlabel('Cross-stream coordinate, y');ylabel('Streamwise velocity, u')
text(0.05,0.5,strcat('Time=',num2str(dt*n),'s'),'Fontname','Times','Fontsize',14)
%print('-dpng','%d.png',num2str(n))
pause(0.001);end
%% 0-------------------0------------------0---------------------0---------------0------------------0

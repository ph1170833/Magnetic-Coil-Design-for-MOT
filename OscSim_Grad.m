tspan = 0:1e-4:1e-2;
clf
vcap=7; %Capture Velocity
x0 = 0; %Initial Position
curr_sm = 0.335; %Current through Coils: Gives 2.5Gauss/cm Gradient 
I = 2e3; %Intensity of Laser
for i=1:8
Current = i*curr_sm;
[t,Y] = ode45(@(t,y) TrapOsc(t,y,I,Current), tspan, [x0; vcap]);
plot(t, Y(:,1))
ylim([-1e-3 3.5e-3])
hold on
end

xlabel('Time (s)')
ylabel('Displacement from Trap Center (m)') 

 Legend=cell(8,1);
 for iter=1:8
   x=2.5*iter;
   %Legend{iter}=strcat('Intensity(mW/cm^2)=', num2str(x));
   Legend{iter}=strcat(num2str(x), 'Gauss/cm');
 end
 title(legend,'Gradient')
 legend(Legend)
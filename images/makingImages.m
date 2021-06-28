% makingImages.m
% Quick example of how to make good images with matlab.
% Written by Rodney Metoyer for the Engineering Mechanics and Space Systems
% Laboatory Research Group. Email rmmetoye@ncsu.edu with questions.

clearvars; close all; clc;

% First lets make some data
x = -pi:0.01:2*pi;
y = cos(x);

% Now lets plot those data

plot(x,y);
% When I call plot it automatically puts the XData and YData on the current
% axes and creates an axes if one does not exist. The axes are
% automatically created as a child of the current figure and if a figure
% does not exist one is created. I am mentioning this because, in general,
% the more you know what is happening when you use a built-in matlab
% function the more you will be able to use that funciton effectively. When
% I ask matlab to plot(x,y), matlab first checks to see if there is an
% instantiated figure object. If there isn't, it creates one. If there is,
% it checks that figure's CurrentAxes property to see if an axes object
% exists. If it does, it sets the Children property to a Line with the
% appropriate XData, YData, and ZData.

% You can do all of what Matlab does yourself, but why bother. Let matlab
% do the work, but know what it is doing. If you are ever unsure of what a
% built in function is doing, set a break point at teh begining of the
% funciton and then step into the function by pressing F11 when you hit the
% break point after you click run. You can also debug from the cmd prompt
% (e.g. dbstep, dbcont, etc.). For more information see:
% https://www.mathworks.com/help/matlab/matlab_prog/debugging-process-and-features.html

% You should always save the handle to the figure in a variable so that you
% can use it later. To do that, create a figure object. The constructor
% returns the handle.
hfig = figure;
% Now the figure you just made is the current figure. If you call plot now,
% matlab will put stuff on the figure that hfig is pointing to.
plot(x,y);

% You will notice that you currently have two figures open. When you create
% a new figure it gets prepended on the graphics root Children array. You
% only have the handle for one of the figures saved right now. You can
% always get handles to other figures, and really do anything with open
% graphics using the graphics root object. Type 'help groot' in the cmd
% prompt for more info.

% For now, let's work with only the figure for which we have the handle and
% ignore the other one. As an axercise, you can try to close the other
% figure from the cmd propmt by pulling its handle from the graphics root
% object. HINT: close(the last object in the root Children array)

% Right now all of the settings are the defaults. You can pass in arguments
% when you call plot that will change some settings. You can also change
% them on teh object directly. Everything (mostly) is public.

% Let's change the color of the line. The Children of the figure are what
% actually get plotted on. In this case we have a single axes object.
ax = hfig.Children;
ln = ax.Children;
ln.Color = 'r';

% Now the line color is red. In the cmd prompt type ln and you will get a
% list of all of the properties of the line (be sure to click show all).

% Let's add a legend and some labels
legend('A Line','Location','best');
xlabel('X Axis');
ylabel('Y Axis');
title('A Figure');

% Now we can save the figure to use in the powerpoint
saveas(hfig,'crappy.png');

% But this figure is kind of crappy. We can make it much better by
% configuing a bunch of stuff. However, that takes time. Luckily, someone
% wrote a wrapper to help us out. Go to
% https://github.com/altmany/export_fig and download export_fig. Once you
% have it downloaded put it in an approproate folder and add the path to
% the search path using the pathtool (help pathtool for more info).

% Let's save a higher-res version of this exact picture
export_fig(hfig,'better.png','-m5');

% Notice that export_fig saves exactly what is on the screen. This is great
% because saveas tries to guess at what you want to save. Sometimes that's
% cool, like when it changes the figure color to white to match the axes,
% but sometimes (usually) it guesses wrong. With export_fig you get exactly
% what you tell it to give you. if you don't tell it anything, you get the
% matlab defaults.

% Let's make the figure clear so that we can put stuff under it in the
% powerpoint. Unfortunately, export_fig has some faults. You can tell it to
% export a figure with a tranparent background, but it won't make the axes
% and other Children transparent. That's ok, we can do that manually.
ax.Legend.Color = 'none';
ax.Color = 'none';
export_fig(hfig,'good','-m5','-transparent','-png','-eps','-pdf');

% Now you have a high resolution, transparent image that you can use in a
% powerpoint. See makingImages.pptx if you have access to it. Notice that I
% also saved two vector formats of the image. The vectors are what you
% would use for publication. 

% You can do a lot with export_fig. See the documentation on github.
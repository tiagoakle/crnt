load NetworkSizeVsItCountSingleLink
nsize = 100:200:5000;
avg = mean(Results);
h=plot(nsize,avg);
xlabel('Number of complexes in network.');
ylabel('Average number of iterations until convergence.');
title('Single class network, average iteration count vs network size');
saveas(h,'SingleNetAvgIterationsVsNetSize.jpg');
load NetworkSizeVsItCountMultipleLink
avg = mean(Results);
h=plot(nsize,avg);
xlabel('Number of complexes in network.');
ylabel('Average number of iterations until convergence.');
title('Two linkage class network, average iteration count vs network size');
saveas(h,'MultipleNetAvgIterationsVsNetSize.jpg');
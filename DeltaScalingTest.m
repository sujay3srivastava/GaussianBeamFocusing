figure;
Deltasort1 = sort(Delta1);
plot(Deltasort1);
hold on;
Deltasort2 = sort(Delta2);
plot(Deltasort2);

xlabel('Meta-Unit points');
ylabel('Delta Values(Post-Scaling');
legend('Delta1 -100um taper,24um BeamWidth','Delta2 -100um taper,24um BeamWidth')

% where are these points located

% Two decisions:
% 1) Either stick with lower size delta and points where max delta are
% there are kept
% 2) Or do processing and remove points based on Optimization criteria
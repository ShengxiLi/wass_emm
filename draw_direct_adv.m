function draw_direct_adv(samples, theta, colorn)
    range = max(max(samples));
    [M,L] = size(theta);
    for i = 1:L
        endpoint = range*theta(:,i);
        x_ax = linspace(-endpoint(1),endpoint(1),1e5);
        y_ax = linspace(-endpoint(2),endpoint(2),1e5);
        plot(x_ax,y_ax,'color',colorn,'LineWidth',5);hold on
    end
load('2.mat'); 
y = y(:); 
u = x(:); 

orders = [2, 3, 4]; 
step_sizes = [0.02]; % Step sizes for LMS algorithm aka mu
num_epochs = length(u); % Number of iterations
N = length(u); 

for p = orders
    figure;
    hold on;
    for mu = step_sizes
        [weights, errors] = lms_fin(u, y, mu, p, num_epochs);
        for i = 1:p
            plot(weights(i, :), 'LineWidth', 0.5, 'DisplayName', ['w[' num2str(i-1) '] for \mu = ' num2str(mu)]);
        end
    end
    hold off;
    title(['Adaptive Filtering via LMS with p = ' num2str(p)], 'FontSize', 14);
    xlabel('Iteration', 'FontSize', 12);
    ylabel('Filter Coefficients', 'FontSize', 12);
    legend('show', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);
    
    % Ta coefficients tou filtrou
    wiener_coeff = weights(:, end);
    disp(['Wiener coefficients for p = ' num2str(p) ':']);
    disp(wiener_coeff');
end

for p = orders
    for mu = step_sizes
        [weights, errors] = lms_fin(u, y, mu, p, num_epochs);
        

        %gia na do an tautizetai to reconstructed shma vazontas to filtro sthn eisodo me thn eksodo 
        output_signal = zeros(length(u), 1);
        for n = p:length(u)
            x = u(n:-1:n-p+1);
            output_signal(n) = weights(:, end)' * x;
        end
        
        figure;
        plot(y, 'r', 'LineWidth', 1.5, 'DisplayName', 'y (Desired Signal)');
        hold on;
        plot(output_signal, 'b', 'LineWidth', 1.5, 'DisplayName', 'w*u (Filtered Signal)');
        plot(u, 'g', 'LineWidth', 1.5, 'DisplayName', 'u (Input Signal)');
        hold off;
        title(['Adaptive Filtering via LMS with p = ' num2str(p) ', \mu = ' num2str(mu)], 'FontSize', 14);
        xlabel('Sample Index', 'FontSize', 12);
        ylabel('Signal Value', 'FontSize', 12);
        legend('show', 'Location', 'best');
        grid on;
        set(gca, 'FontSize', 12);
    end
end


% LMS synarthsh symfwna me tis diafaneies
function [w, mse] = lms_fin(input, output, mu, p, N)
    % Elegxos gia to an p kai N einai >= 2
    if p >= 2 && N >= 2
        % Orismos tou x san u kai tou d san y
        x = input; 
        d = output;
        
        % Arxikopoiisi ton metavliton
        w = zeros(p, N); % Pinakas varon
        w_after_update = zeros(p, 1); % Trexousa timi varon
        e = zeros(N, 1); %pinakas gia th diafora
        
        % loopa gia ton ypologismo ton varon me ti methodo LMS
        for n = p:N
            % Ypologismos tou dianysmatos eisodou
            signal = x(n:-1:n-p+1);
            % Eksodos tou systimatos
            y(n) = w_after_update' * signal;
            % Ypologismos diaforas
            e(n) = d(n) - y(n);
            % Enimerosi varwn me ti methodo LMS symfwna me tis diafaneies
            w_after_update = w_after_update + mu * conj(signal) * e(n);
            % Apothikeusi twn varwn
            w(:, n) = w_after_update;
        end
        
        % Ypologismos tou Mesou tetragonikou sfalmatos
        mse = mean(e(p:N).^2);
        fprintf('The Mean Squared Error is: %f with p %d and mu %f \n', mse, p, mu);
    else
        error('p should be >=2');
    end
end


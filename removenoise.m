data = load('3.mat');
y = data.y;
u = data.u;
N  = length(y);
n  = 0:(N-1);                           
v2 = u;                                             

%gia kathe order kalo th synarthsh
orders = [4,15,30];
for i = 1:length(orders)
    current_M = orders(i); 
    [y_func, w, en] = noise_cancel(y, current_M, v2);   
    
    figure('Position', [100, 100, 800, 600]);
    plot(n, en, 'b', 'LineWidth', 1.5); 
    grid on; % Add grid lines
    xlabel('n', 'FontSize', 12);
    ylabel('Signal', 'FontSize', 12); 
    title(['M = ', num2str(current_M)], 'FontSize', 14); 
end



function[r]=biased_ac_func(x,lags)
    N=length(x);
    for m=1:lags
        for n=1:N+1-m
            x1(m,n)=x(n-1+m);
        end;
    end;
        r1=x1*x';
        r=r1'./N;
end

function[y_func,w,en]=noise_cancel(y,M,v2)
    xn = y;                               
    % disp(['Shape of v2: ', num2str(size(v2))]);
    % disp(['Shape of M: ', num2str(size(M))]);
    %autocorrelation
    v2autoc=biased_ac_func(v2',M);
    Rv2 = toeplitz(v2autoc);
    R   = Rv2;
    %eterosysxetisei
    rd=xcorr(xn,v2,'biased');
    % fprintf('len v2 is %d and len v2 + M - 1 is %d\n', length(v2'), length(v2') + M - 1);
    % disp(['Shape of Rd before: ', num2str(size(rd))]);
     rd = rd(length(v2'):length(v2')+M-1);
    % disp(['Shape of R: ', num2str(size(R))]);
    % disp(['Shape of rd: ', num2str(size(rd))]);
    w=R\(rd);                                                                
    y_func=filter(w,1,v2);                      
    en=xn-y_func;            
end




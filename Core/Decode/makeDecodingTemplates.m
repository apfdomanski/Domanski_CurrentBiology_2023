function templates = makeDecodingTemplates(tb,times,bw,sigma,shoulder,plotYN)
    
if nargin<5
    sigma    = 4;    % kernel standard deviation (s)
    shoulder = 2;    % Time ranges to draw kernel over (sigmas)
    plotYN = false;
end

vector = zeros(length(times),length(tb));
templates.kernel = cell(length(times),1);


if length(sigma)>1
    for i=1:length(times)
        % prepare a kernel:
        edges            =   -shoulder*sigma(i)  :  bw  :  shoulder*sigma(i);
        templates.kernel{i} =   zeros(size(edges));
        templates.kernel{i} =   normpdf(edges, 0, sigma(i));
        % templates.kernel = templates.kernel*kernal.bw*1E-3; % Multiply by bin width so the probabilities sum to 1
        templates.kernel{i}=mat2gray(templates.kernel{i}); % Scale between 0-1
        kernal_center = ceil(length(edges)/2);
        
        idx = find(min(abs(times{i}-tb))==abs(times{i}-tb),1,'first'); % closest bin to peak;
        vector(i,idx) = 1;
        temp = conv(vector(i,:),templates.kernel{i});
        % Trim out the relevant portion of the convolved template
        templates.result(i,:) = temp(kernal_center:length(tb)+kernal_center-1);
    end
     
else
    % prepare a kernel: 
    edges            =   -shoulder*sigma  :  bw  :  shoulder*sigma;
    templates.kernel =   zeros(size(edges));
    templates.kernel =   normpdf(edges, 0, sigma);

    % templates.kernel = templates.kernel*kernal.bw*1E-3; % Multiply by bin width so the probabilities sum to 1
    templates.kernel=mat2gray(templates.kernel); % Scale between 0-1
    % plot(edges*bw,templates.kernel)

    % Find the index of the kernel center
    kernal_center = ceil(length(edges)/2);
    
    for i=1:length(times)
        %idx = find(ismember(tb,times{i})); % bin at centre of peak
        idx = find(min(abs(times{i}-tb))==abs(times{i}-tb),1,'first'); % closest bin to peak;
        vector(i,idx) = 1;
        temp = conv(vector(i,:),templates.kernel);
        
        % Trim out the relevant portion of the convolved template
        templates.result(i,:) = temp(kernal_center:length(tb)+kernal_center-1);
    end
end



templates.bw       = bw;
templates.shoudler = shoulder;
templates.sigma    = sigma;
templates.tb       = tb;

if plotYN
    figure; hold on
    for i=1:length(times)
        plot(tb,i*1.2+templates.result(i,:),'LineWidth',1.5,'color','b')                    
    end
    axis tight
    xlabel('time (s)')
    ylabel('template prototypes')

end

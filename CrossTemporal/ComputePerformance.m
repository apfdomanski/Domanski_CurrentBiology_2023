function [performance] = ComputePerformance(train_data,test_data,training_label,testing_label)
% class is the predicted target label from the test_data
[class,~,~,~,~] = classify(test_data,train_data,training_label);
fitcdiscr
% Checking how different the predicted label is from the actual label for
% all the test trials
perf = class-testing_label';
% Extracting all the correctly predicted target label
perf_ind = find(perf==0);
% The performance of the decoder is computed as the percentege of number of correct
% predictions in the decoding.
performance = length(find(perf==0))*100/length(testing_label);
end
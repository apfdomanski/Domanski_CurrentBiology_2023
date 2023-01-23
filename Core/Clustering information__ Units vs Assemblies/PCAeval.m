function [ScoreTrain,CoeffTrain,ScoreEval,CoeffEval,explainedTrain,explainedTest] =PCAeval(X,X_)

% centre the data around zero using shared mean
mu = mean(X);
% mu = mean([X;X_]);
X  = X-mu;
X_ = X_-mu;

A = (X'*X) / length(X); % compute the covariance matrix (normalised by the num of elements)
[V,latent] = eig(A);    % compute the eigenvectors V and eigenvalues 'latent' -- this results in a increasing ordering
explained = (flipud(100*latent/sum(latent)));
V = fliplr(V);          % flip the eigenvectors so that the most significant come first
% V = V(:,1:D);         % take only those eigenvectors required
Y  = X  * V;            % project the original data into the new coordinate system
Y_ = X_ * V;            % project the unseen data to the new coordinate system
V_ = Y' /  X_';


figure
subplot(1,2,1); hold on
plot(staggerplot(Y(:,1:5),0, 10),'b')
subplot(1,2,2); hold on
plot(staggerplot(Y_(:,1:5),0, 10),'r')

CoeffTrain = V;
ScoreTrain = Y;

CoeffEval = V_;
ScoreEval = Y_;

explainedTrain = explained;


A_ = (X_'*X_) / length(X_); % compute the covariance matrix (normalised by the num of elements)
[V_,latent_] = eig(A_);    % compute the eigenvectors V and eigenvalues 'latent' -- this results in a increasing ordering
explainedTest = (flipud(100*latent_/sum(latent_)));
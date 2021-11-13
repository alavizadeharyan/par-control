tic
%%%%%%%%%%%%%%%%%% Initializations %%%%%%%%%%%%%%%%%%%
N = 201;      % Mesh size of the whole interval + 1
M = 6;        % Number of parallel CPUs + 1
interval_start = 0;
interval_end = 10;
A = [1 0; 0 2];
B = [0 1; 1 0];
y0 = [1;0];
yref = [0;1];
alpha = 1e1;
iterations = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = interval_end - interval_start;
L = (N-1)/(M-1);  % Mesh size in each parallel process.
t = linspace(interval_start, interval_end, N);
tl = linspace(interval_start, interval_end, M);
dt = t(2) - t(1);
betal = T / (tl(2)-tl(1));
alphal = alpha / betal;
c = zeros(1,N-1);
C = zeros(L,M-1);
J = 0;      % Cost function initialized by 0!

% Calculation of initial p and y on the whole interval:
y = zeros(2,N);
y(:,1) = y0;
for n = 1:N-1
  y(:,n+1) = expm(i*A*dt)*expm(i*B*c(n)*dt)*y(:,n);
end

p = zeros(2,N);
p(:,N) = -yref;
for n = N:-1:2
  p(:,n-1) = expm(-i*B*c(n-1)*dt)*expm(-i*A*dt)*p(:,n);
end

%%%%%%%%%%%%%%%%% Starting iterations %%%%%%%%%%%%%%%%%%
for m = 1:iterations
  
  % Calculation of lambdas associated to c:
  lambda = zeros(2,M);
  lambda(:,1) = y(:,1);
  for j=2:M
    lambda(:,j)= (1-(tl(j)/T))*y(:,(j-1)*L+1)-(tl(j)/T)*p(:,(j-1)*L+1);
  end

  % The parallelization process is done here:
  C = reshape(c, L, M-1);
  parfor j=1:M-1
    C(:,j) = monotonic(tl(j), tl(j+1), lambda(:,j), lambda(:,j+1), alphal, C(:,j), A, B, L+1);
  end
  c = reshape(C, 1, N-1);
  
  % Updating y and p on the whole interval:
  for n=1:N-1
    y(:,n+1)=expm(i*A*dt)*expm(i*B*c(n)*dt)*y(:,n);
  end

  for n=N:-1:2
    p(:,n-1)=expm(-i*B*c(n-1)*dt)*expm(-i*A*dt)*p(:,n);
  end
  
  % Updating cost function:
  J = -real(y(:,end)'*(yref))+.5*alpha*dt*c*c';
  fprintf(2,'Iteration %i | J = %f\n',m,J);
  
end
toc

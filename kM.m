%ML k-Means algorithm with insertion sorting in Octave.
%The code runs multiple times (here again = 3). Each time the result may be different, since the algorithm strongly depends on starting points (which are random here). 
%in this code I chose 3 starting points, which further result in 3 clusters being expanded
%my dataset was CTG.xls - https://archive.ics.uci.edu/ml/datasets/Cardiotocography - i chose the NSP to be the objective variable (last column: Normal=1; Suspect=2; Pathologic=3).
%Remember to clean the data from NaNs, empty cells and headings - like in the CTGdata100.xlsx
%k-Means works best with equally numerous and strongly dispersed groups of data

for again = 1:3 %how many times run the code again (do the clustering till reaching the stop condition(s))
  %DATA LOADING
  clear all % remove old variables
              tic
  dataraw = xlsread("CTGdata100.xlsx"); %pkg load io (In command window. Also select the directiory where you have the input file)
  data = dataraw(:,1:21);
  lengthR = length(data(:,1));
  lengthC = length(data(1,:));
  permutt = randperm(lengthR); %permutation for randomization
  
  %WAITBAR
  wb = waitbar (0,sprintf("Iteration: %d  ", 0,0));

  %PICK 3 RANDOM INDEXES and make sure they belong to different classes (not necessarily - for best performance). 
  %These 3 points are the starting points (centroids). They heavily determine the algorithm's effectiveness.
  for p=1:length(permutt)
    if dataraw(permutt(p),end) == 1
      cMean1 = data(permutt(p),:);
      cindx1 = permutt(p);
      break %break the action if found the record with desired class
    endif
  endfor
  for p=1:length(permutt)
    if dataraw(permutt(p),end) == 2
      cMean2 = data(permutt(p),:);
      cindx2 = permutt(p);
      break
    endif
  endfor
  for p=1:length(permutt)
    if dataraw(permutt(p),end) == 3
      cMean3 = data(permutt(p),:);
      cindx3 = permutt(p);
      break
    endif
  endfor
  cindxs = [cindx1, cindx2, cindx3]; %Centroids Indexes
  cMeans = [cMean1; cMean2; cMean3]; %Centroids Means
  %PICKING CENTROIDS - OVER
  
  maxidvardiff = [4200]; %just declaring some higher value
  itrtmax = 100;     % MAXIMALLY itrtmax ITERATIONS TO STOP THE CLUSTERING
  for itrt = 1:itrtmax
    waitbar (itrt/itrtmax,wb,sprintf("Iteration %d of maximally %d iterations", itrt, itrtmax));  
  
      %INDEXES MATRIX
      for j=1:lengthR
        for i=1:length(cindxs(1,:))
          mtrxIndxs(j,i) = i;
        endfor
      endfor
      
      %CALCULATE THE DISTANCES
      differsq(lengthR,length(cindxs(1,:))) = 0; %just zeros
      for r=1:lengthR
        for c=1:length(cindxs(1,:))
          for i=1:lengthC
            differsq(r,c) = differsq(r,c) + (data(r,i) - cMeans(c,i))^2; %0 + difference squared
          endfor
        endfor
      endfor
      dists(lengthR,length(cindxs(1,:))) = 0; %just zeros
      for r=1:lengthR
        for j=1:length(cindxs(1,:))
          dists(r,j) = sqrt(differsq(r,j)); %distances
        endfor
      endfor

      %INSERTION SORT of 2 arrays
      for j=1:length(dists(:,1))
        for i=2:length(dists(1,:))
          valD = dists(j,i);
          valI = mtrxIndxs(j,i);
          while i>1 && valD < dists(j,i-1);
            dists(j,i) = dists(j,i-1);
            mtrxIndxs(j,i) = mtrxIndxs(j,i-1);
            i--;
          endwhile
          dists(j,i) = valD;
          mtrxIndxs(j,i) = valI;
        endfor
      endfor

      % \CLUSTERS - cluster's(1/2/3) Indexes
        clear cIndxs1 cIndxs2 cIndxs3;
      for r=1:lengthR
        if mtrxIndxs(r,1) == 1
          cIndxs1(end+1)=r;
        elseif mtrxIndxs(r,1) == 2
          cIndxs2(end+1)=r;
        elseif mtrxIndxs(r,1) == 3
          cIndxs3(end+1)=r;
        endif
      endfor
      
      % CLUSTERS' MEANS AND VARIANCES
      cMean1 (itrt,:) = mean(data(cIndxs1,:));
      cMean2 (itrt,:) = mean(data(cIndxs2,:));
      cMean3 (itrt,:) = mean(data(cIndxs3,:));
      cMeans = [cMean1(itrt,:); cMean2(itrt,:); cMean3(itrt,:)]; 
      cVar1  (itrt,1) = var(data(cIndxs1));
      cVar2  (itrt,1) = var(data(cIndxs2));
      cVar3  (itrt,1) = var(data(cIndxs3));
      
      %THE CHANGE OF VARIATION BETWEEN ITERATIONS itrt
      if itrt > 1
        dvar(itrt,:) = [abs(cVar1(itrt-1,1)-cVar1(itrt,1)), abs(cVar2(itrt-1,1)-cVar2(itrt,1)), abs(cVar3(itrt-1,1)-cVar3(itrt,1))];
        maxidvar(itrt,1) = max(dvar(itrt,:));
        maxidvardiff(itrt,1) = abs (maxidvar(itrt,1) - maxidvar(itrt-1,1));
      endif
    if maxidvardiff (itrt,1) == 0       % STOP CONDITION (stop when variation between iterations is equal 0 - clusters don't change)
    break
    endif
  endfor
  delete(wb)
  
  %VERIFY
  err = 0;
  for i=1:length(cIndxs1)
    if dataraw(cIndxs1(i),end) != 1
    err++;
    endif
  endfor
  for i=1:length(cIndxs2)
    if dataraw(cIndxs2(i),end) != 2
    err++;
    endif
  endfor
  for i=1:length(cIndxs3)
    if dataraw(cIndxs3(i),end) != 3
    err++;
    endif
  endfor
  
  EFF = 100*((length(data)-err)/length(data)); %effectiveness
  printf("Result in %d iterations, effectiveness %d%%\n", itrt, EFF)
              toc
endfor
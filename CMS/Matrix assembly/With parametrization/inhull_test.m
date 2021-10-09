count=0;
rng('default')
for k=1:100
     indata.theta_k=normrnd(1,.1,4,1); % change the sample point 
     checkhull=check_inhull(indata);
     count=count+checkhull;
    
    
end
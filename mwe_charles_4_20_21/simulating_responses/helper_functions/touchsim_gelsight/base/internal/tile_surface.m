function tdat = tile_surface(dat,xy_size)

if nargin<2
    xy_size = [1000 1000];
end

% make square
[x,y] = size(dat);
if x~=y
    tp = false;
    if x<y
        dat = dat';
        [x,y] = size(dat);
        tp = true;
    end
    diff = x-y;
    dat = dat(1+floor(diff/2):end-ceil(diff/2),:);
    if tp
        dat = dat';
    end
end

% find smallest square without NaNs
if any(isnan(dat(:)))
    center = size(dat,1)/2;
    d1nan = find(isnan(diag(dat)));
    d2nan = find(isnan(diag(dat')));
    dnan = union(d1nan,d2nan);
    dmin = min(abs(center - dnan))-2;
    dat = dat(ceil(center-dmin):floor(center+dmin),ceil(center-dmin):floor(center+dmin));
    dat(isnan(dat)) = 0;
end

% align shifted pattern
if xy_size(1)>size(dat,1) || xy_size(2)>size(dat,1)
    max_overlap = ceil(0.1*size(dat,1));
    
    % find best vertical overlap
    for i=50:max_overlap
        cc_tmp = corrcoef(reshape(dat(1:i,:),[],1),reshape(dat(end-i+1:end,:),[],1));
        cc(i) = cc_tmp(2);
    end
    [~,v_shift] = max(cc);
    
    % find best horizontal overlap
    for i=50:max_overlap
        cc_tmp = corrcoef(reshape(dat(:,end-i+1:end),[],1),reshape(dat(:,1:i),[],1));
        cc(i) = cc_tmp(2);
    end
    [~,h_shift] = max(cc);
    
    v_overlap = dat(1:v_shift,:);
    h_overlap = dat(:,1:h_shift);
    
    dat(end-v_shift+1:end,:) = dat(end-v_shift+1:end,:)+v_overlap;
    dat(end-v_shift+1:end,:) = dat(end-v_shift+1:end,:)/2;
    dat(:,end-h_shift+1:end) = dat(:,end-h_shift+1:end)+h_overlap;
    dat(:,end-h_shift+1:end) = dat(:,end-h_shift+1:end)/2;
    
    dat = dat(v_shift+1:end,h_shift+1:end);
    
    v_reps = ceil(xy_size(1)/size(dat,1));
    h_reps = ceil(xy_size(2)/size(dat,2));
    
    dat = repmat(dat,v_reps,h_reps);
end

tdat = dat(1:xy_size(1),1:xy_size(2));

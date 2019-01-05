function imgFilenamePattern = recoverFilenameFromPattern(imgFilenamePattern, frame)



for ii = 10:-1:1
    
    pp = strfind(imgFilenamePattern,repmat('?',[1 ii]));
    precis = ['%.' num2str(ii) 'd'];
    for kk = 1:length(pp)
        imgFilenamePattern(pp(kk):pp(kk)+ii-1) = num2str(frame,precis);
    end
    
end
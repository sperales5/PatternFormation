list = {'GM Model','Allen Cat Model','Allen Snake Model'}; 
[indx,tf] = listdlg('ListString',list);

if ismember(1,indx)
    disp('GM Model Selected')
end

if ismember(2,indx)
    disp('Allen Cat Model Selected')
end

if ismember(3,indx)
    disp('Allen Snake Model Selected')
end
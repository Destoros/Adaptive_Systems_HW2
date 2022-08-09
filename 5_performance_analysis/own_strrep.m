function new_string = own_strrep(f)
       
    new_string = strrep(f.Title.String,' ', '_');
    new_string = strrep(new_string,'(', '');
    new_string = strrep(new_string,')', '');     
    new_string = strrep(new_string,'\', '');
    new_string = strrep(new_string,'.', '_');

end
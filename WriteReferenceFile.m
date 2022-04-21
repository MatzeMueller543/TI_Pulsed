function WriteReferenceFile(Reference_10min,FileName)

Table2Write 	= table(cellstr(Reference_10min.Timestamp),Reference_10min.LOS_TI_N',Reference_10min.LOS_TI_S',...
            'VariableNames',{'Timestamp','LOS_TI_North','LOS_TI_South'});
writetable(Table2Write,FileName)

end




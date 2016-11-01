function [P_info, Q1, Map] = calculate_P_info(Q0, cutoff)

[NumSeq NumPos] = size(Q0);
Q1 = Q0;

cutoff_freq = (NumSeq-cutoff)/NumSeq;

P_info = repmat(0, [6 NumPos]);
Map = [];
A65 = repmat(65, [NumSeq 1]);
A67 = repmat(67, [NumSeq 1]);
A71 = repmat(71, [NumSeq 1]);
A84 = repmat(84, [NumSeq 1]);
A45 = repmat(45, [NumSeq 1]);

amb = [];
amb_count = [];
for i = NumPos : -1: 1
    
   P_info(1, i) = sum( A65 == Q1(:, i));            %A 
   P_info(2, i) = sum( A67 == Q1(:, i));            %C
   P_info(3, i) = sum( A71 == Q1(:, i));            %G
   P_info(4, i) = sum( A84 == Q1(:, i));            %T
   P_info(5, i) = sum( A45 == Q1(:, i));            %gaps
   P_info(6, i) = max( P_info(1:5, i))/NumSeq;      %conservation  
   
   if sum(P_info(1:5, i)) < NumSeq
     amb = [amb i];  
     amb_count = [amb_count NumSeq-sum(P_info(1:5, i))];
     [ i]
   end  
   if sum(P_info(1:5, i)) < NumSeq || P_info(6, i) >= cutoff_freq || P_info(5, i) > 0 
      Q1(:, i) = [];
     % [i P_info(6, i) P_info(5, i)]
   else 
      Map = [i Map];
   end     
end   
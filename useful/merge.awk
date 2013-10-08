BEGIN{ FS = "\t";}
{
 # populate array with customer as index and sum of balances as values
 balance[$1] = balance[$1] + $2
}
END{ # sort the indices of the array and put them in array indices
     n = asorti(balance,indices)

     # print the sorted array
     for (i = 1; i <= n; i++)
         printf "%s\t%d\n", indices[i], balance[indices[i]]
}

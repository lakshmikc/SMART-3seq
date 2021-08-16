import os
import csv
import glob


## after some google https://mail.python.org/pipermail/tutor/2004-November/033475.html
## The idea is to keep the count column into a list.


files = glob.glob('*-gene-counts.txt')
files.sort(key=lambda f: int(filter(str.isdigit, f)))

list_column = []
n = 1
for file in files:

    print n,file
    column_data = []
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter = "\t")
            # skip the comment line
        comment = next(reader)
        if n <= 1:
            for row in reader:
                # for the first file, keep the gene column as well
                column_data.append(row[0] + '\t' + row[7])
        else:
            for row in reader:
                column_data.append(row[7])
        n = n + 1
    list_column.append(column_data)


# This creates a list of row lists from the list of column lists
# If any of the column lists are too short they will be padded with None
# map function is a gem :)
rows = map(None, *list_column)

with open('PR1643-mm-genes-counts.txt','w') as f_out:
     for row in rows:
         f_out.write('\t'.join(row))
         f_out.write('\n')



files = glob.glob("*-exon-counts.txt")
files.sort(key=lambda f: int(filter(str.isdigit, f)))

list_column = []
n = 1
for file in files:
    print n,file
    column_data = []
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter = "\t")
            # skip the comment line
        comment = next(reader)
        if n <= 1:
            for row in reader:
                # for the first file, keep the gene column as well
                column_data.append(row[0] + '\t' + row[7])
        else:
            for row in reader:
                column_data.append(row[7])
        n = n + 1
    list_column.append(column_data)


# This creates a list of row lists from the list of column lists
# If any of the column lists are too short they will be padded with None
# map function is a gem :)
rows = map(None, *list_column)

with open('PR1643-mm-exons-counts.txt','w') as f_out:
     for row in rows:
         f_out.write('\t'.join(row))
         f_out.write('\n')








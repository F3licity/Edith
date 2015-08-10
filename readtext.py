#Compiling
import csv

#This way you can write in a new file that you can open with excel
ifile = open('eftychia.csv', 'rb')
reader = csv.reader(ifile)
ofile = open('test.csv','wb')
writer = csv.writer(ofile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_ALL)
#So basicly is created a copy of eftychia.csv into test.csv

for row in reader:
    writer.writerow(row)
    print(row)
ifile.close()
ofile.close()
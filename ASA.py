import csv
import math
import matplotlib.pyplot as plt
s=input('Enter path name for ASA folder: ')
z=input('Enter path name for output folder: ')
z=z.replace('"', '')
s=s.replace('"', '')
amino, val, l, l1, avg, mvav=[[], [], []], [[], [], []], [], [], [], []
window_s=12
for a in range(3):
    a1=str(a+1)
    name=f'{s}/{a1}.txt'
    f=open(name, 'r+')
    s1=f.readlines()
    for i in range(0, len(s1), 4):
        s1[i+1]=s1[i+1].replace('\n', '')
        s1[i+1]=s1[i+1].strip()
        s1[i+2]=s1[i+2].replace('\n', '')
        s1[i+2]=s1[i+2].strip()
        for u in s1[i+1].split(' '):
            if u!='':
                l.append(u)
        for u in s1[i+2].split(' '):
            if u!='':
                l1.append(u)
        for j,k in zip(l, l1):
            amino[a].append(j)
            val[a].append(k)
        l, l1=[], []
    f.close()
p=len(amino[0])

for i in range(p):
    sum1=0
    for j in range(3):
        sum1=sum1+int(val[j][i])
    mean=sum1/3
    avg.append(mean)

with open(f'{z}/records.csv', 'a+', newline='') as csvFile:
    writer=csv.writer(csvFile)
    writer.writerow(['Position', '1st Seq.', '2nd Seq.', '3rd Seq.', 'Average'])
    writer.writerow(['', '', '', '', ''])
    for i in range(p):
        writer.writerow([i+1, val[0][i], val[1][i], val[2][i], avg[i]])
csvFile.close()

for i in range(p-(window_s-1)):
    q=avg[i:(i+window_s)]
    mvav.append((sum(q)/window_s))
p1=len(mvav)
    
with open(f'{z}/moving-average_coord.csv', 'a+', newline='') as csvFile1:
    writer=csv.writer(csvFile1)
    for i in range(p1):
        writer.writerow([i+1, mvav[i]])
csvFile1.close()

diff, relation, pr1=[], [], []
with open(f'{s}/coord.csv', 'r', newline='') as csvFile:
    reader=csv.reader(csvFile)
    print(reader)
    i=0
    for row in reader:
        pr=int(row[1])
        pr1.append(int(row[1]))
        diff.append(abs(mvav[i]-pr))
        #relation.append((abs(mvav[i]-pr)*abs(mvav[i]-pr)*abs(mvav[i]-pr))-((8*(mvav[i]+pr)*(mvav[i]+pr)*(mvav[i]+pr))/27))
        #relation.append((mvav[i]*mvav[i]*abs(mvav[i]-pr)*abs(mvav[i]-pr))/math.sqrt(pr))
        #relation.append((mvav[i]*abs(mvav[i]-pr))/pr)
        #relation.append(abs(mvav[i]-pr)/(mvav[i]*pr))
        #relation.append(abs(mvav[i]-pr)*(mvav[i]+(1/pr)))
        #relation.append(abs(math.pow(mvav[i], 1.05) - math.sqrt(pr)))
        #relation.append((mvav[i]-pr)/pr)
        #relation.append((mvav[i]-pr)/mvav[i])
        #relation.append((mvav[i]-pr)/(mvav[i] * pr))
        a1=math.sqrt(math.pow(mvav[i], 2) + math.pow(1/pr, 2))
        b1=math.sqrt(math.pow(1/pr, 2) + math.pow(mvav[i]-pr, 2))
        c1=math.sqrt(math.pow(mvav[i]-pr, 2) + math.pow(mvav[i], 2))
        sep1=(a1+b1+c1)/2
        relation.append(math.sqrt(sep1*(sep1-a1)*(sep1-b1)*(sep1-c1)))
        i+=1
csvFile.close()

with open(f'{z}/final_plotting_coord.csv', 'a+', newline='') as csvFile:
    writer=csv.writer(csvFile)
    for i in range(len(diff)):
        writer.writerow([i+1, diff[i], relation[i]])
csvFile.close()
sum2=sum(pr1)
mean1=sum2/len(pr1)
plt.figure()
plt.plot([x for x in range(1, len(diff)+1)], diff, 'ro', markersize=0.6)
plt.savefig(f'{z}/difference_plot.png')

plt.figure()
plt.plot([x for x in range(1, len(relation)+1)], relation, 'ro', markersize=0.6)
plt.savefig(f'{z}/relation_plot.png')

relation1=relation[:]
rel=sorted(relation1, reverse=True)
g=open(f'{z}/best_results.txt', 'a+')
coun=0
for i in range(len(rel)):
    p=relation1.index(rel[i])
    relation1[p]=-1
    if pr1[p]<(mean1):
        coun+=1
        g.write(str(coun) + '    ' + str(p+1) + '    ' + str(rel[i]) + '    ' + str(pr1[p]))
        for a in range(3):
            s2=''
            for l in range(window_s):
                s2+=amino[a][p+l]
            g.write('   ' + s2)
        g.write('\n')
g.close()
print(relation1)
max1=rel[0]
start=relation.index(rel[0])
for a in range(3):
    s2=amino[a][start: start+window_s]
    print(s2)

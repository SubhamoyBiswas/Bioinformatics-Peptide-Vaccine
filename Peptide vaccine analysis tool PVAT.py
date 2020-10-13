import tkinter as tk
import tkinter.ttk
import matplotlib.pyplot as plt
import csv
import math
import os
win=tk.Tk()
win.title('PVAT')
win.resizable(0, 0)
win.configure(background='white', cursor='top_right_corner')
win.geometry('685x500')
style=tkinter.ttk.Style()
style.theme_use('default')
style.configure("black.Horizontal.TProgressbar", background='green')
bar = tkinter.ttk.Progressbar(win, length=220, style='black.Horizontal.TProgressbar')
main_path, out_path='', ''
window_s=12
window_s1=12

def asa_find(s, z, folder_len):
    amino, val, l, l1, avg, mvav, folder_len=[[], [], []], [[], [], []], [], [], [], [], str(folder_len)
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
    with open(f'{z}/{folder_len}/coord.csv', 'r', newline='') as csvFile:
        reader=csv.reader(csvFile)
        print(reader)
        i=0
        for row in reader:
            pr=float(row[1])
            pr1.append(float(row[1]))
            diff.append(abs(mvav[i]-pr))
            x1=float(mvav[i])
            y1=float(1/pr)
            z1=float(mvav[i]-pr)
            a1=math.sqrt(y1*y1 + z1*z1 + z1*y1)
            b1=math.sqrt(z1*z1 + x1*x1 + z1*x1)
            c1=math.sqrt(x1*x1 + y1*y1 + x1*y1)
            s11=(a1+b1+c1)/2
            s=math.sqrt(s11*(s11-a1)*(s11-b1)*(s11-c1))
            relation.append(s)
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
    sts=[]
    rel=sorted(relation, reverse=True)
    with open(f'{z}/best_results.csv', 'a+', newline='') as csvFile:
        writer=csv.writer(csvFile)
        writer.writerow(['Rank', 'Start', 'Score', 'PV', 'Peptide 1', 'Peptide 2', 'Peptide 3'])
        for i in range(len(rel)):
            p=relation.index(rel[i])
            relation[p]=-1
            if pr1[p]<=mean1 and mvav[p] >= pr1[p]:
            #if mvav[p] >= pr1[p]:
                for a in range(3):
                    s2=''
                    for l in range(window_s):
                        s2+=amino[a][p+l]
                    sts.append(s2)
                writer.writerow([str(i+1), str(p+1), str(rel[i]), str(pr1[p]), sts[0], sts[1], sts[2]])
                sts=[]
    csvFile.close()
    max1=rel[0]
    print(len(mvav))
    bar['value']=100
    tk.Label(win, text='Complete! See your output folder', background='white', foreground='blue').grid(row=25)
    
def func_call():
    main_path=e1.get()
    asa_path=e2.get()
    out_path=e3.get()
    main_path=main_path.replace('"', '')
    asa_path=asa_path.replace('"', '')
    out_path=out_path.replace('"', '')
    length=retrieval(main_path, out_path, asa_path)
    if length!=0:
        asa_find(asa_path, out_path, length)
    else:
        tk.Label(win, text='ERROR!! NO REFSEQ FOUND', background='white', foreground='red', font=('Verdana', '20', 'bold')).grid(row=27)
    
def retrieval(main_path, out_path, asa_path):
    f=open(main_path, 'r')
    tk.Label(win, text='*********** Processing...(0%) ***********', background='white', foreground='blue').grid(row=25)
    s=f.readline()
    s1=''
    c, h, r, v, e, w=-1, 0, [], [], '', []
    while s!='':
        if s[0]=='>':
            r.append(s.split(' ', 1))
            r[h][0]=r[h][0].replace('>', '')
            
            h+=1
            v.append(e)
            e=''
        else:
            e+=s
        s=f.readline()
        if s=='':
            v.append(e)
    v.remove('')

    a=''
    err, error=0, []
    p=len(v)
    for i in range(p):
        v[i]=v[i].replace('\n','') #replacing the newlines by empty strings
        #remove the comment notation from one or more of the following four comment lines accordingly if you need to study a part of the sequence for all cases and not whole
        #v[i]=v[i][291:794]
        #start_index=v[i].find('IRCI')
        #end_index=v[i].find('AVSA')
        #v[i]=v[i][start_index:(end_index+4)]
        k=len(v[i])
        v[i]=v[i].upper()
        for j in range(k):
            if v[i][j] in ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']:
                a=a+v[i][j]
                #removing any such situations where any other wrong base apart from a, t, c, g is present due to publishing mistake
            else:
                err+=1
        error.append(err)
        err=0
        v[i]=a
        a=''
        w.append(len(v[i]))
    
    h=len(r)
    pr_record=[]
    length=[]
    w2=w
    w=set(w)
    w=list(w)
    dict1={}
    for i in range(len(w)):
        nm=str(w[i])
        os.mkdir(f'{out_path}/{nm}')
        dict1[f'{nm}']=[]
        file_n=f'{out_path}/{nm}/record.csv'
        with open(file_n, 'a+', newline='') as csvFile:
            writer=csv.writer(csvFile)
            r1=['Acc. ID', 'Details', 'rejected a. acids', 'pr']
            r1=list(r1)
            for lp in range(w[i]-(window_s-1)):
                r1.append(str(lp+1) + ' to ' + str(lp+window_s))
            writer.writerow(r1)
            r1=[]
            r11=['']*7
            r11=list(r11)
            writer.writerow(r11)
        csvFile.close()
    
    r1=0
    w1=[]
    bar['value']=5
    tk.Label(win, text='*********** Processing...(5%) ***********', background='white', foreground='blue').grid(row=25)
    for i in range(p):
        w1=len(v[i])
        a1=a2=a3=a4=a4=a5=a6=a7=a8=a9=a10=a11=a12=a13=a14=a15=a16=a17=a18=a19=a20=ua1=ua2=ua3=ua4=ua4=ua5=ua6=ua7=ua8=ua9=ua10=ua11=ua12=ua13=ua14=ua15=ua16=ua17=ua18=ua19=ua20=0
        for d in v[i]:
            if d=='A':
                a1+=1
            elif d=='R':
                a2+=1
            elif d=='N':
                a3+=1
            elif d=='D':
                a4+=1
            elif d=='C':
                a5+=1
            elif d=='Q':
                a6+=1
            elif d=='E':
                a7+=1
            elif d=='G':
                a8+=1
            elif d=='H':
                a9+=1
            elif d=='I':
                a10+=1
            elif d=='L':
                a11+=1
            elif d=='K':
                a12+=1
            elif d=='M':
                a13+=1
            elif d=='F':
                a14+=1
            elif d=='P':
                a15+=1
            elif d=='S':
                a16+=1
            elif d=='T':
                a17+=1
            elif d=='W':
                a18+=1
            elif d=='Y':
                a19+=1
            elif d=='V':
                a20+=1
            ua1+=a1
            ua2+=a2
            ua3+=a3
            ua4+=a4
            ua5+=a5
            ua6+=a6
            ua7+=a7
            ua8+=a8
            ua9+=a9
            ua10+=a10
            ua11+=a11
            ua12+=a12
            ua13+=a13
            ua14+=a14
            ua15+=a15
            ua16+=a16
            ua17+=a17
            ua18+=a18
            ua19+=a19
            ua20+=a20
            
        if w1!=0:
            ua1, ua2, ua3, ua4, ua5, ua6, ua7, ua8, ua9, ua10, ua11, ua12, ua13, ua14, ua15, ua16, ua17, ua18, ua19, ua20= ua1/w1, ua2/w1, ua3/w1, ua4/w1, ua5/w1, ua6/w1, ua7/w1, ua8/w1, ua9/w1, ua10/w1, ua11/w1, ua12/w1, ua13/w1, ua14/w1, ua15/w1, ua16/w1, ua17/w1, ua18/w1, ua19/w1, ua20/w1
            pr=math.sqrt(ua1*ua1 + ua2*ua2 + ua3*ua3 + ua4*ua4 + ua5*ua5 + ua6*ua6 + ua7*ua7 + ua8*ua8 + ua9*ua9 + ua10*ua10 + ua11*ua11 + ua12*ua12 + ua13*ua13 + ua14*ua14 + ua15*ua15 + ua16*ua16 + ua17*ua17 + ua18*ua18 + ua19*ua19 + ua20*ua20)
            pr_record.append(pr)
            r[i].append(error[i])
            r[i].append(pr)
            name=str(w1)
            file_n=f'{out_path}/{name}/record.csv'
            with open(file_n, 'a+', newline='') as csvFile:
                writer=csv.writer(csvFile)
                for lp in range(w1-(window_s-1)):
                    stg=v[i][lp:(lp+window_s)]
                    a1=a2=a3=a4=a4=a5=a6=a7=a8=a9=a10=a11=a12=a13=a14=a15=a16=a17=a18=a19=a20=ua1=ua2=ua3=ua4=ua4=ua5=ua6=ua7=ua8=ua9=ua10=ua11=ua12=ua13=ua14=ua15=ua16=ua17=ua18=ua19=ua20=0
                    for d in stg:
                        if d=='A':
                            a1+=1
                        elif d=='R':
                            a2+=1
                        elif d=='N':
                            a3+=1
                        elif d=='D':
                            a4+=1
                        elif d=='C':
                            a5+=1
                        elif d=='Q':
                            a6+=1
                        elif d=='E':
                            a7+=1
                        elif d=='G':
                            a8+=1
                        elif d=='H':
                            a9+=1
                        elif d=='I':
                            a10+=1
                        elif d=='L':
                            a11+=1
                        elif d=='K':
                            a12+=1
                        elif d=='M':
                            a13+=1
                        elif d=='F':
                            a14+=1
                        elif d=='P':
                            a15+=1
                        elif d=='S':
                            a16+=1
                        elif d=='T':
                            a17+=1
                        elif d=='W':
                            a18+=1
                        elif d=='Y':
                            a19+=1
                        elif d=='V':
                            a20+=1
                        ua1+=a1
                        ua2+=a2
                        ua3+=a3
                        ua4+=a4
                        ua5+=a5
                        ua6+=a6
                        ua7+=a7
                        ua8+=a8
                        ua9+=a9
                        ua10+=a10
                        ua11+=a11
                        ua12+=a12
                        ua13+=a13
                        ua14+=a14
                        ua15+=a15
                        ua16+=a16
                        ua17+=a17
                        ua18+=a18
                        ua19+=a19
                        ua20+=a20
                    ua1, ua2, ua3, ua4, ua5, ua6, ua7, ua8, ua9, ua10, ua11, ua12, ua13, ua14, ua15, ua16, ua17, ua18, ua19, ua20= ua1/w1, ua2/w1, ua3/w1, ua4/w1, ua5/w1, ua6/w1, ua7/w1, ua8/w1, ua9/w1, ua10/w1, ua11/w1, ua12/w1, ua13/w1, ua14/w1, ua15/w1, ua16/w1, ua17/w1, ua18/w1, ua19/w1, ua20/w1
                    pr=math.sqrt(ua1*ua1 + ua2*ua2 + ua3*ua3 + ua4*ua4 + ua5*ua5 + ua6*ua6 + ua7*ua7 + ua8*ua8 + ua9*ua9 + ua10*ua10 + ua11*ua11 + ua12*ua12 + ua13*ua13 + ua14*ua14 + ua15*ua15 + ua16*ua16 + ua17*ua17 + ua18*ua18 + ua19*ua19 + ua20*ua20)
                    r[i].append(pr)
                lis1=r[i][4:len(r[i])]
                lis1=list(lis1)
                dict1[f'{name}'].append(lis1)
                writer.writerow(r[i])
            csvFile.close()
    bar['value']=45
    tk.Label(win, text='Processing...(45%)', background='white', foreground='blue').grid(row=25)
    for j1 in w:
        name=str(j1)
        pv=[]
        file_n=f'{out_path}/{name}/record.csv'
        with open(file_n, 'a+', newline='') as csvFile1:
            writer=csv.writer(csvFile1)
            rw=['']*4
            for lp in range(len(dict1[f'{name}'][0])):
                rw1=[]
                for lines in dict1[f'{name}']:
                    rw1.append(lines[lp])
                rw1=set(rw1)
                rw1=list(rw1)
                rw.append(len(rw1))
                pv.append(len(rw1))
            writer.writerow(rw)
            pv2=[]
            for i in range(len(pv)):
                #pv1=pv[i:(i+window_s1)]
                #mean=sum(pv1)/window_s1
                pv2.append(pv[i])
        csvFile1.close()
        file_n=f'{out_path}/{name}/coord.csv'
        with open(file_n, 'a+', newline='') as csvFile:
            writer=csv.writer(csvFile)
            for i in range(len(pv2)):
                writer.writerow([(i+1), pv[i]])
        csvFile.close()
    bar['value']=65
    tk.Label(win, text='Processing...(65%)', background='white', foreground='blue').grid(row=25)
    for j1 in w:
        name=str(j1)
        file_n=f'{out_path}/{name}/record.csv'
        with open(file_n, 'a+', newline='') as csvFile1:
            writer=csv.writer(csvFile1)
            r2=[]
            for lp in range(j1):
                r2.append(' ')
            r2=list(r2)
            writer.writerow(r2)
            basep=w2.count(j1)
            sum1=0
            for i in range(len(w2)):
                if w2[i]==j1:
                    sum1+=pr_record[i]
            mean, variance, sd=0, 0, 0
            if basep!=0:
                mean=sum1/basep
            for i in range(len(w2)):
                if w2[i]==j1:
                    variance+=(pr_record[i]-mean)*(pr_record[i]-mean)
            if basep!=0:
                variance=variance/basep
            sd=math.sqrt(variance)
            writer.writerow(['  ', '    '])
            writer.writerow(['mean', mean])
            writer.writerow(['Std. dev.', sd])
        csvFile1.close()
    length=0
    bar['value']=80
    tk.Label(win, text='Processing...(80%)', background='white', foreground='blue').grid(row=25)
    for i in range(p):
        if '_' in r[i][0]:
            length=len(v[i])
            break
    bar['value']=90
    tk.Label(win, text='Processing...(90%)', background='white', foreground='blue').grid(row=25)
    print(length)
    return length
tk.Label(win, text='      .........................................................................................................................................                     ', background='yellow', foreground='yellow').grid(row=1)
tk.Label(win, text='             PEPTIDE VACCINE ANALYSIS TOOL - PVAT              ', background='yellow', foreground='black', font=('Verdana', '22')).grid(row=2)
tk.Label(win, text='_______________________________________________________________________________________________________________', background='yellow', foreground='yellow').grid(row=6)         
tk.Label(win, text='                 ', background='white', foreground='blue').grid(row=7)         
tk.Label(win, text='Enter input pathname for the fasta file: ', background='white', foreground='black', font=('Verdana', '14', 'bold')).grid(row=8)
e1=tk.Entry(win, width=40, font=('Courier', '12'), cursor='pencil')
e1.grid(row=9)
tk.Label(win, text='       ', background='white', foreground='black', font=('Verdana', '14', 'bold')).grid(row=10)
tk.Label(win, text='Enter input folder pathname containing surface accessibility (ASA) values:', background='white', foreground='brown', font=('Verdana', '14', 'bold')).grid(row=11)
e2=tk.Entry(win, width=40, font=('Courier', '12'), cursor='pencil')
e2.grid(row=12)
tk.Label(win, text='    ', background='white', foreground='black', font=('Verdana', '14', 'bold')).grid(row=16)
tk.Label(win, text='Enter output pathname: ', background='white', foreground='brown', font=('Verdana', '14', 'bold')).grid(row=17)
e3=tk.Entry(win, width=40, font=('Courier', '12'), cursor='pencil')
e3.grid(row=18)
tk.Label(win, text='      ', background='white', foreground='black').grid(row=24)
#img=tk.PhotoImage(file="protein.png")
#panel=tk.Label(win, image=img, background='grey').grid(row=25)
tk.Label(win, text='*Your progress will be shown here', background='white', foreground='red').grid(row=25)
bar.grid(row=26)
bar['value']=0
#tk.Label(win, text='      ', background='white', foreground='white').grid(row=26)
tk.Label(win, text='      ', background='white', foreground='black').grid(row=27)
tk.Button(win, width=12, text='Quit', cursor='pirate', command=win.destroy, foreground='black', relief='flat', font=('Verdana', '14'), highlightthickness=4, highlightbackground='#FF7676').place(x=179, y=410)
tk.Button(win, width=12, text='Submit', cursor='hand1', command=func_call, foreground='black', relief='flat', font=('Verdana', '14'), highlightthickness=4, highlightbackground='#66CC00').place(x=389.50, y=410)
tk.mainloop()

import tkinter as tk
import matplotlib.pyplot as plt
win=tk.Tk()
win.title('SPLIT AND PR (NCBI)')
win.resizable(0, 0)

win.configure(background='black')

main_path, out_path='', ''
window_s=12
window_s1=12
def retrieval():
    import csv
    import math
    import os
    
    main_path=e1.get()
    out_path=e2.get()
    main_path=main_path.replace('"', '')
    out_path=out_path.replace('"', '')
    f=open(main_path, 'r')
    
    s=f.readline()
    s1=''
    c, h, r, v, e, w=-1, 0, [], [], '', []
    while s!='':
        s=s.replace('\n', '')
        if s.isalpha():
            lk=list(set(list(s)))
            if sorted(lk)==['A', 'C', 'G', 'T']:
                e+=s
        else:
            s=s.replace('>', '')
            s=s.replace('\n', ' ')
            r.append([])
            c+=1
            r[c].append(s)
            v.append(e)
            e=''
        s=f.readline()
        if s=='':
            v.append(e)
    ec=v.count('')
    for i in range(ec):
        v.remove('')
    a=''
    err, error=0, []
    p=len(v)
    for i in range(p):
        v[i]=v[i].replace('\n','') #replacing the newlines by empty strings
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
                    x1=y1=ux1=uy1=0
                    for d1 in stg:
                        if d1=='A':
                            x1-=1
                        elif d1=='G':
                            x1+=1
                        elif d1=='C':
                            y1+=1
                        elif d1=='T':
                            y1-=1
                        ux1+=x1
                        uy1+=y1
                    ux1, uy1=ux1/window_s, uy1/window_s
                    gr1=math.sqrt(ux1*ux1 + uy1*uy1)
                    r[i].append(gr1)
                lis1=r[i][4:len(r[i])]
                lis1=list(lis1)
                dict1[f'{name}'].append(lis1)
                writer.writerow(r[i])
            csvFile.close()
                
    w2.remove(0)
    w.remove(0)
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
            for i in range(len(pv)-(window_s1-1)):
                pv1=pv[i:(i+window_s1)]
                mean=sum(pv1)/window_s1
                pv2.append(mean)
        csvFile1.close()
        file_n=f'{out_path}/{name}/coord.csv'
        with open(file_n, 'a+', newline='') as csvFile:
            writer=csv.writer(csvFile)
            for i in range(len(pv2)):
                writer.writerow([(i+1), pv2[i]])
        csvFile.close()
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
tk.Label(win, text='                                                                                                                      ', background='black', foreground='yellow').grid(row=1)
tk.Label(win, text='   WELCOME TO PR CALCULATOR !!   ', background='black', foreground='yellow', font=('Jokerman', '22')).grid(row=2)
tk.Label(win, text='...........................................................................................................', background='black', foreground='yellow').grid(row=5)
tk.Label(win, text='', background='black', foreground='white').grid(row=6)         
tk.Label(win, text='Enter input path: ', background='black', foreground='cyan', font=('Rockwell', '14', 'bold')).grid(row=7)
tk.Label(win, text='Enter output path: ', background='black', foreground='cyan', font=('Rockwell', '14', 'bold')).grid(row=9)
win.iconbitmap("Iconsmind-Outline-DNA.ico")
e1=tk.Entry(win, width=40, font=('Courier', '12'))
e2=tk.Entry(win, width=40, font=('Courier', '12'))
    
e1.grid(row=8)
e2.grid(row=10)

tk.Label(win, text='', background='black', foreground='white').grid(row=16)
img=tk.PhotoImage(file="protein.png")
panel=tk.Label(win, image=img, background='black').grid(row=17)
tk.Label(win, text='', background='black', foreground='white').grid(row=18)
tk.Label(win, text='', background='black', foreground='white').grid(row=19)
tk.Label(win, text='', background='black', foreground='white').grid(row=20)
tk.Label(win, text='', background='black', foreground='white').grid(row=21)
tk.Button(win, width=12, text='Quit', command=win.destroy, foreground='brown', relief='ridge', font=('Segoe Script', '14'), highlightthickness=4, highlightbackground='red', highlightcolor='white').place(x=70.75, y=595)
tk.Button(win, width=12, text='Submit', command=retrieval, foreground='dark green', relief='ridge', font=('Segoe Script', '14'), highlightthickness=4, highlightbackground='green', highlightcolor='white').place(x=304.25, y=595)

tk.mainloop()

import tkinter as tk
import matplotlib.pyplot as plt
win=tk.Tk()
win.title('SPLIT AND PR')
win.resizable(0, 0)

win.configure(background='black')

main_path, out_path='', ''
window_s=''
window_s1=''
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

    amino=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    a=''
    err, error=0, []
    p=len(v)
    for i in range(p):
        v[i]=v[i].replace('\n','') 
        k=len(v[i])
        v[i]=v[i].upper()
        for j in range(k):
            if v[i][j] in amino:
                a=a+v[i][j]
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
            r1=['Details', 'rejected a. acids', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
            r1=list(r1)
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
        r[i].append(error[i])
        for j in amino:
            lp=v[i].count(j)
            r[i].append(lp)
        name=str(w1)
        file_n=f'{out_path}/{name}/record.csv'
        with open(file_n, 'a+', newline='') as csvFile:
            writer=csv.writer(csvFile)
            writer.writerow(r[i])
        csvFile.close()
                
tk.Label(win, text='                                                                                                                      ', background='black', foreground='yellow').grid(row=1)
tk.Label(win, text='   AMINO ACID FREQUENCY CALCULATOR !!   ', background='black', foreground='yellow', font=('Jokerman', '22')).grid(row=2)
tk.Label(win, text='.............................................................................................................................................................................................', background='black', foreground='yellow').grid(row=5)
tk.Label(win, text='', background='black', foreground='white').grid(row=6)         
tk.Label(win, text='Enter input path: ', background='black', foreground='cyan', font=('Rockwell', '14', 'bold')).grid(row=7)
tk.Label(win, text='Enter output path: ', background='black', foreground='cyan', font=('Rockwell', '14', 'bold')).grid(row=9)
win.iconbitmap("Iconsmind-Outline-DNA.ico")
e1=tk.Entry(win, width=45, font=('Courier', '12'))
e2=tk.Entry(win, width=45, font=('Courier', '12'))
    
e1.grid(row=8)
e2.grid(row=10)
tk.Label(win, text='', background='black', foreground='white').grid(row=11)

img=tk.PhotoImage(file="protein.png")
panel=tk.Label(win, image=img, background='black').grid(row=17)
tk.Label(win, text='', background='black', foreground='white').grid(row=18)
tk.Label(win, text='', background='black', foreground='white').grid(row=19)
tk.Label(win, text='', background='black', foreground='white').grid(row=20)
tk.Label(win, text='', background='black', foreground='white').grid(row=21)

tk.Button(win, width=12, text='Quit', command=win.destroy, foreground='brown', relief='ridge', font=('Segoe Script', '14'), highlightthickness=4, highlightbackground='red', highlightcolor='white').place(x=47.75, y=595)
tk.Button(win, width=12, text='Submit', command=retrieval, foreground='dark green', relief='ridge', font=('Segoe Script', '14'), highlightthickness=4, highlightbackground='green', highlightcolor='white').place(x=451.25, y=595)

tk.mainloop()

import tkinter as tk
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
win=tk.Tk()
win.title('SPLIT AND GR')
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
    window_s=e3.get()
    window_s=int(window_s)
    main_path=main_path.replace('"', '')
    out_path=out_path.replace('"', '')
    f=open(main_path, 'r')
    
    s=f.readline()
    s1=''
    c, h, r, v, e, w=-1, 0, [], [], '', []
    while s!='':
        if s[0]=='>':
            r.append(s.split('|'))
            r[h][0]=r[h][0].replace('>', '')
            hw=r[h][0].split(':')
            r[h][0]=hw[0]
            r[h].insert(1, hw[1])
            
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
        k=len(v[i])
        v[i]=v[i].upper()
        for j in range(k):
            if v[i][j] in ['A', 'T', 'C', 'G']:
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
    gr_record=[]
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
            r1=['Acc. ID', 'CDS region', 'Host', 'Country', 'Date of coll.', 'rejected bases', 'gr', 'slope', 'y-intercept']
            for lp in range(w[i]-(window_s-1)):
                r1.append(str(lp+1) + ' to ' + str(lp+window_s))
            r1=list(r1)
            writer.writerow(r1)
            r1=[]
            r11=['']*(w[i] - (window_s-8))
            r11=list(r11)
            writer.writerow(r11)
        csvFile.close()
    
    r1=0
    w1=[]
    cdx, cdy, slavg, yincavg, flag=[], [], 0, 0, 0
    for j1 in w:
        name=str(j1)
        file_n=f'{out_path}/{name}/quadratic_coeff.csv'
        with open(file_n, 'a+', newline='') as csvFile:
            writer=csv.writer(csvFile)
            writer.writerow(['ACC ID', 'Coeff. of x^2', 'Coeff. of x', 'Constant', '', 'ux', 'uy', 'gr', 'angle'])
            writer.writerow(['', '', '', '', '', '', '', '', ''])
        csvFile.close()
    for i in range(p):
        cdx.append([])
        cdy.append([])
        w1=len(v[i])
        x=y=ux=uy=0
        coef1, coef2, coef3, coef4, coef5, ds1, ds2, ds3 = 0, 0, 0, 0, 0, 0, 0, 0
        cordx, cordy=[], []
        for d in v[i]:
            if d=='A':
                x-=1
                coef1 += x*x*x*x
                coef2 += x*x*x
                coef3 += x*x
                coef4 += x
                coef5 += 1
                ds1 += x*x*y
                ds2 += x*y
                ds3 += y
            elif d=='G':
                x+=1
                coef1 += x*x*x*x
                coef2 += x*x*x
                coef3 += x*x
                coef4 += x
                coef5 += 1
                ds1 += x*x*y
                ds2 += x*y
                ds3 += y
            elif d=='C':
                y+=1
                coef1 += x*x*x*x
                coef2 += x*x*x
                coef3 += x*x
                coef4 += x
                coef5 += 1
                ds1 += x*x*y
                ds2 += x*y
                ds3 += y
            elif d=='T':
                y-=1
                coef1 += x*x*x*x
                coef2 += x*x*x
                coef3 += x*x
                coef4 += x
                coef5 += 1
                ds1 += x*x*y
                ds2 += x*y
                ds3 += y
            cordx.append(x)
            cordy.append(y)
            cdx[i].append(x)
            cdy[i].append(y)
            ux+=x
            uy+=y
        if w1!=0:
            sl = ((coef5*ds2) - (ds3*coef4))/((coef5*coef3) - (coef4*coef4))
            yinc = (ds3 - (sl*coef4))/coef5
            sumx, sumy=ux, uy
            ux, uy=ux/w1, uy/w1
            gr=math.sqrt(ux*ux + uy*uy)
            gr_record.append(gr)
            r[i].append(error[i])
            r[i].append(gr)
            r[i].append(sl)
            r[i].append(yinc)
            r[i][3]=r[i][3].replace('|', '')
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
            acc1=str(r[i][0])
            simul1=np.array([[coef1, coef2, coef3], [coef2, coef3, coef4], [coef3, coef4, coef5]])
            simul2=np.array([ds1, ds2, ds3])
            simul3=np.linalg.solve(simul1, simul2)
            file_n=f'{out_path}/{name}/coordinates_{acc1}.csv'
            with open(file_n, 'a+', newline='') as csvFile:
                writer=csv.writer(csvFile)
                writer.writerow(['X', 'Y'])
                for lp in range(len(cordx)):
                    writer.writerow([cordx[lp], cordy[lp]])
                writer.writerow(['', ''])
                writer.writerow([sumx, sumy])
                writer.writerow(['', ''])
                writer.writerow(['A', v[i].count('A')])
                writer.writerow(['T', v[i].count('T')])
                writer.writerow(['G', v[i].count('C')])
                writer.writerow(['C', v[i].count('G')])
                writer.writerow(['Unknown', error[i]])
                writer.writerow(['', ''])
                writer.writerow(['GR', gr])
            csvFile.close()
            file_n=f'{out_path}/{name}/quadratic_coeff.csv'
            with open(file_n, 'a+', newline='') as csvFile:
                writer=csv.writer(csvFile)
                if ux==0:
                    if uy==0:
                        atan1=0
                    elif uy>0:
                        atan1=90
                    elif uy<0:
                        atan1=270
                elif uy==0:
                    if ux>=0:
                        atan1=0
                    elif ux<0:
                        atan1=180
                elif ux>0 and uy>0:
                    atan1=math.atan(uy/ux)
                    atan1=(180/3.1415926535897)*atan1
                elif ux>0 and uy<0:
                    atan1=math.atan(abs(uy/ux))
                    atan1=(180/3.1415926535897)*atan1
                    atan1=360-atan1
                elif ux<0 and uy>0:
                    atan1=math.atan(abs(uy/ux))
                    atan1=(180/3.1415926535897)*atan1
                    atan1=180-atan1
                else:
                    atan1=math.atan(abs(uy/ux))
                    atan1=(180/3.1415926535897)*atan1
                    atan1=180+atan1
                writer.writerow([acc1, simul3[0], simul3[1], simul3[2], '', ux, uy, gr, atan1])
            csvFile.close()

    for j1 in w:
        name=str(j1)
        file_n=f'{out_path}/{name}/record.csv'
        slope1, yinc1, time1=[], [], []
        for lines in dict1[f'{name}']:
            s=str(lines[0])
            s=s.split('/')
            if s[0].isnumeric():
                time1.append(int(s[0]))
                slope1.append(lines[3])
                yinc1.append(lines[4])
        fig=plt.figure()
        ax=fig.add_subplot(111, projection='3d')
        ax.scatter(slope1, time1, yinc1, c='r', marker='o', s=0.8)
        ax.view_init(elev=20, azim=45)
        ax.set_xlabel('Slope')
        ax.set_zlabel('Y-intercept')
        ax.set_ylabel('Time in Years')
        #plt.show()
        #plt.plot(slope1, yinc1, 'ro', markersize=1)
        #plt.axis([0.5, 0.7, -37, -15])
        fig.savefig(f'{out_path}/{name}/similarity_plot.png')
        plt.close()
        with open(file_n, 'a+', newline='') as csvFile1:
            writer=csv.writer(csvFile1)
            rw=['']*8
            for lp in range(5, len(dict1[f'{name}'][0])):
                rw1=[]
                for lines in dict1[f'{name}']:
                    rw1.append(lines[lp])
                rw1=set(rw1)
                rw1=list(rw1)
                rw.append(len(rw1))
            writer.writerow(rw)
        csvFile1.close()
    for j1 in w:
        name=str(j1)
        plt.figure()
        for i in range(p):
            if len(v[i])==j1:
                plt.plot(list(cdx[i]), list(cdy[i]), linewidth=0.5)
        plt.savefig(f'{out_path}/{name}/plot_{j1}.png')
        plt.close()
    for j1 in w:
        name=str(j1)
        file_n=f'{out_path}/{name}/record.csv'
        with open(file_n, 'a+', newline='') as csvFile1:
            writer=csv.writer(csvFile1)
            r2=[]
            for lp in range(j1-(window_s-7)):
                r2.append(' ')
            r2=list(r2)
            writer.writerow(r2)
            basep=w2.count(j1)
            sum1=0
            for i in range(len(w2)):
                if w2[i]==j1:
                    sum1+=gr_record[i]
            mean, variance, sd=0, 0, 0
            if basep!=0:
                mean=sum1/basep
            for i in range(len(w2)):
                if w2[i]==j1:
                    variance+=(gr_record[i]-mean)*(gr_record[i]-mean)
            if basep!=0:
                variance=variance/basep
            sd=math.sqrt(variance)
            writer.writerow(['  ', '    '])
            writer.writerow(['mean', mean])
            writer.writerow(['Std. dev.', sd])
        csvFile1.close()
tk.Label(win, text='                                                                                                                      ', background='black', foreground='yellow').grid(row=1)
tk.Label(win, text=' WELCOME TO GR CALCULATOR !! ', background='black', foreground='yellow', font=('Jokerman', '22')).grid(row=2)
tk.Label(win, text='..................................................................................................................................', background='black', foreground='yellow').grid(row=5)
tk.Label(win, text='', background='black', foreground='white').grid(row=6)         
tk.Label(win, text='Enter input path: ', background='black', foreground='cyan', font=('Rockwell', '14', 'bold')).grid(row=7)
tk.Label(win, text='Enter output path: ', background='black', foreground='cyan', font=('Rockwell', '14', 'bold')).grid(row=9)
win.iconbitmap("Iconsmind-Outline-DNA.ico")
e1=tk.Entry(win, width=53, font=('Courier', '12'))
e2=tk.Entry(win, width=53, font=('Courier', '12'))
    
e1.grid(row=8)
e2.grid(row=10)
tk.Label(win, text='', background='black', foreground='white').grid(row=11)
tk.Label(win, text='Enter window length for window_wise gr: ', background='black', foreground='pink', font=('Rockwell', '14', 'bold')).grid(row=12)
e3=tk.Spinbox(win, from_=2, to = 100, width=10, font=('Courier'))

e3.grid(row=13)
tk.Label(win, text='', background='black', foreground='white').grid(row=16)
img=tk.PhotoImage(file="giphy.gif")
panel=tk.Label(win, image=img, background='black').grid(row=17)
tk.Label(win, text='', background='black', foreground='white').grid(row=18)
tk.Label(win, text='', background='black', foreground='white').grid(row=19)
tk.Label(win, text='', background='black', foreground='white').grid(row=20)
tk.Button(win, width=12, text='Quit', command=win.destroy, foreground='brown', relief='ridge', font=('Segoe Script', '14'), highlightthickness=4, highlightbackground='red', highlightcolor='white').place(x=38.75, y=590)
tk.Button(win, width=12, text='Submit', command=retrieval, foreground='dark green', relief='ridge', font=('Segoe Script', '14'), highlightthickness=4, highlightbackground='green', highlightcolor='white').place(x=375.25, y=590)

tk.mainloop()

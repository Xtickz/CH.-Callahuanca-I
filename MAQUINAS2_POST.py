# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:42:06 2021

@author: USER
"""

import tkinter as tk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
#DATOS DE LA CENTRAL
S=44*10**6
V=8*10**3
P=30.80*10**6
#por unidad
s=1
v=1
i=1
p=P/S
xd=1.450
xq=0.850
#FUNCIONES
def vect(M,the):
    x_d=M*np.cos(the)
    y_d=M*np.sin(the)
    return x_d,y_d
def Eq(fp,V=v,I=i,Xq=xq):
    fpa=np.arccos(fp)
    return np.sqrt((I*Xq*np.cos(fpa))**2 +(V+I*Xq*np.sin(fpa))**2)
def dt(fp,V=v,I=i,Xq=xq):
    fpa=np.arccos(fp)
    return np.arctan((I*Xq*np.cos(fpa))/(V+I*Xq*np.sin(fpa)))
def EF(fp,eq,delta,I=i):
    fpa=np.arccos(fp)
    return eq+(xd-xq)*I*np.sin(fpa+delta)
def P(Ef,g,V=v,Xd=xd,Xq=xq):
    return (Ef*V*np.sin(g))/(Xd),((V**2)/2)*(1/Xq -1/Xd)*np.sin(2*g)
def Q(Ef,g,V=v,Xd=xd,Xq=xq):
    return (Ef*V*np.cos(g))/(Xd),((V**2)/2)*(1/Xq -1/Xd)*np.cos(2*g),-((V**2)/2)*(1/Xq +1/Xd)
def preg_1():
    #miFrame_2.delete()
    g=np.linspace(-np.pi,np.pi)
    fig, axes = plt.subplots(ncols=3, nrows=2,figsize=(12.5,7.5),sharey=True)
    for i in [0,1]:
        for j in [0,1,2]:
            axes[i,j].spines['right'].set_color('none') 
            axes[i,j].spines['top'].set_color('none')         
            axes[i,j].xaxis.set_ticks_position('bottom')   
            axes[i,j].yaxis.set_ticks_position('left')          
            axes[i,j].spines['bottom'].set_position(('data', 0))   
            axes[i,j].spines['left'].set_position(('data', 0))
    #1.1
    Ef1=float(var1.get())
    P1p,P1pp=P(Ef1,g)
    Q1p,Q1pp,Q1ppp=Q(Ef1,g)
    axes[0,0].plot(g,P1p+P1pp,'r')
    axes[1,0].plot(g,Q1p+Q1pp+Q1ppp*np.ones(50),'r')
    axes[0,1].plot(g,P1p,'black')
    axes[0,1].plot(g,P1pp,'black')
    axes[0,1].plot(g,P1p+P1pp,'r')
    axes[1,1].plot(g,Q1p,'black')
    axes[1,1].plot(g,Q1pp,'black')
    axes[1,1].plot(g,Q1ppp*np.ones(50),'black')
    axes[1,1].plot(g,Q1p+Q1pp+Q1ppp*np.ones(50),'r')
    #1.2
    Ef2=float(var2.get())
    P2p,P2pp=P(Ef2,g)
    Q2p,Q2pp,Q2ppp=Q(Ef2,g)
    axes[0,0].plot(g,P2p+P2pp,'b')
    axes[1,0].plot(g,Q2p+Q2pp+Q2ppp*np.ones(50),'b')
    axes[0,2].plot(g,P2p,'black')
    axes[0,2].plot(g,P2pp,'black')
    axes[0,2].plot(g,P2p+P2pp,'b')
    axes[1,2].plot(g,Q2p,'black')
    axes[1,2].plot(g,Q2pp,'black')
    axes[1,2].plot(g,Q2ppp*np.ones(50),'black')
    axes[1,2].plot(g,Q2p+Q2pp+Q2ppp*np.ones(50),'b')
    axes[0,0].set_title("P1 y P2")
    axes[1,0].set_title("Q1 Y Q2")
    axes[0,1].set_title("Descomposicion de P1")
    axes[0,2].set_title("Descomposicion de P2")
    axes[1,1].set_title("Descomposicion de Q1")
    axes[1,2].set_title("Descomposicion de Q2")
    canvas=FigureCanvasTkAgg(fig,master=miFrame_2)
    canvas.draw()
    canvas.get_tk_widget().place(x=0,y=0)
    #toolbar=NavigationToolbar2Tk(canvas,root)
    #toolbar.update()
    #canvas.get_tk_widget().place(x=0,y=0)
    texto.delete(1.0,"end")
    CONT_1=1.0
    for EF_H in [Ef1,Ef2]:
        H=np.array([np.pi/6,np.pi/4,np.pi/3,np.pi/2,(2*np.pi)/3])
        PH_p,PH_pp=P(EF_H,H)
        QH_p,QH_pp,QH_ppp=Q(EF_H,H)
        texto.insert(CONT_1," Para Ef= {}: g={}\n ".format(EF_H,H))
        texto.insert(CONT_1+1,"P:{} \n".format(PH_p+PH_pp))
        texto.insert(CONT_1+2,"Q:{} \n".format(QH_p+QH_pp+QH_ppp))
        CONT_1=CONT_1+3.0
def preg_2():
    fig1, axes1 = plt.subplots(ncols=3, nrows=2,figsize=(12.5,7.5),sharey=True,sharex=True)
#ax1 = fig.gca()
    for i in [0,1]:
        for j in [0,1,2]:
            axes1[i,j].spines['right'].set_color('none') 
            axes1[i,j].spines['top'].set_color('none')         
            axes1[i,j].xaxis.set_ticks_position('bottom')   
            axes1[i,j].yaxis.set_ticks_position('left')          
            axes1[i,j].spines['bottom'].set_position(('data', 0))   
            axes1[i,j].spines['left'].set_position(('data', 0))
    j=0
    for k in [0,0.6,0.8,0,0.6,0.8]:
        fp_k=np.arccos(k)
        eq_k=Eq(k)
        delta_k=dt(k)
        ef_k=EF(k,eq_k,delta_k)
        xd_k=[0,0,0]
        yd_k=[0,0,0]
        if j<3:
            xd_k[0],yd_k[0]=vect(ef_k,np.pi/2 +delta_k)
            xd_k[1],yd_k[1]=vect(1,np.pi/2)
            xd_k[2],yd_k[2]=vect(1,np.pi/2 - fp_k)
            axes1[0,j].quiver([0,0],[0,0],[xd_k[1],xd_k[2]],[yd_k[1],yd_k[2]],scale_units='xy',angles='xy',scale=1)
            axes1[0,j].quiver(0,0,xd_k[0],yd_k[0],scale_units='xy',angles='xy',scale=1,color='red')
        else:
            xd_k[0],yd_k[0]=vect(ef_k,np.pi/2 -delta_k)
            xd_k[1],yd_k[1]=vect(1,np.pi/2)
            xd_k[2],yd_k[2]=vect(1,np.pi/2 + fp_k)
            axes1[1,j-3].quiver([0,0],[0,0],[xd_k[1],xd_k[2]],[yd_k[1],yd_k[2]],scale_units='xy',angles='xy',scale=1)
            axes1[1,j-3].quiver(0,0,xd_k[0],yd_k[0],scale_units='xy',angles='xy',scale=1,color='red')
        plt.ylim(-0.5,2.5)
        plt.xlim(-1.5,1.5)
        j=j+1
    axes1[0,0].set_title("f.d.p=0 inductivo")
    axes1[1,0].set_title("f.d.p=0 capacitivo")
    axes1[0,1].set_title("f.d.p=0.6 inductivo")
    axes1[0,2].set_title("f.d.p=0.8 inductivo")
    axes1[1,1].set_title("f.d.p=0.6 capacitivo")
    axes1[1,2].set_title("f.d.p=0.8 capacitivo")
    canvas=FigureCanvasTkAgg(fig1,master=miFrame_2)
    canvas.draw()
    canvas.get_tk_widget().place(x=0,y=0)
    fp=np.array([0,0.6,0.8,1,0.8,0.6,0])
    eq=Eq(fp)
    delta=dt(fp)
    ef=EF(fp,eq,delta)
    texto.delete(1.0,"end")
    texto.insert(1.0,"g=[0,0.6,0.8,1] inductivo y g=[0.8,0.6,0] capacitivo\n")
    texto.insert(2.0,"If: {}\n".format(ef))
    texto.insert(4.0,"delta (g): {}\n".format(delta*180/np.pi))
#
def preg_3():
    i=1.0
    fp=np.array([0,0.6,0.8,1,0.8,0.6,0])
    texto.delete(1.0,"end")
    for k in [1.25,0.75,0.5,0.25]:
        i_k=1*k
        eq_k=Eq(fp,I=i_k)
        delta_k=dt(fp,I=i_k)
        ef_k=EF(fp,eq_k,delta_k,I=i_k)
        texto.insert(i,"para k={} : fdp={}\n".format(k,fp))
        texto.insert(i+1,"Ef:{}\n".format(ef_k*V))
        texto.insert(i+3,"delta:{}\n".format(delta_k*180/np.pi))
        i=i+5
#
root=tk.Tk()
root.title("CALLAHUANCA I")
root.geometry("1300x675")
root.config(bg="midnight blue")

var1= tk.StringVar()
var2= tk.StringVar()
consola= tk.StringVar()
# Frame de los datos y las pregs
miFrame=tk.Frame(root)
miFrame.place(x=10,y=25)
miFrame.config(width=400,height=600,bg="gray")
## Label del titulo
Label_1=tk.Label(miFrame,text="C.H CALLAHUANCA I")
Label_1.config(justify="center",bg="firebrick",width=53,height=2,fg="white")
Label_1.place(x=10,y=25)
## Label de los datos
Label_2=tk.Label(miFrame,text="Sn= 44MVA  V=8KV  P=30.80MW  Xd= 1.450 pu  Xq=0.850 pu")
Label_2.config(justify="center",bg="firebrick",width=53,height=2,fg="white")
Label_2.place(x=10,y=65)
## Label de la primera preg
Label_3=tk.Label(miFrame,text="PARTE 1: Graficar P-g y Q-g con cada uno de sus componentes")
Label_3.config(justify="center",bg="tan",width=53,height=2,fg="white")
Label_3.place(x=10,y=125)

Label_4=tk.Label(miFrame,text=" V = 1.0 pu y Ef = 1.40 pu\n\n V = 1.0 pu y Ef = 0.85 pu")
Label_4.config(justify="center",bg="firebrick",width=20,height=5,fg="white")
Label_4.place(x=240,y=175)
## Creando las entradas y los botones para la pregunta 1
Label_5=tk.Label(miFrame,text=" Ef1\n\n\n Ef2")
Label_5.config(justify="center",bg="firebrick",width=3,height=5,fg="white")
Label_5.place(x=10,y=175)

entra_1=tk.Entry(miFrame)
entra_1.config(justify="left",bd=6,textvariable=var1)
entra_1.place(x=40,y=175)

entra_2=tk.Entry(miFrame)
entra_2.config(justify="left",bd=6,textvariable=var2)
entra_2.place(x=40,y=225)

boton_1=tk.Button(miFrame,text="enviar",command=preg_1)
boton_1.place(x=180,y=275)
#Creando el espacio para la grafica
miFrame_2=tk.Frame(root)
miFrame_2.place(x=420,y=25)
miFrame_2.config(width=900,height=500,bg="gray")
#Salida de consola
texto=tk.Text(root)
texto.config(width=110,height=8,padx=10,pady=10,state="normal")
texto.place(x=420,y=550)
##label y boton para la preg 2.1
Label_6=tk.Label(miFrame,text="PARTE 2.1: Calcular la corriente de excitacion pu y los angulos g")
Label_6.config(justify="center",bg="tan",width=53,height=2,fg="white")
Label_6.place(x=10,y=310)

Label_7=tk.Label(miFrame,text="Para g=[1,0.8,0.6,0] inductivo y g=[0.8,0.6,0] capacitivo")
Label_7.config(justify="center",bg="firebrick",width=53,height=2,fg="white")
Label_7.place(x=10,y=355)

boton_2=tk.Button(miFrame,text="enviar",command=preg_2)
boton_2.place(x=180,y=400)
##label y boton para la preg 2.2
Label_8=tk.Label(miFrame,text='''PARTE 2.2: Para los factores de potencia dados en 2.1 y si la corriente \n es I = K*In Calcular la excitación Ef y el ángulo delta''')
Label_8.config(justify="center",bg="tan",width=53,height=3,fg="white")
Label_8.place(x=10,y=435)

Label_9=tk.Label(miFrame,text='''para K = 1,25; 0,75; 0,50; 0,25''')
Label_9.config(justify="center",bg="firebrick",width=53,height=2,fg="white")
Label_9.place(x=10,y=495)

boton_3=tk.Button(miFrame,text="enviar",command=preg_3)
boton_3.place(x=180,y=545)
# firma
Label_10=tk.Label(root,text='''By: Gianfranco Yucra''')
Label_10.config(justify="center",bg="midnight blue",width=23,height=2,fg="yellow")
Label_10.config(font=("Curier",16))
Label_10.place(x=40,y=630)
root.mainloop()
import numpy as np
from scipy.integrate import odeint
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import time
import tkinter as tk

#dictionary 
dict = {
        
    "taxa": 0,
    "populacao": 1

}

# function that returns dz/dt
def model(z,t,u):
    x = z
    dxdt = u[dict.get("taxa")]*x*(1 - x/u[dict.get("populacao")])
    return dxdt

def tscrypt(t, u):
    if t > 50:
        u[dict.get("populacao")] = 1
        u[dict.get("taxa")] = 0.002
        
    if t > 100:
        u[dict.get("populacao")] = 200
        u[dict.get("taxa")] = 0.2
        
    if t > 200:
        u[dict.get("populacao")] = 500
        
    if t > 300:
        u[dict.get("populacao")] = 2000
        u[dict.get("taxa")] = 0.02
    
    if t > 400:
        u[dict.get("populacao")] = 200
        u[dict.get("taxa")] = 0.002           
        
    if t > 500:
        u[dict.get("populacao")] = 2000
        u[dict.get("taxa")] = 0.05
    return u

def loop():

  global i, n, z0, u, x, t, up, figure, root, plot, prange;

  while ((i in range(1,n)) & enable == 1):
  # span for next time step
    tspan = [t[i-1],t[i]]
    # solve for next step
    z = odeint(model,z0,tspan,args= (u,))
    # store solution for plotting
    x[i] = z[1]
    up[i] = u[dict.get("populacao")]
    # next initial condition
    z0 = z[1]
    #u = tscrypt(t[i], u)
    i += 1;
        
    if i%10 == 0:
      break;

  i -= 1;
   
  in_range = (i-prange)
  if in_range < 0:
      in_range = 0;

  max_of_y = np.amax(x[in_range:i]);
  if max_of_y < 5000:
    max_of_y = 5000;

  plt.clf()
  figure = Figure(figsize=(5, 4), dpi=100)
  plot = figure.add_subplot(1, 1, 1)
  plt.ylim(-(max_of_y*.05), (max_of_y*1.05))
  plt.plot(t[i-prange:i],up[i-prange:i],'g:',label='u(t)')
  plt.plot(t[i-prange:i],x[i-prange:i],'b-',label='x(t)')
  plt.ylabel('values')
  plt.xlabel('time')
  plt.legend(loc='best')
  plt.show()
  #time.sleep(.1)
  #print('\33[2J')
  i += 1;

  root.after(50, loop)



#==============================================================================================================================================================

class Page(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
    def show(self):
        self.lift()

class Page1(Page):
    def __init__(self, *args, **kwargs):
      Page.__init__(self, *args, **kwargs)
      self.label = tk.Label(self, text="This is page 1")
      self.label.pack(side="top", fill="both", expand=True)

      self.entry1 = tk.Entry(self)
      self.entry1.pack(side="top")

      btn1 = tk.Button(self, text="taxa", command=self.muda_par1)
      btn1.pack(side="top")

      self.entry2 = tk.Entry(self)
      self.entry2.pack(side="top")

      btn2 = tk.Button(self, text="populacao", command=self.muda_par2)
      btn2.pack(side="top")

      self.btnenable = tk.Button(self, text="enable", command=self.muda_enable)
      self.btnenable.pack(side="top")
      self.btnenable.bind("<Return>", self.muda_enable)


      global figure

      #self.canvas = FigureCanvasTkAgg(figure, self)
      #self.canvas.get_tk_widget().pack(side="top")

    def muda_text(self):
      newstr = "mudou";
      self.label["text"] = newstr;

    def update_canvas(self):
      self.canvas = FigureCanvasTkAgg(figure, self)

    def muda_par1(self):
      global u 

      u[dict.get("taxa")] = float(self.entry1.get())

    def muda_par2(self):
      global u 

      u[dict.get("populacao")] = float(self.entry2.get())

    def muda_enable(self, event=None):
      global enable

      if enable == 1:
        enable = 0;
      else:
        enable = 1;



class Page2(Page):
   def __init__(self, *args, **kwargs):
       Page.__init__(self, *args, **kwargs)
       label = tk.Label(self, text="This is page 2")
       label.pack(side="top", fill="both", expand=True)

class Page3(Page):
   def __init__(self, *args, **kwargs):
       Page.__init__(self, *args, **kwargs)
       label = tk.Label(self, text="This is page 3")
       label.pack(side="top", fill="both", expand=True)

class MainView(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        self.master.title("Mitochondria Graph Plotter")
        p1 = Page1(self)
        p2 = Page2(self)
        p3 = Page3(self)

        buttonframe = tk.Frame(self)
        container = tk.Frame(self)
        buttonframe.pack(side="top", fill="x", expand=False)
        container.pack(side="top", fill="both", expand=True)

        p1.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p2.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p3.place(in_=container, x=0, y=0, relwidth=1, relheight=1)

        b1 = tk.Button(buttonframe, text="Page 1", command=p1.lift)
        b2 = tk.Button(buttonframe, text="Page 2", command=p2.lift)
        b3 = tk.Button(buttonframe, text="Page 3", command=p3.lift)

        b1.pack(side="left")
        b2.pack(side="left")
        b3.pack(side="left")

        p1.show()



if __name__ == "__main__":

  figure = Figure(figsize=(5, 4), dpi=100)
  plot = figure.add_subplot(1, 1, 1)

  # initial condition
  enable = 1;
  z0 = 100

  u = {};
  u[dict.get("taxa")] = 0.2
  u[dict.get("populacao")] = 500

  # number of time points
  n = 30000
  tf = 7000
  factor = int(n/tf)
  # time points
  t = np.linspace(0,tf,n)

  # store solution
  x = np.empty_like(t)
  up = np.empty_like(t)
  # record initial conditions
  x[0] = z0

  prange = 400

  i = 1

  # solve ODE
  #loop(i, n, z0, u, x, t, up);








  root = tk.Tk()
  main = MainView(root)
  main.pack(side="top", fill="both", expand=True)
  root.wm_geometry("700x500+300+100")

  root.after(50, loop)
  root.mainloop()






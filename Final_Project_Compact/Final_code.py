import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D    
import pylab
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
from scipy import *
import os
from PIL import Image,ImageTk
#Change working directory
os.chdir("C:/PH2150/Final_Project_Compact")
# Fonts used in GUI defined here
titlefont=("Times New Roman",15,"bold")
contentfont=("Times New Roman",13)
INSTRUCTION_FONT = ("Times New Roman",12, "bold")


class MainMenu(Tk):
    """Making all the pages accessible from the buttons on the home page"""
    def __init__(self,*args,**kwargs):
        Tk.__init__(self,*args,**kwargs)
        
        
        
        allframeshere= Frame(self)
        allframeshere.pack(side="top",fill="both",expand=True)
        allframeshere.grid_rowconfigure(0,weight=1)
        allframeshere.grid_columnconfigure(0,weight=1)
        self.frames = {}
        for A in (Homepage,FourierMain,PageOne,PageTwo,PageThree):
            Pageloc= A(allframeshere,self)
            self.frames[A]=Pageloc 
            Pageloc.grid(row=0,column=0,sticky="WENS")           
        self.viewframe(Homepage)
    def viewframe(self,s):
        Pageloc=self.frames[s]
        Pageloc.tkraise()

class Homepage(Frame):
    """Buttons to take the user from the home page to different parts of the GUI. Everything in this
    class contributes to what is displayed on the home page i.e. the page that first appears when
    the program is run."""
    def __init__(self, parent, controller):
        Frame.__init__(self,parent)
        self.controller=controller
        label=Label(self,text="Fourier Series Simulation Home Page",font=titlefont)
        label.pack(side="top",fill="x",pady=10)
        label = Label(self, text="""Welcome to the Fourier approximation simulation!""",font=("Times New Roman", 24, "bold"))
        label.pack(fill="x")
        label.pack(side="top",fill="x",pady=10)
        label = Label(self, text="""Use the buttons to navigate""",font=("Times New Roman", 20))
        label.pack(fill="x")
        b1=Button(self,text="Fourier Plots",font="titlefont",command=lambda:controller.viewframe(FourierMain))
        b1.pack()
        b2=Button(self,text="What does changing the time period do?",font="titlefont",command=lambda:controller.viewframe(PageOne))
        b2.pack()
        b3=Button(self,text="What does changing the number of harmonics do?",font="titlefont",command=lambda:controller.viewframe(PageTwo))
        b3.pack()
        b4=Button(self,text="Theory of Fourier series",font="titlefont",command=lambda:controller.viewframe(PageThree))
        b4.pack()
        b5=Button(self,text="PDF report",font="titlefont",command=lambda: os.startfile("final_project_report.pdf"))
        b5.pack()
        # Joseph Fourier himself greets the user
        frame=Frame(self) 
        frame.pack()
        self.placeImage(frame)
    def placeImage(self,frame): 
        pilImage = Image.open("joseph_fourier.png")
        self.photo=ImageTk.PhotoImage(pilImage)
        self.imagejoe = Label(frame,image=self.photo)
        self.imagejoe.pack()        
       
class FourierMain(Frame):
    """This is the class that defines the Fourier functions which are plotted on one page"""
    def __init__(self,parent,controller):
        Frame.__init__(self,parent)
        self.controller=controller
        label = Label(self, text="Fourier Plots",font=titlefont)
        label.pack(side="top",fill="x", pady=10)
        b=Button(self,text="Back",command=lambda: controller.viewframe(Homepage)) # Button to return user to homepage
        b.pack()
               
        frame=Frame(self)       
        frame.pack()
        self.makePlot(frame)
        self.makeInputs(frame)
    
    def makePlot(self, frame):
        self.fig = Figure(figsize=(6,5), dpi=100) # Makes a figure object
        self.figPlot = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig,frame)
        self.canvas.show() # Makes canvas visible
        self.canvas.get_tk_widget().grid(row = 0, column = 1)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, frame) # Toolbar
        self.toolbar.update()
        self.toolbar.grid(row = 7, column = 1)
        
        label = Label(self, text="""Use the sliders to choose values for the time period and number of harmonics for the Fourier series approximation of a square wave signal. 
        Press PLOT to display the result!""",font=INSTRUCTION_FONT) # Instructions that tell the user how to produce a plot
        label.pack(fill="x")

        
    
    def makeInputs(self, frame):
        InputFrame = Frame(frame)
        InputFrame.grid(column = 0, row = 0)
        
        #Arbitrary labels have to be included
        self.T = Label(InputFrame, text = "", justify=LEFT)
        self.T.grid(column= 0,row= 1)
        
        self.n_ = Label(InputFrame, text = "", justify=LEFT)
        self.n_.grid(column= 0,row= 2)
        
       
        #Sliders to produce input values. This makes the program more robust so values that are large enough
        #to break the program are inacessible.

        self.T = Scale(from_=10, to=40,length=400, tickinterval=5, orient=HORIZONTAL, command = self.T, label='Time period',font=contentfont,bg="purple", fg="yellow",borderwidth=5)
        self.T.pack(side="left")
        
        self.n_ = Scale(from_=0, to=50,length=400,tickinterval=5, orient=HORIZONTAL, command = self.n_, label='Number of harmonics',font=contentfont,bg="purple",fg="yellow",borderwidth=5) 
        self.n_.pack(side="right")
        
       
        # Button to plot the graph
        self.butPlot = Button(InputFrame, text ='PLOT', command = self.calcPattern, font=("Calibri",16,"bold"),foreground="red")# pressing the button runs the method calcPattern
        self.butPlot.grid(column=1,row=5)
        
        # Button to clear the axes
        self.butQuit = Button(InputFrame,  text="CLEAR",command=self.clearPlot, font=("Calibri",16,"bold"),foreground="blue")
        self.butQuit.grid(column=1,row=6)  
    
    def clearPlot(self): # This is the function that is called when the axes are to be cleared
        self.figPlot.cla()
        self.canvas.show()    
        
    def calcPattern(self):
        # x range and integer values for T and n
        x_ = np.linspace(-20,20,1000)
        T = int(self.T.get())
        n_ = int(self.n_.get())
        
       
        def Square(x): 
            bottom_left = (-T/2) #Defining positions of square wave vertices
            bottom_right = 0
            top_left = 0
            top_right = (T/2)
            # Two equations that construct a square wave:    
            maximum = 1 # f(x) = 1
            minimum = -1 # f(x) = -1
        
            while True: # Defining which values of x correspond to which intervals of T
                
                if (x >= bottom_left) and (x <= bottom_right):
                    return minimum 
                elif (x >= top_left) and (x <= top_right):
                    return maximum
        
                else:
                    bottom_left -= T/2
                    bottom_right -= T/2
                    top_left += T/2
                    top_right += T/2
        
                    if maximum == 1:
                        maximum = -1
                        minimum = 1
        
                    else:
                        maximum = 1
                        minimum = -1
        
        # Since this is an odd function, all a_n coefficients = 0 so don't need to be calculated
        
        def bsubn(n): # Calculating b_n coefficients
            n = int(n) # n must be a positive integer
            
            if (n%2 != 0):
                return 4/(np.pi*n)
        
            else:
                return 0
            
        def approximate(n,x): # Defining the Fourier approximation
            a0 = 0 # The average value of the set of functions is zero
            series = a0
        
            for n in range(1,n):
        
                try:
                    series = series + bsubn(n)*np.sin((2*np.pi*n*x)/T) # This is the expression for the Fourier wave approximation
        
                except:
                    pass
        
            return series
        
        y = []
        f = []
        
        for i in x_: # Calling previous functions to retrieve the relevant values
            y.append(Square(i)) 
            f.append(approximate(n_,i))
        
        # Title, legend and axes labels
        
        imPlot = self.figPlot.plot(x_,y,color='green',label='Square wave')
        imPlot = self.figPlot.plot(x_,f,color='orange',label='Fourier approximation')
        imPlot = self.figPlot.legend(loc=0)
        imPlot = self.figPlot.set_ylabel(r'$f(x)$')
        imPLot = self.figPlot.set_xlabel(r'$x$')
        imPlot = self.figPlot.set_title('Fourier series approximation of a square wave for n = '+str(n_)) 
        # The value of n is updated to the value chosen from the slider and put in the title of the plot
        
        self.canvas.show()  # Makes the plot appear
   
    
        

class PageOne(Frame):
    """Separate page displaying information on changing T values"""
    def __init__(self,parent,controller):
        Frame.__init__(self,parent)
        self.controller=controller
        label = Label(self, text="Changing the time period",font=titlefont)
        label.pack(side="top",fill="x", pady=10)
        b=Button(self,text="Back",command=lambda: controller.viewframe(Homepage)) #Button returns use to homepage
        b.pack()
       
        frame=Frame(self) 
        frame.pack()
        self.placeImage(frame)
        
    def placeImage(self,frame): # Imports image of the text 
        pilImage = Image.open("change_t_content.png")
        self.photo=ImageTk.PhotoImage(pilImage)
        self.imageCT = Label(frame,image=self.photo)
        self.imageCT.pack()
    
     
        
        
class PageTwo(Frame):
    """Separate page displaying information on changing n values"""
    def __init__(self,parent,controller):
        Frame.__init__(self,parent)
        self.controller=controller
        label = Label(self, text="Changing the number of harmonics",font=titlefont)
        label.pack(side="top",fill="x", pady=10)
        b=Button(self,text="Back",command=lambda: controller.viewframe(Homepage))
        b.pack()
       
        frame=Frame(self)       
        frame.pack()
        self.placeImage(frame)
    def placeImage(self,frame):# Imports image of the text
        pilImage = Image.open("change_n_content.png")
        self.photo=ImageTk.PhotoImage(pilImage)
        self.imageCn = Label(frame,image=self.photo)
        self.imageCn.pack()                                      

class PageThree(Frame):
    """Separate page displaying information on background Fourier theory"""
    def __init__(self,parent,controller):
        Frame.__init__(self,parent)
        self.controller=controller
        label = Label(self, text="Theory of Fourier series",font=titlefont)
        label.pack(side="top",fill="x", pady=10)
        b=Button(self,text="Back",command=lambda: controller.viewframe(Homepage))
        b.pack()
     
        frame=Frame(self)       
        frame.pack()
        self.placeImage(frame)
    def placeImage(self,frame):# Imports an image of LaTeX formatted equations and text
        pilImage = Image.open("fourier_theory_content.png")
        self.photo=ImageTk.PhotoImage(pilImage)
        self.imageFT = Label(frame,image=self.photo)
        self.imageFT.pack()

    

root = MainMenu()
root.title("Fourier Series GUI")
root.mainloop() #This produces the GUI


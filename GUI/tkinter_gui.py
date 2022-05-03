import tkinter as tk
from tkinter import ttk
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import pyplot as plt

# set some font
LARGE_FONT = ('Verdana', 14)

# class for creating main and other windows
class ContactFEM(tk.Tk):

    def __init__(self, graph, *args, **kwargs):
        self.graph = graph
        self.fig1, self.ax = plt.subplots(3, 2)
        self.fig2, self.ax2 = None, None
        self.i = 0  # number of printed step of Lemke's algorithm
        self.steps = graph.lemke.steps  # total number of steps to solve LCP

        self.plot_initial_visuals()  # plot schemes on fig1

        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.iconbitmap(self, default=os.path.realpath(__file__).split("GUI", 1)[0]+'contactFEM.ico')
        tk.Tk.wm_title(self, 'Contact FEM program')

        container = tk.Frame(self)  # initialize container
        container.pack(side='top', fill='both', expand=True)
        container.grid_rowconfigure(0, weight=1)  # 0 is minimum size, weight - priority in window
        container.grid_columnconfigure(0, weight=1)

        # Creating menu
        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label='save project', command=lambda: tk.messagebox.showinfo(message='Not supported just yet!'))
        menubar.add_cascade(label='File', menu=filemenu)
        tk.Tk.config(self, menu=menubar)

        # Creating multiple frames
        self.frames = {}
        for frame_obj in (StartPage, Page1, Page2, Page3):

            frame = frame_obj(container, self)
            self.frames[frame_obj] = frame

            frame.grid(row=0, column=0, sticky='nsew')  # 0 and 0 number of rows and col. nsew = North South East Wets
        # plot canvas
        self.mpl_canvas(container, self)  # plot

        self.show_frame(StartPage)
        # this is below everything (need to finish the process after closing window)
        self.protocol('WM_DELETE_WINDOW', self.close_app)

    def show_frame(self, controller):
        """
        Shows the cont frame
        :param controller: We choose from self.frames {} dictionary the frame which we need to show.
        :return:
        """
        frame = self.frames[controller]
        frame.tkraise()

    def mpl_canvas(self, parent, controller):
        self.canvas = FigureCanvasTkAgg(self.fig1, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, self)
        toolbar.update()
        self.canvas.mpl_connect('key_press_event', self.on_key_press)

    def on_key_press(self, event):
        print("you pressed {}".format(event.key))
        if event.key == 'left':
            self.left()
        if event.key == 'right':
            self.right()

    def right(self):
        if self.i == self.steps:
            return
        self.i += 1
        print('Step =', self.i)

        # new
        self.canvas.flush_events()  # do not know how it works and if i need this
        self.ax[0, 1].clear()
        self.ax[1, 1].clear(), self.ax[2, 1].clear()
        self.ax[2, 0].clear()
        self.graph.plot_scheme_i_step_lemke(plt_def_contact=self.ax[0, 1],
                                            plt_n=self.ax[1, 1], plt_t =self.ax[2, 1], i=self.i)
        self.ax[2, 0].text(0, 0, 'step:'+str(self.i), color='white')
        self.canvas.draw_idle()

    def left(self):
        if self.i == 0:
            return
        self.i -= 1
        print('Step =', self.i)
        self.ax[0, 1].clear()
        self.ax[1, 1].clear(), self.ax[2, 1].clear()
        self.ax[2, 0].clear()
        self.graph.plot_scheme_i_step_lemke(plt_def_contact=self.ax[0, 1],
                                            plt_n=self.ax[1, 1], plt_t =self.ax[2, 1], i=self.i)
        self.ax[2, 0].text(0, 0, 'step:' + str(self.i), color='white')
        self.canvas.draw_idle()

    # to end the process after closing window
    def close_app(self):
        exit()

    def set_parameters(self, fig, axes):
        """
        Setting parameters to figure
        :param fig:
        :param axes:
        :return:
        """
        fig.patch.set_facecolor('black')
        # set size
        fig.set_size_inches(13, 6, forward=True)
        for ax in axes.reshape(-1):
            #ax.set_aspect('equal')
            #ax.set_aspect('equal', adjustable='box')
            ax.set_facecolor('black')
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('white')
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            # ax.axhline(y=0, color='k')
            # ax.axvline(x=0, color='k')

    def plot_initial_visuals(self):
        # set parameters to axes
        self.set_parameters(fig=self.fig1, axes=self.ax)
        # plot design scheme
        self.graph.plot_supports(self.ax[0, 0], self.graph.sm.supports)
        self.graph.plot_nodes(self.ax[0, 0], deformed=False)
        self.graph.plot_list_of_elements(self.ax[0, 0], deformed=False)
        self.graph.plot_external_forces(self.ax[0, 0])
        # plot deformed linear scheme
        self.graph.plot_nodes(self.ax[1, 0], deformed=True)
        self.graph.plot_deformed_scheme(self.ax[1, 0], self.graph.u_linear_const)
        # plot deformed nonlinear scheme
        self.graph.plot_nodes(self.ax[0, 1], deformed=True, linear=False)
        self.graph.plot_def_i_step_lemke(plt_def_contact=self.ax[0, 1], i=self.steps)
        self.graph.plot_lcp_nt_i_step_lemke(plt_n=self.ax[1, 1], plt_t=self.ax[2, 1], i=self.steps)


class StartPage(tk.Frame):

    def __init__(self, parent, controller, *args, **kwargs):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text='Start Page', font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text='Visit Page 1', command=lambda: controller.show_frame(Page1))
        button1.pack()
        button2 = ttk.Button(self, text='Visit Page 2', command=lambda: controller.show_frame(Page2))
        button2.pack()
        button3 = ttk.Button(self, text='Visit Page 3', command=lambda: controller.show_frame(Page3))
        button3.pack()

class Page1(tk.Frame):
    """
    frame = frame_obj(container, self) - that's how class is crated
    :param container:  = tk.Frame(tk.Tk)
    :param self: =  tk.Tk
    in the end we get
    famer = frame_obj(tk.Frame(tk.Tk), tk.Tk, ax)
    """
    def __init__(self, parent, controller, *args, **kwargs):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text='Page One', font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text='Back to Start page', command=lambda: controller.show_frame(StartPage))
        button1.pack()
        button2 = ttk.Button(self, text='Visit page 2', command=lambda: controller.show_frame(Page2))
        button2.pack()


class Page2(tk.Frame):

    def __init__(self, parent, controller, *args, **kwargs):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text='Page Two', font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text='Back to Start page', command=lambda: controller.show_frame(StartPage))
        button1.pack()
        button2 = ttk.Button(self, text='To page 1', command=lambda: controller.show_frame(Page1))
        button2.pack()

        if controller.fig2 is not None:
            canvas = FigureCanvasTkAgg(controller.fig2, self)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            toolbar = NavigationToolbar2Tk(canvas, self)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class Page3(tk.Frame):
    """
    frame = frame_obj(container, self) - that's how class is crated
    :param container:  = tk.Frame(tk.Tk)
    :param self: =  tk.Tk
    in the end we get
    famer = frame_obj(tk.Frame(tk.Tk), tk.Tk, ax)
    """
    def __init__(self, parent, controller,  *args, **kwargs):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text='Page Three', font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text='Back to Start page', command=lambda: controller.show_frame(StartPage))
        button1.pack()

        self.ax = controller.ax
        self.fig1 = controller.fig1

        #self.mpl_canvas(parent, controller)

    def mpl_canvas(self, parent, controller):
        self.canvas = FigureCanvasTkAgg(self.fig1, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, self)
        toolbar.update()
        self.canvas.mpl_connect('key_press_event', self.on_key_press)

    def on_key_press(self, event):
        print("you pressed {}".format(event.key))
        if event.key == 'left':
            print('left')
        if event.key == 'right':
            print('right')




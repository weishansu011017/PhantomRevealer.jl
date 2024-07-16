import numpy as np
import matplotlib.figure as mfg
import matplotlib.axes as maxe
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gspec

#Include LaTeX rendering
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 13,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsfonts}"
})

# Set nan if divided by zero
np.seterr(divide='ignore')

def openinteractive():
    '''
    Open the interactive interface for plotting
    '''
    if not plt.isinteractive():
        plt.ion()
            
def Get_vminmax(array:np.ndarray,minzero:bool = False):
    '''
    Calculate the minimum and maximum value for colorbar.
    '''
    median = np.nanmedian(array)
    std = np.nanstd(array)
    vmax = median + 3*std
    vmin = median - 3*std
    if minzero and vmin<0.0:
        vmin = 1e-15
        
    print(f"Warning: Automatically calculate (vmin, vmax) = {(vmin,vmax)}")
    return np.array([vmin,vmax])
        
def polar2cart(theta:float, radius:float):
    '''
    Transfer the polar coordinate into the cartisian coordinate.
    '''
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    return x, y

def Basic_properties_array(arr,array_name):
    max = np.nanmax(arr)
    min = np.nanmin(arr)
    mean = np.nanmean(arr)
    median = np.nanmedian(arr)
    std = np.nanstd(arr)
    print(f'----------------Properties of Array {array_name}----------------')
    print(f'Maximum: {max}')
    print(f'Minimum: {min}')
    print(f'Mean: {mean}')
    print(f'Median: {median}')
    print(f'STD: {std}')
    print('----------------------------------------------------------------')

class figure_ax:
    '''
    The plotting object in the code.
    
    Field in Class:
        fig: The figure/Canvas.
        ax: The axes of plotting/colorbar
        proj: The projection of plotting
        
    Usage:
        1. Generate an object
            fax = figure_ax()
            
        2. Setup the desired figure
            fax.setup_fig(ncols=2, nrows=2, figsize=(10,6))
            
        3. Make the plot at each axes
            fax.ax[0].plot(x,y,*plotting_params)
            cont = fax.ax[1].pcolor(x,y,z,*plotting_params)
            
        4. Setup the colorbar for each map (if necessary). Allowing giving the colorbar a specific label so that the new colorbar can be drawn on the old ax rather then the new one.
            colorbar = fax.setup_colorbar(cont,'LABEL_OF_COLORBAR')
            
        5. Save the figure
            fax.fig.savefig(filepath, dpi=dpi,bbox_inches='tight',transparent=False)
        
        6. Close the figure
            fax.close_fig()

        
    '''
    def __init__(self,proj=None):
        '''
        Two kinds of projection
        "None": Normal cartisian plotting
        "polar": Polar plotting
        '''
        self.fig: mfg.Figure = None
        self.ax: maxe._axes.Axes = None
        self.ncols = None
        self.nrows = None
        self.proj = proj
        
    def checking_initialized(self):
        if self.fig is None or self.ax is None:
            raise NameError(f"No figure or ax has found! You should initialized the figure by the 'setup_figure()' method.")
        
    def setup_fig(self, nrows=1,ncols=1,figsize = (10, 6)):
        openinteractive()
        if self.fig is None or self.ax is None:
            self.ncols = ncols
            self.nrows = nrows
            self.fig = plt.figure(figsize=figsize, constrained_layout=True)
            gs = gspec.GridSpec(nrows, ncols, figure=self.fig, hspace=0.0, wspace=0.0)

            ax_list = []

            for c in range(ncols):
                for r in range(nrows):
                    if ncols > 1 and nrows > 1:
                        sharex = 'all'
                        sharey = 'all'
                    elif ncols > 1:
                        sharex = None
                        sharey = 'all'
                    elif nrows > 1:
                        sharex = 'all'
                        sharey = None
                    else:
                        sharex = None
                        sharey = None

                    if c == 0 and r == 0:
                        ax = self.fig.add_subplot(gs[r, c], projection=self.proj)
                    else:
                        ax = self.fig.add_subplot(gs[r, c], projection=self.proj, sharex=ax_list[0] if sharex else None, sharey=ax_list[0] if sharey else None)
                    ax_list.append(ax)
                    
                    if self.proj != 'polar':
                        if c != 0:  
                            plt.setp(ax.get_yticklabels(), visible=False)
                        
                        if r != nrows - 1:  
                            plt.setp(ax.get_xticklabels(), visible=False)
                    else:
                        ax.set_xticklabels([]) 
                        ax.set_yticklabels([])  
                        ax.xaxis.set_visible(False)  
                        ax.yaxis.set_visible(False)  
                        

            self.ax = ax_list[0] if (ncols == 1 and nrows == 1) else ax_list
        
            if self.proj != 'polar':
                plt.tight_layout()
                
    def get_number_of_ax(self):
        ax = self.ax
        if type(ax) != list:
            return 1
        else:
            return len(self.ax)
            
    def setup_colorbar(self,cont,colorbar_ax_label = None):
        ''' Setup the colobar for the specific heatmap.
        
        The heatmap can be selected by specifing the "colorbar_ax_label" variable.
        '''
        if colorbar_ax_label is None:
            colarbar_ax = None
            for ax in self.fig.axes:
                if ax.get_label() == '<colorbar>':
                    colarbar_ax = ax
            if colarbar_ax is None:
                colorbar = self.fig.colorbar(cont)
            else:
                colorbar = self.fig.colorbar(cont,cax = colarbar_ax)
        else:
            colarbar_ax = None
            for ax in self.fig.axes:
                if ax.get_label() == colorbar_ax_label:
                    colarbar_ax = ax
            if colarbar_ax is None:
                colorbar = self.fig.colorbar(cont)
                colorbar.ax.set_label(colorbar_ax_label)
            else:
                colorbar = self.fig.colorbar(cont,cax = colarbar_ax)
        return colorbar
    
    def draw_fig(self):
        if not self.fig is None:
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
    
    def save_fig(self,filepath='figure_1.png',dpi=450):
        if self.fig:
            self.fig.savefig(filepath, dpi=dpi,bbox_inches='tight',transparent=False)
        else:
            print("No figure to save.")
    
    def close_fig(self):
        if self.fig is not None and self.ax is not None:
            plt.close(self.fig) 
            self.fig = None
            self.ax = None
            
    def reset_fig(self):
        ncols = self.ncols
        nrows = self.nrows
        figsize = tuple(self.fig.get_size_inches())
        self.close_fig()
        self.setup_fig(nrows,ncols,figsize)
            
class two_axes_plot(figure_ax):
    props = dict(boxstyle='round', facecolor='black')
    anato_text_position=(-0.02,0.95)
    def __init__(self, x:np.ndarray, y:np.ndarray, xlabel:str, ylabel:str, proj=None):
        super().__init__(proj)
        self._x : np.ndarray = np.array(x)
        self._y : np.ndarray = np.array(y)
        self._xlabel : str = xlabel
        self._ylabel : str = ylabel
        self.meshgrid()
    
    def grid_shape(self):
        return (self._y.shape[0], self._x.shape[0])
    
    def meshgrid(self):
        self._X,self._Y = np.meshgrid(self._x, self._y)
        
    def get_colorbar_pos_indies(self):
        indies = np.zeros(self.nrows)
        for i in range(0,self.nrows):
            indies[i] = (i+1)*self.ncols-1
        return indies
    
    def pcolor_draw(self,images:list,colormaps:list, clabels:list, Log_flags:list, anatonate_labels:list, minzero_flags = None, vlims = None,draw=True):
        '''
        The clims parameters should be the following format:
        
        [[vmin1,vmax1], [vmin2,vmax2], [vmin3,vmax3], ... ]
        
        '''
        cls = self.__class__
        # Checking whether the number of figures is matching.
        self.checking_initialized()
        if len(images) != self.get_number_of_ax():
            raise ValueError(f"Mismatching number of arrays! The number of figures was set to be {self.get_number_of_ax()}, but {len(images)} elements in images was given.")
        if len(colormaps) != self.nrows:
            raise ValueError(f"Mismatching number of colormaps! The number of rows of figures was set to be {self.nrows}, but {len(colormaps)} elements in colormaps was given.")
        if len(Log_flags) != self.nrows:
            raise ValueError(f"Mismatching number of Log_flags! The number of rows of figures was set to be {self.nrows}, but {len(Log_flags)} elements in Log_flags was given.")
        if len(clabels) != self.nrows:
            raise ValueError(f"Mismatching number of clabels! The number of rows of figures was set to be {self.nrows}, but {len(clabels)} elements in clabels was given.")
        if len(anatonate_labels) != self.ncols:
            raise ValueError(f"Mismatching number of anatonate_labels! The number of columns of figures was set to be {self.ncols}, but {len(anatonate_labels)} elements in anatonate_labels was given.")
        if (not (minzero_flags is None)) and (len(minzero_flags) != self.nrows):
            raise ValueError(f"Mismatching number of minzero_flags! The number of rows of figures was set to be {self.nrows}, but {len(minzero_flags)} elements in minzero_flags was given.")
        if (minzero_flags is None):
            minzero_flags = np.zeros(self.nrows)
        
        # Preparing Color limit
        colorbar_pos_indies = self.get_colorbar_pos_indies()
        if vlims is None:
            vlims = np.empty(self.nrows, dtype=np.ndarray)
            j = 0
            for i in range(self.get_number_of_ax()):
                if i in colorbar_pos_indies:
                    image = images[i]
                    vlims[j] = Get_vminmax(image, minzero_flags[j])
                    j += 1
        else:
            if len(vlims) != self.nrows:
                raise ValueError(f"Mismatching number of colormap limits! The number of rows of figures was set to be {self.nrows}, but {len(vlims)} pairs of colormap limit was given.")
            else:
                for vlim in vlims:
                    if len(vlim) != 2:
                        raise ValueError(f"Wrong format of 'vlims' parameters! The clims parameters should be [[vmin1,vmax1], [vmin2,vmax2], [vmin3,vmax3], ... ]")
            
        # Preparing Color Norm
        Norms = np.empty(self.nrows, dtype=object)
        Norms[:] = [mcolors.LogNorm(*vlims[j]) if Log_flags[j] else mcolors.Normalize(*vlims[j]) for j in range(len(Log_flags))]
        
        # Preparing interior colormap ax labels
        cax_label = [f"<colorbar_col{n+1}>" for n in range(self.nrows)]
        # Preparing grid
        X = self._X
        Y = self._Y
        ax = self.ax
        
        # Checking the grid size is mathing, and then drawing.
        j = 0
        for i,image in enumerate(images):
            image = np.array(image)
            k = i % self.ncols
            if image.shape != self.grid_shape():
                raise ValueError(f"Mismatching shape of array! The grid size of figures was set to be ({self.grid_shape()}), but {image.shape} was given.")
            else:
                if type(ax) == list:
                    cont = ax[i].pcolor(X,Y,image, cmap=colormaps[j], norm=Norms[j])
                    ax[i].text(*cls.anato_text_position,anatonate_labels[k] ,transform=ax[i].transAxes,c='white', fontsize=14, verticalalignment='top', bbox=cls.props)
                    if i in colorbar_pos_indies:
                        colorbar = self.setup_colorbar(cont,colorbar_ax_label=cax_label[j])
                        colorbar.set_label(clabels[j])
                        j += 1
                else:
                    cont = ax.pcolor(X,Y,image, cmap=colormaps[0], norm=Norms[0])
                    ax.text(*cls.anato_text_position,anatonate_labels[0] ,transform=ax.transAxes,c='white', fontsize=14, verticalalignment='top', bbox=cls.props)
                    colorbar = self.setup_colorbar(cont,colorbar_ax_label=cax_label[0])
                    colorbar.set_label(clabels[0])
                    
        if draw:
            self.draw_fig()
            
    def streamplot_draw(self,vxgrids:list,vzgrids:list,color="black",density=2.0,draw=True):
        self.checking_initialized()
        if len(vxgrids) != self.get_number_of_ax():
            raise ValueError(f"Mismatching number of arrays! The number of figures was set to be {self.get_number_of_ax()}, but {len(vxgrids)} elements in images was given.")
        if len(vzgrids) != self.get_number_of_ax():
            raise ValueError(f"Mismatching number of arrays! The number of figures was set to be {self.get_number_of_ax()}, but {len(vzgrids)} elements in images was given.")
        X = self._X
        Y = self._Y
        ax = self.ax
        
        for i in range(self.get_number_of_ax()):
            vxgrid = np.array(vxgrids[i])
            vzgrid = np.array(vzgrids[i])
            
            if vxgrid.shape != self.grid_shape():
                raise ValueError(f"Mismatching shape of array! The grid size of figures was set to be ({self.grid_shape()}), but {vxgrid.shape} was given.")
            elif vzgrid.shape != self.grid_shape():
                raise ValueError(f"Mismatching shape of array! The grid size of figures was set to be ({self.grid_shape()}), but {vzgrid.shape} was given.")
            else:
                if type(ax) == list:
                    ax[i].streamplot(X,Y,vxgrid,vzgrid,color=color, density=density)
                else:
                    ax.streamplot(X,Y,vxgrid,vzgrid,color=color, density=density)
        if draw:
            self.draw_fig()
        
class polar_plot(two_axes_plot):
    anato_text_position=(-0.02,0.95)
    def __init__(self, s: np.ndarray, phi: np.ndarray, slabel: str, philabel: str):
        super().__init__(phi, s, philabel, slabel,  'polar')

    @property
    def s(self):
        return self._y
    
    @property
    def phi(self):
        return self._x
    
    @property
    def slabel(self):
        return self._ylabel
    
    @property
    def philabel(self):
        return self._xlabel
    
    @property
    def S(self):
        return self._Y
    
    @property
    def PHI(self):
        return self._X
            
class cart_plot(two_axes_plot):
    anato_text_position=(0.0065,0.98)
    def __init__(self, x: np.ndarray, y: np.ndarray, xlabel: str, ylabel: str):
        super().__init__(x, y, xlabel, ylabel)
    
    @property
    def y(self):
        return self._y
    
    @property
    def x(self):
        return self._x
    
    @property
    def ylabel(self):
        return self._ylabel
    
    @property
    def xlabel(self):
        return self._xlabel
    
    @property
    def Y(self):
        return self._Y
    
    @property
    def X(self):
        return self._X
    
    def set_ylabel(self,axid):
        self.ax[axid].set_ylabel(self.ylabel)
    
    def set_xlabel(self,axid):
        self.ax[axid].set_xlabel(self.xlabel)